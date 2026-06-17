source("R/packages.R")
source("R/eval_utils.R")
source("R/fda_utils.R")

y <- readRDS("datasets/y_train.rds")
y_test <- readRDS("datasets/y_test.rds")
folds <- readRDS("datasets/folds.rds")
raw_train_by_layer <- readRDS("datasets/raw_train_by_layer.rds")

df_raw_test <- readRDS("datasets/df_raw_test.rds")
raw_test_by_layer <- list()
for (l in unique(df_raw_test$names)) {
  raw_test_by_layer[[l]] = df_raw_test %>% filter(names == l)
  raw_test_by_layer[[l]] = as.data.frame(raw_test_by_layer[[l]])
}

hc <- readRDS("results/class_ord_from_healthy.rds")   

basis <- create.bspline.basis(c(0,1), 60, norder=4)

lambda_smooth <- 1e-04   # pooled-GCV optimal, fixed before nested CV


# ------------------------------------------------------------------------------
# Utils for the model

fit_mfpca <- function(fd_layers, M = 20, k_uniexp = 30) {
  # fd_layers: list of 6 funData objects (one per layer), as produced by
  #            smooth_all_layers() or build_funData_layer()
  # M:         number of multivariate PCs to compute
  # k_uniexp:  dimension of the basis used in the per-layer univariate FPCA
  #            (should not exceed the smoothing basis dimension)
  
  # Combine the 6 univariate funData objects into a multiFunData
  mfd <- funData::multiFunData(fd_layers)
  
  # One uniExpansions entry per layer 
  # Splines1D matches the smoothing basis type and is efficient
  uni_exp <- lapply(seq_along(fd_layers), function(i) {
    list(type = "splines1D", k = k_uniexp)
  })
  
  MFPCA::MFPCA(
    mFData        = mfd,
    M             = M,
    uniExpansions = uni_exp,
    fit           = TRUE 
  )
}

project_onto_mfpca <- function(fd_layers_new, mfpca_fit) {
  # fd_layers_new: list of 6 funData objects for new patients,
  #                evaluated on the SAME grid as the training data
  # mfpca_fit:     MFPCA object from fit_mfpca() on training data
  #
  # Returns: n_new × M score matrix where row i is patient i's
  # projection onto the training eigenfunctions psi_1, ..., psi_M
  
  M <- length(mfpca_fit$values)             # number of components
  t_grid  <- fd_layers_new[[1]]@argvals[[1]]       # shared grid across layers
  n_new <- nrow(fd_layers_new[[1]]@X)            # number of new patients
  n_lyr <- length(fd_layers_new)
  
  # Numerical integration weights for the trapezoidal rule on t_grid
  # trapezoidal weights: half-step at the endpoints, full step in the middle
  dt      <- diff(t_grid)
  weights <- c(dt[1], dt[-1] + dt[-length(dt)], dt[length(dt)]) / 2
  
  # Population mean function per layer (from training)
  # Used to center new curves before projecting (eigenfunctions are
  # defined relative to the training mean, not zero)
  mu <- lapply(seq_len(n_lyr), function(l) {
    mfpca_fit$meanFunction[[l]]@X[1, ]   # 1 × n_grid 
  })
  
  # Compute scores
  # For each new patient i and component k, sum across layers l of the
  # weighted inner product between centered curve and k-th eigenfunction.
  scores <- matrix(0, nrow = n_new, ncol = M)
  
  for (k in 1:M) {
    for (l in 1:n_lyr) {
      # Layer-l component of the k-th multivariate eigenfunction
      psi_kl <- mfpca_fit$functions[[l]]@X[k, ]   # length n_grid
      
      # Centered curves: n_new × n_grid
      X_centered <- sweep(fd_layers_new[[l]]@X, 2, mu[[l]], FUN = "-")
      
      # Weighted inner product per patient:
      # X_centered %*% (psi_kl * weights) computes the integral for all patients at once
      scores[, k] <- scores[, k] + as.numeric(X_centered %*% (psi_kl * weights))
    }
  }
  
  scores
}

align_signs <- function(mfpca_fold, mfpca_ref) {
  n_lyr   <- length(mfpca_fold$functions)
  M_local <- nrow(mfpca_fold$functions[[1]]@X)
  t_grid  <- mfpca_fold$functions[[1]]@argvals[[1]]
  dt      <- diff(t_grid)
  w       <- c(dt[1], dt[-1] + dt[-length(dt)], dt[length(dt)]) / 2
  
  for (k in 1:M_local) {
    ip <- sum(sapply(1:n_lyr, function(l) {
      sum(mfpca_fold$functions[[l]]@X[k, ] *
            mfpca_ref$functions[[l]]@X[k, ] * w)
    }))
    if (ip < 0) {
      for (l in 1:n_lyr) {
        mfpca_fold$functions[[l]]@X[k, ] <- -mfpca_fold$functions[[l]]@X[k, ]
      }
      mfpca_fold$scores[, k] <- -mfpca_fold$scores[, k]
    }
  }
  mfpca_fold
}

smooth_all_layers <- function(raw_bylyr, lambda, basis) {
  lapply(seq_along(raw_bylyr), function(l) {
    results <- lapply(seq_len(nrow(raw_bylyr[[l]])), function(i) {
      smooth_patient_single_lam(raw_bylyr[[l]][i, , drop = FALSE],
                                k = NULL, lam = lambda, basis = basis)
    })
    t_grid <- seq(0, 1, length.out = 750)
    mat    <- t(sapply(results, function(obs) {
      as.numeric(eval.fd(t_grid, obs$smooth_res))
    }))
    funData::funData(argvals = t_grid, X = mat)
  })
}

# ------------------------------------------------------------------------------
# (Unpenalized) multinomial logreg via nnet::multinom on MFPCA scores
# Performs single-level CV: MFPCA is refit on each training fold
# validation patients are projected onto the training eigenfunctions 
# (sign-aligned to the reference fit on full training data)
# predictions are stored for evaluation
# do MFPCA's variance-optimal directions carry the
# disease discriminative information? 

fit_logreg_mfpca <- function(x_raw_bylyr, y, folds, basis, K,
                             M = 20,
                             lambda_smooth_fixed = 1e-04,
                             reference_class = "Healthy",
                             verbose = TRUE) {
  
  if (reference_class %in% levels(y)) y <- relevel(y, ref = reference_class)
  
  n     <- length(y)
  n_k   <- length(folds)
  testy <- y
  predy <- factor(rep(NA, n), levels = levels(y))
  
  mn_log_loss <- function(true_labels, prob_matrix) {
    prob_matrix <- pmax(prob_matrix, 1e-15)
    idx <- cbind(seq_along(true_labels), as.integer(true_labels))
    -mean(log(prob_matrix[idx]))
  }
  
  # Reference MFPCA: anchor for sign alignment across folds
  if (verbose) message("Fitting reference MFPCA on full training set...")
  fd_train_full     <- smooth_all_layers(x_raw_bylyr, lambda_smooth_fixed, basis)
  mfpca_ref         <- fit_mfpca(fd_train_full, M = M)
  train_scores_full <- mfpca_ref$scores[, 1:K, drop = FALSE]
  
  if (verbose) {
    cat(sprintf("Variance explained by first %d MFPCs: %.2f%%\n\n",
                K, sum(mfpca_ref$values[1:K]) / sum(mfpca_ref$values) * 100))
  }
  
  fold_logloss <- numeric(n_k)
  
  for (k_idx in seq_along(folds)) {
    val_idx   <- folds[[k_idx]]
    train_raw <- lapply(x_raw_bylyr, function(df) df[-val_idx, , drop = FALSE])
    val_raw   <- lapply(x_raw_bylyr, function(df) df[ val_idx, , drop = FALSE])
    y_train   <- y[-val_idx]
    y_val     <- y[ val_idx]
    
    fd_train <- smooth_all_layers(train_raw, lambda_smooth_fixed, basis)
    fd_val   <- smooth_all_layers(val_raw,   lambda_smooth_fixed, basis)
    mfpca_k  <- fit_mfpca(fd_train, M = M)
    mfpca_k  <- align_signs(mfpca_k, mfpca_ref)
    
    train_scores <- mfpca_k$scores[, 1:K, drop = FALSE]
    val_scores   <- project_onto_mfpca(fd_val, mfpca_k)[, 1:K, drop = FALSE]
    
    colnames(train_scores) <- paste0("MFPC", 1:K)
    colnames(val_scores)   <- paste0("MFPC", 1:K)
    
    df_train <- data.frame(label = y_train, train_scores)
    df_val   <- data.frame(val_scores)
    
    fit_k <- nnet::multinom(label ~ ., data = df_train, trace = FALSE, maxit = 500)
    
    probs_val <- predict(fit_k, newdata = df_val, type = "probs")
    preds_val <- predict(fit_k, newdata = df_val, type = "class")
    predy[val_idx] <- factor(as.character(preds_val), levels = levels(y))
    
    fold_logloss[k_idx] <- mn_log_loss(y_val, probs_val)
    
    if (verbose)
      cat(sprintf("Fold %d: mnLogLoss = %.4f\n", k_idx, fold_logloss[k_idx]))
  }
  
  if (verbose)
    cat(sprintf("\nMean CV mnLogLoss = %.4f (sd = %.4f)\n",
                mean(fold_logloss), sd(fold_logloss)))
  
  # Final model on full training set
  colnames(train_scores_full) <- paste0("MFPC", 1:K)
  df_full   <- data.frame(label = y, train_scores_full)
  final_fit <- nnet::multinom(label ~ ., data = df_full, trace = FALSE, maxit = 500)
  
  list(
    output       = list(testy = testy, predy = predy),
    final_fit    = final_fit,
    mfpca_full   = mfpca_ref,
    K            = K,
    fold_logloss = fold_logloss
  )
}


# Main fit is with K = 8 
# used for the coefficient plot and test-set evaluation
res_logreg <- fit_logreg_mfpca(raw_train_by_layer, y, folds, basis,
                               K = 8, M = 20,
                               lambda_smooth_fixed = lambda_smooth)

eval_logreg <- nestcv_modeval(res_logreg)
print(eval_logreg$metrics)

# final_mods <- readRDS("results/final_models.rds")
# final_mods[["logreg_mfpc_scores"]] <- res_logreg
# saveRDS(final_mods, file = "results/final_models.rds")
# saveRDS(res_logreg, file = "results/cv_logreg_mfpca.rds")


# ------------------------------------------------------------------------------
# coefficient heatmap
# nnet::multinom returns coefficients in reference-class parameterisation.
# Healthy is the reference (implicit zero coefficients0

coef_matrix <- coef(res_logreg$final_fit)
coef_matrix <- coef_matrix[, colnames(coef_matrix) != "(Intercept)", drop = FALSE]

ref_row <- matrix(0, nrow = 1, ncol = ncol(coef_matrix),
                  dimnames = list("Healthy", colnames(coef_matrix)))
coef_matrix <- rbind(ref_row, coef_matrix)

# Reorder rows to match the AIRM dendrogram
# coef_matrix <- coef_matrix[hc, , drop = FALSE]

pheatmap(
  coef_matrix,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  scale           = "none",
  color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks          = seq(-max(abs(coef_matrix)), max(abs(coef_matrix)), length.out = 101),
  border_color    = "grey85",
  cellwidth       = 44, cellheight = 38,
  fontsize        = 12, fontsize_row = 13, fontsize_col = 11,
  angle_col       = 45,
  display_numbers = TRUE,
  number_format   = "%.2f",
  number_color    = "grey20",
  main            = "Multinomial logreg coefficients on MFPCA scores (vs Healthy)",
  # treeheight_row  = 90,
  margins         = c(80, 10)
)


# ------------------------------------------------------------------------------
# log - odds plots
# Build coefficient long table from the multinom fit 
coef_mat <- coef(res_logreg$final_fit)   # (K_classes - 1) × (K + 1)
coef_mat <- coef_mat[, colnames(coef_mat) != "(Intercept)", drop = FALSE]

# Add Healthy as the zero reference row
coef_mat <- rbind(Healthy = rep(0, ncol(coef_mat)), coef_mat)

coef_long <- as.data.frame(coef_mat) %>%
  tibble::rownames_to_column("class") %>%
  tidyr::pivot_longer(-class, names_to = "feature", values_to = "log_odds") %>%
  dplyr::filter(class != "Healthy")   # drop the reference row from the bar plot

coef_long$feature <- factor(coef_long$feature,
                            levels = paste0("MFPC", res_logreg$K:1))

disease_palette <- c(
  "MH"     = "#F28E2B", 
  "DR"     = "mediumpurple", 
  "CSR"    = "#E15759", 
  "AMD"    = "seagreen3"
)

ggplot(coef_long, aes(x = log_odds, y = feature, fill = class)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ class, nrow = 1) +
  scale_fill_manual(values = disease_palette, guide = "none") +
  labs(
    title    = "Log-odds versus Healthy on MFPCA scores",
    subtitle = "Each bar: change in log-odds (disease vs. Healthy) per unit increase in the MFPC",
    x        = "Log-odds vs. Healthy (per unit MFPC)",
    y        = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text       = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y      = element_text(face = "bold")
  )


# ------------------------------------------------------------------------------
# test set evaluation
# Project the held-out test patients onto the training-set MFPCA eigenfunctions
# (inside res_logreg$mfpca_full), then predict with the final model

fd_test <- smooth_all_layers(raw_test_by_layer, lambda_smooth, basis)
test_scores <- project_onto_mfpca(fd_test, res_logreg$mfpca_full)[, 1:res_logreg$K, drop = FALSE]
colnames(test_scores) <- paste0("MFPC", 1:res_logreg$K)
df_test <- as.data.frame(test_scores)

preds_test <- predict(res_logreg$final_fit, newdata = df_test, type = "class")
preds_test <- factor(as.character(preds_test), levels = levels(y_test))

cm_test <- caret::confusionMatrix(preds_test, y_test)
cat("\n── Test-set evaluation ──\n")
print(cm_test)

test_metrics <- list(
  macro_F1     = mean(cm_test$byClass[, "F1"],                na.rm = TRUE),
  balanced_acc = mean(cm_test$byClass[, "Balanced Accuracy"], na.rm = TRUE),
  per_class_F1 = cm_test$byClass[, "F1"],
  confusion    = cm_test$table
)
cat(sprintf("\nTest macro_F1     = %.3f\n", test_metrics$macro_F1))
cat(sprintf("Test balanced_acc = %.3f\n", test_metrics$balanced_acc))

saveRDS(test_metrics, "results/test_metrics_logreg_mfpca.rds")
