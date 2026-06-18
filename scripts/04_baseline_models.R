source("R/packages.R", echo=FALSE)
source("R/model_defs.R")
source("R/eval_utils.R")

datasets = readRDS("datasets/train_datasets.rds")

y = readRDS("datasets/y_train.rds")

# the k (external) folds need to be common between models
set.seed(2026)
folds = createFolds(y, k = 10, returnTrain = FALSE)
# save them since they will need to be used also by other models
saveRDS(folds, "datasets/folds.rds")

folds = readRDS("datasets/folds.rds")

# if one wants to add a model (from library(nestedcv)), he just has to write the fit_* function inside R/model_defs.R
# and add the model to this list
models <- list(
  logreg = fit_logreg,
  rf = fit_rf,
  svm = fit_svm,
  knn = fit_knn
)

# list where to save the results of each model fitted on a specific dataset
results = list()

# actual loop for each model on each dataset
for (ds_name in names(datasets)) {
  for (mod_name in names(models)) {
    key <- paste(mod_name, ds_name, sep = "_") # key to lookup for results, ex: logreg_means
    message("Fitting ", key, " ...")
    set.seed(2026) # reproducibility
    fit <- models[[mod_name]](datasets[[ds_name]], y, folds) # models[[mod_name]] selects the correct fit function
    results[[key]] <- list(dataset = ds_name,
                           model   = mod_name,
                           fit = fit, 
                           eval = nestcv_modeval(fit, folds))
  }
}

saveRDS(results, "results/baseline_results.rds")

results = readRDS("results/baseline_results.rds")

# summary of the metrics used for the comparison
# rbind the results of the second argument 
summary_cv_df = summary_cv(results)

# order by F1 score
print(summary_cv_df |>
        dplyr::arrange(dplyr::desc(macro_F1), dplyr::desc(balanced_acc)))

# check if order is consistent with more features added
print(summary_cv_df |>
        dplyr::group_by(dataset) |>
        dplyr::arrange(dplyr::desc(macro_F1), .by_group = TRUE))

# F1 score per class 
# (order of classes is based on distances of the covariance matrices from class healthy)
perclassF1_df = perclass_F1_cv(results, class_ord_path = "results/class_ord_from_healthy.rds")

print(perclassF1_df)

summary_tbl[which.max(summary_tbl$macro_F1), ]


# ------------------------------------------------------------------------------
# Analysis of penalized multinomial logistic regression on means_logCV_ratios dataset
logreg_fit<- results[["logreg_means_logCV_ratios"]]$fit   
# note that glmnet fits the multinomial logistic regression in the "symmetric form"
# you do not pick a baseline class 
# you get estimates of coefficients for each class 

# heatmap of coefficents of the penalized logistic regression
all_features <- colnames(datasets[["means_logCV_ratios"]])

coef_matrix <- t(sapply(coef(logreg_fit), function(ck) {
  ck <- ck[names(ck) != "(Intercept)"]   # drop intercept
  out <- setNames(rep(0, length(all_features)), all_features)
  out[names(ck)] <- ck                   # fill in the non-zero values
  out
}))

# standardize coefficients
feature_sds <- apply(datasets$means_logCV_ratios, 2, sd)
sd_vec <- feature_sds[colnames(coef_matrix)]
coef_matrix_std <- sweep(coef_matrix, MARGIN = 2, STATS = sd_vec, FUN = "*")

# make healthy the reference
healthy_row <- coef_matrix_std["Healthy", ]
coef_matrix_std_ref <- sweep(coef_matrix_std, MARGIN = 2, STATS = healthy_row, FUN = "-")
coef_matrix_std_ref <- coef_matrix_std_ref[rownames(coef_matrix_std_ref) != "Healthy", ]

coefs_logreg <- pheatmap(
  coef_matrix_std_ref,
  # classes are clustered according to the correlation matrix distances
  cluster_rows    = F, # hc if you want also to plot Healthy 
  cluster_cols    = FALSE,
  scale           = "none",
  color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks          = seq(-max(abs(coef_matrix_std_ref)), max(abs(coef_matrix_std_ref)), length.out = 101),
  border_color    = "grey85",
  cellwidth       = 44,    
  cellheight      = 38,
  fontsize        = 12,
  fontsize_row    = 13,
  fontsize_col    = 11,
  angle_col       = 45,
  display_numbers = TRUE,
  number_format   = "%.2f",
  number_color    = "grey20",
  main            = "Standardized coefficients vs. Healthy\n(log-odds change per 1 SD increase in feature)",
  treeheight_row  = 90,
  margins         = c(80, 10) 
) 
ggsave("figures/classification/coefs_logreg.png", plot = coefs_logreg, width = 12, height = 8, 
       units = "in", dpi = 300)

# confusion matrix
cf_mat_logreg <- plot_confusion_matrix(results[["logreg_means_logCV_ratios"]])
cf_mat_logreg


# odds ratios
final_model <- logreg_fit$final_fit
lambda_opt  <- logreg_fit$final_param["lambda"]
coefficients <- coef(final_model, s = lambda_opt)

# Calcola gli Odds Ratio per ogni classe 
odds_ratios_logreg <- lapply(coefficients, function(class_coef) {
  coef_matrix <- as.matrix(class_coef)
  data.frame(
    Log_Odds = coef_matrix[, 1],
    Odds_Ratio = exp(coef_matrix[, 1])
  )
})

#print(odds_ratios_logreg)

# Represent odds ratios
class_labels <- names(coefficients)

coef_list = coefficients

coef_long <- map_dfr(class_labels, function(cls) {
  cm <- as.matrix(coef_list[[cls]])   # (p+1) x 1
  tibble(
    feature  = rownames(cm),
    estimate = cm[, 1],
    class    = cls
  )
}) %>%
  filter(feature != "(Intercept)") %>% 
  mutate(
    odds_ratio = exp(estimate),
    class      = factor(class, levels = class_labels)
  )

# Problem: since features are on different scales, a 1 unit increase can't have the same
# meaning for each feature. 
# Odds ratio considering betas not standardized can't be directly compared
# ex: 
summary(as.data.frame(datasets$means_logCV_ratios)$ratio_l4l5)
summary(as.data.frame(datasets$means_logCV_ratios)$layer5)

# Consider changes in 1 standard deviation for the features
feature_sds <- apply(datasets$means_logCV_ratios, 2, sd)
sd_table <- tibble(feature = names(feature_sds), sd = feature_sds)

coef_long_std <- coef_long %>%
  left_join(tibble(feature = names(feature_sds), sd = feature_sds), by = "feature") %>%
  mutate(estimate_std = estimate * sd)

# The formulation in glmnet does not have a reference class
# you need to reconduct to the healthy class as reference setting
healthy_coefs <- coef_long %>%
  filter(class == "Healthy") %>%
  dplyr::select(feature, healthy_estimate = estimate)

coef_long_ref_std <- coef_long %>%
  left_join(healthy_coefs, by = "feature") %>%
  left_join(sd_table, by = "feature") %>%
  mutate(
    estimate_vs_healthy     = estimate - healthy_estimate,  # log-odds vs Healthy
    estimate_vs_healthy_std = estimate_vs_healthy * sd    # standardized
  )

# Filter out Healthy (would be all zeros by construction)
coef_long_ref_std_plot <- coef_long_ref_std %>%
  filter(class != "Healthy") %>%
  mutate(class = droplevels(class))

my_cols = c(
  "MH"     = "#F28E2B", 
  "DR"     = "mediumpurple", 
  "CSR"    = "#E15759", 
  "AMD"    = "seagreen3",
  "Healthy" = "dodgerblue"
)

log_odds_plot <- ggplot(coef_long_ref_std_plot,
                          aes(x = estimate_vs_healthy_std, y = feature, fill = class)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.4) +
  geom_col(width = 0.7, alpha = 0.85) +
  facet_wrap(~ class, ncol = 4) +
  scale_fill_manual(values = my_cols[names(my_cols) != "Healthy"],
                    guide = "none") +
  labs(
    title    = "Standardized log-odds versus Healthy",
    subtitle = "Each bar: change in log-odds (disease vs. Healthy) per 1 SD increase in the feature",
    x = "Standardized log-odds vs. Healthy (per 1 SD)",
    y = NULL
    # caption = "Bars are directly comparable across features. Positive: feature raises odds of disease over Healthy."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

print(log_odds_plot)
ggsave("figures/classification/log_odds.png", plot = log_odds_plot, width = 8, height = 6, 
       units = "in", dpi = 300)


# ------------------------------------------------------------------------------
# Variable importance plot for the best - performing random forest
rf_var_importance <- plot_rf_importance(final_mods[["rf_means_logCV_ratios"]])
rf_var_importance
ggsave("figures/classification/rf_var_importance.png", plot = rf_var_importance, width = 8, height = 6, 
       units = "in", dpi = 300)

# ------------------------------------------------------------------------------
# Test set prediction

# For each (model, dataset) combination in "results', we extract:
#   - the final fitted model 
#   - the tuning parameter at which it was fit
#   - dataset name, model name, the original key
#

# Everything is  stored under named keys in a list called "final_models", 
# mirroring the keys in "results"

final_models <- list()

for (key in names(results)) {
  entry <- results[[key]]
  fit   <- entry$fit
  
  # nestcv objects store the final-fit object differently depending on the
  # underlying model
  # class(fit)
  
  if (inherits(fit, "nestcv.glmnet")) {
    final_models[[key]] <- list(
      dataset       = entry$dataset,
      model         = entry$model,
      type          = "glmnet",
      final_fit     = fit$final_fit,
      final_param   = fit$final_param,
      lambda        = fit$final_param["lambda"],
      alpha         = fit$final_param["alpha"]
    )
    
  } else if (inherits(fit, "nestcv.train")) {
    # nestcv.train wraps caret::train; final_fit is a caret train object
    final_models[[key]] <- list(
      dataset       = entry$dataset,
      model         = entry$model,
      type          = "caret",
      final_fit     = fit$final_fit,
      final_param   = fit$final_param
    )
  } 
}

saveRDS(final_models, file = "results/final_models.rds")

final_mods <- readRDS("results/final_models.rds")

y_test <- readRDS("datasets/y_test.rds")
test_datasets <- readRDS("datasets/test_datasets.rds")
class_ord_path <- readRDS("results/class_ord_from_healthy.rds")

# Dispatches to the correct predict method based on the model object class
predict_final_model <- function(entry, x_test, y_test) {
  
  fit <- entry$final_fit
  
  # caret-trained models
  if (inherits(fit, "train")) {
    preds <- predict(fit, newdata = x_test)
  }
  
  # nestcv.train wrappers (from RF, SVM, kNN)
  # nestcv.train stores the deployable model in $final_fit (which is
  # a caret train object)
  else if (inherits(fit, "nestcv.train")) {
    preds <- predict(fit$final_fit, newdata = x_test)
  }
  
  # cv.glmnet (from nestcv.glmnet penalized logistic regression) 
  # cv.glmnet has its own predict method 
  # pass s = "lambda.min" for the cross-validated optimal lambda 
  # chosen during the final CV
  else if (inherits(fit, "cv.glmnet")) {
    raw   <- predict(fit, newx = x_test, s = "lambda.min", type = "class")
    preds <- factor(raw[, 1], levels = levels(y_test))
  }
  
  preds
}

# Evaluate every model on its corresponding test dataset
test_results <- lapply(names(final_mods), function(key) {
  entry <- final_mods[[key]]
  
  # Select the matching test set features
  if (entry$dataset %in% names(test_datasets)) {
    x_test <- test_datasets[[entry$dataset]]}
  
  list(
    dataset = entry$dataset,
    model   = entry$model,
    eval    = evaluate_final_model(entry, x_test, y_test)
  )
})
names(test_results) <- names(final_mods)

# Returns predictions, truth, confusion matrix, and metrics
evaluate_final_model <- function(entry, x_test, y_test) {
  
  preds <- predict_final_model(entry, x_test, y_test)
  truth <- y_test
  cm    <- caret::confusionMatrix(preds, truth)
  
  list(
    predictions   = preds,
    truth         = truth,
    confusion_mat = cm,
    metrics = list(
      # Aggregate metrics
      macro_F1            = mean(cm$byClass[, "F1"],                na.rm = TRUE),
      balanced_acc        = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
      accuracy            = cm$overall["Accuracy"],
      kappa               = cm$overall["Kappa"],
      
      # Per-class — useful for diagnostic tables
      per_class_F1        = cm$byClass[, "F1"],
      per_class_precision = cm$byClass[, "Pos Pred Value"],   # caret's name for precision
      per_class_recall    = cm$byClass[, "Sensitivity"]        # caret's name for recall
    )
  )
}

summary_test <- function(results) {
  do.call(rbind, lapply(results, function(r) {
    m <- r$eval$metrics
    tibble::tibble(
      dataset         = r$dataset,
      model           = r$model,
      macro_F1        = m$macro_F1,
      balanced_acc    = m$balanced_acc,
      accuracy        = unname(m$accuracy),
      kappa           = unname(m$kappa)
    )
  }))
}

perclass_metric <- function(results, metric_name = "per_class_F1",
                               class_ord_path = "results/class_ord_from_healthy.rds") {
  class_order <- readRDS(class_ord_path)
  
  tbl <- do.call(rbind, lapply(results, function(r) {
    vec <- r$eval$metrics[[metric_name]]
    names(vec) <- sub("Class: ", "", names(vec))
    tibble::tibble(dataset = r$dataset, model = r$model, !!!vec)
  }))
  
  tbl %>% dplyr::select(dataset, model, dplyr::all_of(class_order))
}

summary_df <- summary_test(test_results)

cat("\n── HEADLINE TABLE: aggregate metrics ──\n")
print(summary_df |>
        dplyr::arrange(dplyr::desc(macro_F1), dplyr::desc(balanced_acc)) |>
        print(n = Inf))

cat("\n── BY DATASET: how does feature complexity affect each model? ──\n")
print(summary_df |>
        dplyr::group_by(dataset) |>
        dplyr::arrange(dplyr::desc(macro_F1), .by_group = TRUE) |>
        print(n = Inf))

cat("\n── PER-CLASS F1: which diseases are easy/hard? ──\n")
perclassF1_df  <- perclass_metric(test_results, "per_class_F1")
print(perclassF1_df, n = Inf)

cat("\n── PER-CLASS PRECISION: when the model predicts X, is it right? ──\n")
perclassPrec_df <- perclass_metric(test_results, "per_class_precision")
print(perclassPrec_df, n = Inf)

cat("\n── PER-CLASS RECALL: of the true X cases, how many are caught? ──\n")
perclassRec_df  <- perclass_metric(test_results, "per_class_recall")
print(perclassRec_df, n = Inf)








