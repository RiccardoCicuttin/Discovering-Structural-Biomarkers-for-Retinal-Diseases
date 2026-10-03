source("R/packages.R")

# Macro-averaged F1 from a caret confusionMatrix.
# caret returns NA for a class it never predicts, because precision is 0/0.
# The limiting value of F1 there is 0 (recall is 0), so such a class scores 0
# rather than being dropped from the average: a model that ignores a class must
# be penalised for it. A class absent from the reference altogether has no F1 to
# score and is still dropped.
macro_f1 <- function(cm) {
  f1      <- cm$byClass[, "F1"]
  support <- colSums(cm$table)          # reference counts, cols = truth
  f1[is.na(f1) & support > 0] <- 0
  mean(f1, na.rm = TRUE)
}

# Per-class F1 on the same convention as macro_f1(), for diagnostic tables.
per_class_f1 <- function(cm) {
  f1      <- cm$byClass[, "F1"]
  support <- colSums(cm$table)
  f1[is.na(f1) & support > 0] <- 0
  f1
}

# function to obtain evaluation metrics from results of a nested cv
# this gives an estimate of the performance of the procedure used to 
# choose hyperparaeters and fit the model
nestcv_modeval <- function(obj, folds) {
  
  truth <- obj$output$testy
  # those are the actual labels, just permutated due to the procedure
  
  preds <- factor(obj$output$predy, levels = levels(truth))
  # Outer CV predictions (unbiased)
  # After the loop, every patient has exactly one prediction, 
  # made by the model whose training fold did not include them. 
  # obj$output is a synthetic dataset of size n in which every
  # prediction is out of sample for the model that produced it
  
  # Pooled confusion matrix and metrics across all out-of-fold predictions
  cm <- caret::confusionMatrix(preds, truth)
  
  # Per-fold metrics: compute on each fold separately, then take mean and sd
  # these measures are not of a single model that actually exists
  # they're estimated by pooling 10 procedure instances 
  fold_metrics <- do.call(rbind, lapply(seq_along(folds), function(k) {
    fold_idx <- folds[[k]]
    preds_f  <- preds[fold_idx]
    truth_f  <- truth[fold_idx]
    cm_f     <- caret::confusionMatrix(preds_f, truth_f)
    
    tibble::tibble(
      fold         = k,
      macro_F1     = macro_f1(cm_f),
      balanced_acc = mean(cm_f$byClass[, "Balanced Accuracy"], na.rm = TRUE)
    )
  }))
  
  list(
    predictions   = preds,
    truth         = truth,
    confusion_mat = cm,
    fold_metrics  = fold_metrics,
    metrics = list(
      # Aggregate metrics: pooled across all predictions
      macro_F1        = macro_f1(cm),
      balanced_acc    = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
      per_class_F1    = per_class_f1(cm),
      # Variability across folds: meaningful uncertainty estimate for CV
      macro_F1_sd     = sd(fold_metrics$macro_F1),
      balanced_acc_sd = sd(fold_metrics$balanced_acc)
    )
  )
}


# Create tibble with nested cv metrics
summary_cv <- function(results) {
  
  tbl <- do.call(rbind, lapply(results, function(r) {
    m <- r$eval$metrics
    tibble::tibble(
      dataset         = r$dataset,
      model           = r$model,
      macro_F1        = m$macro_F1,
      macro_F1_sd     = m$macro_F1_sd,
      balanced_acc    = m$balanced_acc,
      balanced_acc_sd = m$balanced_acc_sd
      # add further metrics here: kappa = m$kappa, etc.
    )
  }))
}


# Create tibble with nested F1 scores per class
perclass_F1_cv <- function(results,
                            class_ord_path = "results/class_ord_from_healthy.rds") {
  
  # Load class ordering: Healthy first, then diseases by distance from Healthy.
  # This ordering is meaningful: classes further from Healthy in covariance
  # space should in principle be easier to classify — the table makes this
  # pattern immediately visible.
  class_order <- readRDS(class_ord_path)
  
  tbl <- do.call(rbind, lapply(results, function(r) {
    f1_vec        <- r$eval$metrics$per_class_F1
    names(f1_vec) <- sub("Class: ", "", names(f1_vec))   # strip caret's prefix
    tibble::tibble(dataset = r$dataset, model = r$model, !!!f1_vec)
  }))
  
  # Reorder class columns according to the covariance-based ordering
  tbl %>% dplyr::select(dataset, model, dplyr::all_of(class_order))
}


# Confusion matrix plot
plot_confusion_matrix <- function(result_entry, title = NULL) {
  
  preds <- result_entry$eval$predictions
  truth <- result_entry$eval$truth
  cm    <- result_entry$eval$confusion_mat
  
  # Auto-generate title from dataset and model name if not supplied
  if (is.null(title))
    title <- paste0(toupper(result_entry$model), " — ", result_entry$dataset)
  
  cm_df <- as.data.frame(cm$table) %>%
    dplyr::group_by(Reference) %>%
    dplyr::mutate(
      pct     = Freq / sum(Freq),
      pct_lbl = ifelse(Freq > 0,
                       paste0(Freq, "\n(", scales::percent(pct, accuracy = 1), ")"),
                       "0")
    ) %>%
    dplyr::ungroup()
  
  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = pct_lbl), size = 3.8, lineheight = 1.2) +
    geom_tile(
      data = dplyr::filter(cm_df, Prediction == Reference),
      aes(x = Reference, y = Prediction),
      fill = NA, color = "#2166AC", linewidth = 1.5
    ) +
    scale_fill_gradient2(
      low      = "white",
      mid      = "#FDDBC7",
      high     = "#B2182B",
      midpoint = 0.3,
      limits   = c(0, 1),
      labels   = scales::percent,
      name     = "Row %"
    ) +
    labs(
      title    = title,
      subtitle = paste0(
        "Macro F1 = ",      round(result_entry$eval$metrics$macro_F1,     3),
        "   |   Balanced accuracy = ", round(result_entry$eval$metrics$balanced_acc, 3)
      ),
      x = "True class",
      y = "Predicted class"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid      = element_blank(),
      axis.text       = element_text(size = 12),
      plot.title      = element_text(face = "bold", size = 15),
      plot.subtitle   = element_text(size = 11, color = "grey40"),
      legend.position = "right"
    ) +
    coord_fixed() + scale_y_discrete(limits = rev)
}


# variable importance plot for random forests
plot_rf_importance <- function(final_mods_entry, top_n = NULL,
                               title = NULL) {
  
  fit <- final_mods_entry$final_fit
  
  # nestcv.train wraps a caret::train inside $final_fit 
  # can be unwrapped if needed
  if (inherits(fit, "nestcv.train")) fit <- fit$final_fit
  
  # caret::varImp returns importance as $importance (data.frame).
  # For multinomial classification, columns are per-class importance plus an
  # overall column. We use the overall mean across classes as a single score.
  imp <- caret::varImp(fit)$importance
  
  imp_df <- tibble::tibble(
    feature    = rownames(imp),
    importance = rowMeans(imp) # avg across classes        
  ) %>%
    dplyr::arrange(dplyr::desc(importance))
  
  # Optionally keep only the top N features
  if (!is.null(top_n)) imp_df <- imp_df %>% dplyr::slice_head(n = top_n)
  
  # Normalise to 0-100 (relative to maximum)
  # (Introduction to statistical learning uses this convention)
  imp_df <- imp_df %>%
    dplyr::mutate(importance = 100 * importance / max(importance))
  
  if (is.null(title)) {
    title <- paste0("Random forest variable importance — ",
                    final_mods_entry$dataset)
  }
  
  ggplot(imp_df, aes(x = importance,
                     y = forcats::fct_reorder(feature, importance))) +
    geom_col(fill = "#B2182B", color = "black", width = 0.7) +
    labs(
      title = title,
      x     = "Variable importance",
      y     = NULL
    ) +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    theme_minimal(base_size = 13) +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(face = "bold")
    )
}








