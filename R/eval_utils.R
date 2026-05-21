source("R/packages.R")

# function to obtain evaluation metrics from results of a nested cv
# this gives an estimate of the performance of the procedure used to 
# choose hyperparaeters and fit the model
nestcv_modeval <- function(obj) {
  
  truth = obj$output$testy
  # those are the actual labels, just permutated due to the procedure
  
  preds = factor(obj$output$predy, levels = levels(truth))
  # Outer CV predictions (unbiased)
  # After the loop, every patient has exactly one prediction, 
  # made by the model whose training fold did not include them. 
  # obj$output is a synthetic dataset of size n in which every
  # prediction is out of sample for the model that produced it
  
  cm = caret::confusionMatrix(preds, truth)
  
  # these measures are not of a single model that actually exists
  # they're estimated by pooling 10 procedure instances 
  list(
    predictions = preds,
    truth = truth,
    confusion_mat = cm,
    metrics = list(
      macro_F1 = mean(cm$byClass[, "F1"], na.rm = TRUE),
      macro_F1_sd = sd(cm$byClass[, "F1"], na.rm = TRUE),
      per_class_F1 = cm$byClass[, "F1"],
      balanced_acc = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
      balanced_acc_sd = sd(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE)
    )
  )
  # other ones can be computed since the confusion matrix gets saved
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