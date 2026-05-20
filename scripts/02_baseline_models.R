library(caret)
library(nestedcv)
library(MLmetrics)
library(glmnet)
library(pbapply)
library(pheatmap)

# 3 different sets of features
# saved them in a list for practicality
datasets = list(
  # nest cv functions require a matrix in input
  means = as.matrix(retina_df_means[, -ncol(retina_df_means)]),
  means_sd = as.matrix(retina_df[, -c(13:18)]),
  means_sd_ratios = as.matrix(retina_df[, -ncol(retina_df)]),
  means_logcv_ratios = as.matrix(retina_df_logCV[, -ncol(retina_df)])
)

y = retina_df$label

set.seed(2026)

# the k (external) folds need to be common between models
folds = createFolds(y, k = 10, returnTrain = FALSE)

# function to obtain evaluation metrics from results of a nested cv
# this gives an estimate of the performance of the procedure used to 
# choose hyperparaeters and fit the model
nestcv_modeval <- function(obj) {
  
  # Outer CV predictions (unbiased)
  preds = as.factor(obj$output$predy)
  # After the loop, every patient has exactly one prediction, 
  # made by the model whose training fold did not include them. 
  # obj$output is a synthetic dataset of size n in which every
  # prediction is out of sample for the model that produced it
  truth = obj$output$testy
  # those are the actual labels, just permutated due to the procedure
  
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

# Each model will have a specific function that performs nested cv
# Fit functions: 

# 1. Multiclass penalized logistic regression
fit_logreg <- function(x, y, folds) {
  
  # hyperparameter to be tuned in the elastic net
  alphas = c(0, 0.1, 0.25, 0.5, 0.75, 1)
  
  nestcv.glmnet(
    y = y, x = x,
    family = "multinomial",
    alphaSet = alphas,
    standardize = TRUE, 
    outer_folds = folds,
    n_outer_folds = 10,
    n_inner_folds = 10,
    pass_outer_folds = TRUE, 
    cv.cores = parallel::detectCores(logical = FALSE))
    # pass_outer_folds uses the outer folds computed externally to fit the final model
    # instead of new ones, useful for reproducibility.
    # final CV (on the entire dataset, runs once at the end) mirrors the inner CVs 
    # you needs n_outer_folds = n_inner_folds otherwise the metrics you compute from 
    # the nested procedure do not refer to the actual method you used to fit the model
    # (nested cv is just a performance assessment)
}

# Setting object to pass to the caret::train function
# which is used inside the inner loop of nested cv
inner_ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE, # required for logLoss
  # the model must give P(Y=k|x) for each class
  summaryFunction = mnLogLoss, # multinomial logLoss
  savePredictions = "final",
  allowParallel = FALSE
)

# 2. Support vector machines
fit_svm <- function(x, y, folds) {
  
  # this is a commonly used heuristic to tune svm hyperparameters
  tg <- expand.grid(
    C = 2^seq(-3, 7, by = 2),
    sigma = 2^seq(-5, 1, by = 2)
    # sigma since caret "svmRadial" uses
    # K(x, x') = exp( -sigma * ||x - x'||^2 )
  )
  
  nestcv.train(
    y = y, x = x,
    method = "svmRadial",
    tuneGrid = tg,
    trControl = inner_ctrl,
    metric = "logLoss",
    preProcess = c("center", "scale"),
    outer_folds = folds,
    n_outer_folds = 10,
    n_inner_folds = 10,
    pass_outer_folds = TRUE, 
    cv.cores = parallel::detectCores(logical = FALSE)
  )
}

# 3. KNN
fit_knn = function(x, y, folds){
  
  nestcv.train(
    y = y, 
    x = x,
    method = "knn", 
    metric = "logLoss",             
    preProcess = c("center", "scale"), 
    tuneLength = 3,                     
    trControl = inner_ctrl,
    outer_folds = folds,
    n_outer_folds = 10,
    n_inner_folds = 10,
    tuneGrid = data.frame(k = seq(5, 21, by=2)),
    pass_outer_folds = TRUE, 
    cv.cores = parallel::detectCores(logical = FALSE)
  )
}

# 4. Random forest
fit_rf <- function(x, y, folds) {
  
  tg <- expand.grid(
    # number of randomly selected covariates for each tree
    mtry = sort(unique(c(2, 4, 6, max(2, floor(sqrt(ncol(x))))))), # sqrt(p) usually the default
    # splits are chosen to maximize the decrease in Gini impurity
    splitrule = "gini",
    # minimum number of observations a node must contain to be eligible for further splitting
    min.node.size = c(1, 5) # 1 is the default
  )
  
  nestcv.train(
    y = y, x = x,
    # ranger is a faster alternative to "rf" 
    method = "ranger", 
    tuneGrid = tg,
    trControl = inner_ctrl,
    metric = "logLoss",
    outer_folds  = folds,
    n_outer_folds = 10,
    num.trees = 500,
    importance = "impurity",
    cv.cores = parallel::detectCores(logical = FALSE)
  )
}

# if one wants to add a model, he just has to write the fit_ function and add the model to this list
models <- list(
  logreg = fit_logreg,
  rf = fit_rf,
  svm = fit_svm,
  knn = fit_knn
)

# list where to save the results of each model fitted on a specific dataset
results = list()

class(models[[mod_name]](datasets[[ds_name]], y, folds))

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
                           eval = nestcv_modeval(fit))
  }
}

# summary of the metrics used for the comparison
# rbind the resuts of the second argument 
summary_tbl <- do.call(rbind, lapply(results, function(r) {
  m <- r$eval$metrics
  tibble::tibble(
    dataset = r$dataset,
    model = r$model,
    macro_F1 = m$macro_F1,
    macro_F1_sd = m$macro_F1_sd,
    balanced_acc = m$balanced_acc,
    balanced_acc_sd = m$balanced_acc_sd
    # other_metric = m$*
  )
}))

#print(summary_tbl)

# order by F1 score
print(summary_tbl |>
        dplyr::arrange(dplyr::desc(macro_F1), dplyr::desc(balanced_acc)) |>
        print(n = Inf))

# check if order is consistent with more features added
print(summary_tbl |>
        dplyr::group_by(dataset) |>
        dplyr::arrange(dplyr::desc(macro_F1), .by_group = TRUE) |>
        print(n = Inf))
 
# F1 score per class
perclassF1_tbl <- do.call(rbind, lapply(results, function(r) {
  f1_vec <- r$eval$metrics$per_class_F1
  names(f1_vec) <- sub("Class: ", "", names(f1_vec))
  tibble::tibble(dataset = r$dataset, model = r$model, !!!f1_vec)
}))

class_order <- c("Healthy", names(ranking_AI))
perclassF1_tbl <- perclassF1_tbl %>%
  dplyr::select(dataset, model, all_of(class_order))

print(perclassF1_tbl)

summary_tbl[which.max(summary_tbl$macro_F1) , ]


# Analysis of selected model
logreg_fit<- results[["knn_means_logcv_ratios"]]$fit   
# note that glment fits the multinomial logistic regression in the "symmetric form"
# you do not pick a baseline class 
# you get estimates of coefficients for each class 

# heatmap of coefficents of the penalized logistic regression
all_features <- colnames(datasets[["means_sd_ratios"]])

coef_list <- coef(logreg_fit)

coef_matrix <- t(sapply(coef_list, function(ck) {
  ck <- ck[names(ck) != "(Intercept)"]   # drop intercept
  out <- setNames(rep(0, length(all_features)), all_features)
  out[names(ck)] <- ck                   # fill in the non-zero values
  out
}))

pheatmap(
  coef_matrix,
  # classes are clustered according to the correlation matrix distances
  cluster_rows    = hc, 
  cluster_cols    = FALSE,
  scale           = "none",
  color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks          = seq(-max(abs(coef_matrix)), max(abs(coef_matrix)), length.out = 101),
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
  main            = "Multinomial logistic regression coefficients",
  treeheight_row  = 90,
  margins         = c(80, 10) 
)

# confusion matrix
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

plot_confusion_matrix(results[["knn_means_logcv_ratios"]])

