source("R/packages.R")
source("R/model_defs.R")
source("R/eval_utils.R")

# 3 different sets of features
# saved them in a list for practicality
datasets = readRDS("datasets/train_datasets.rds")

y = readRDS("datasets/y_train.rds")

# the k (external) folds need to be common between models
set.seed(2026)
folds = createFolds(y, k = 10, returnTrain = FALSE)
# save them since they will need to be used also by other models
saveRDS(folds, "datasets/folds.rds")

# if one wants to add a model, he just has to write the fit_* function inside R/model_defs.R
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
                           eval = nestcv_modeval(fit))
  }
}

saveRDS(results, "results/baseline_results.rds")

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


# Analysis of selected model - penalized logistic regression

logreg_fit<- results[["logreg_means_logcv_ratios"]]$fit   
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
plot_confusion_matrix(results[["logreg_means_logcv_ratios"]])

#oddsratio
fit_logreg_best <- results[["logreg_means_logCV_ratios"]]$fit
final_model <- fit_logreg_best$final_fit

# 3. Estrai i coefficienti usando il lambda ottimale sintonizzato da nestcv
lambda_opt <- fit_logreg_best$final_param["lambda"]
coefficients <- coef(final_model, s = lambda_opt)

# 4. Calcola gli Odds Ratio per ogni classe (multinomiale)
odds_ratios_logreg <- lapply(coefficients, function(class_coef) {
  coef_matrix <- as.matrix(class_coef)
  data.frame(
    Log_Odds = coef_matrix[, 1],
    Odds_Ratio = exp(coef_matrix[, 1])
  )
})

print(odds_ratios_logreg)



























