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
# Final fitted models

# For each (model, dataset) combination in "results', we extract:
# - the final fitted model
# - the tuning parameter at which it was fit
# - dataset name, model name, the original key
#
# Everything is  stored under named keys in a list called "final_models",
# mirroring the keys in "results"
#
# extract_final_model() lives in R/eval_utils.R: it knows how nestcv.glmnet and
# nestcv.train each store their final fit. It is shared with the test-set
# evaluation of scripts/07_restricted_models.R, so both chapters deploy the same
# object in the same way

final_models <- Filter(Negate(is.null), lapply(results, extract_final_model))

saveRDS(final_models, file = "results/final_models.rds")

final_mods <- readRDS("results/final_models.rds")


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
