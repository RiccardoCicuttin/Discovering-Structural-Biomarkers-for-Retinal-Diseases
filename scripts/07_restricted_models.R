source("R/packages.R", echo=FALSE)
source("R/feature_utils.R")
source("R/model_defs.R")
source("R/eval_utils.R")

# ------------------------------------------------------------------------------
# Restricted feature sets and evaluation on the held-out test set
#
# Three progressively stronger restrictions of the feature set selected in
# scripts/04_baseline_models.R (means_logCV_ratios) are built
# (and then evaluated on the test set):
#   full              all 6 layers, full domain (baseline)
#   central           all 6 layers, central region only
#   layers56          layers 5 and 6, full domain
#   central_layers56  layers 5 and 6, central region only
#
# Every configuration is refitted, not re-applied: the same model classes go
# through the same nested-CV tuning on the training set (same outer folds, same
# seed, same grids), so each one gets its own hyperparameters, and the final fit
# of that procedure is predicted once on the test set
#
# Thhe test set is used only to report performance
# With 86 held-out patients (9 AMD, 10 CSR, 15 DR, 38 Healthy, 14 MH) 

central_range <- c(0.25, 0.75)   # central retina, in normalized domain units
target_layers <- c(5, 6)

df_raw_train <- readRDS("datasets/df_raw_train.rds")
df_raw_test  <- readRDS("datasets/df_raw_test.rds")

train_datasets <- readRDS("datasets/train_datasets.rds")
test_datasets  <- readRDS("datasets/test_datasets.rds")

y_train <- readRDS("datasets/y_train.rds")
y_test  <- readRDS("datasets/y_test.rds")

# the outer folds of the baseline chapter, so that the nested CV is the same one
folds <- readRDS("datasets/folds.rds")

mean_cols  <- paste0("layer",       1:6)
ratio_cols <- paste0("ratio_l",     c("1l2","2l3","3l4","4l5","5l6"))
logCV_cols <- paste0("logCV_layer", 1:6)


# ------------------------------------------------------------------------------
# Central region: features recomputed on [0.25, 0.75]
#
# restrict_domain() sets to NA the measurements whose normalized position falls
# outside the window, so the feature functions of scripts/01_preprocessing.R now
# summarise the central region only. sd, and hence logCV, is computed afresh
# over the retained positions

df_central_train <- restrict_domain(df_raw_train, central_range)
df_central_test  <- restrict_domain(df_raw_test,  central_range)

train_logCV_central <- add_logCV(compute_features(df_central_train))
test_logCV_central  <- add_logCV(compute_features(df_central_test))

# the list name is kept identical to the full-domain one, so that the model keys
# (e.g. "logreg_means_logCV_ratios") match and the two chapters stay comparable
train_datasets_central <- list(
  means_logCV_ratios = as_x(train_logCV_central, c(mean_cols, logCV_cols, ratio_cols))
)

test_datasets_central <- list(
  means_logCV_ratios = as_x(test_logCV_central,  c(mean_cols, logCV_cols, ratio_cols))
)


# ------------------------------------------------------------------------------
# Layers 5 and 6: the features that involve one of the target layers
#
# A feature involves a layer when that layer's index appears in its name, so
# this keeps "layer5", "layer6", their logCV, and the ratios that touch either
# one ("ratio_l4l5", "ratio_l5l6")

involves <- function(cols, layers) {
  cols[sapply(cols, function(cl) {
    idx <- as.integer(unlist(regmatches(cl, gregexpr("[0-9]", cl))))
    any(idx %in% layers)
  })]
}

cols_56 <- c(involves(mean_cols,  target_layers),
             involves(logCV_cols, target_layers),
             involves(ratio_cols, target_layers))

# plus the patient-specific correlation between the two layers, one column of
# the flattened 6x6 inter-layer correlation matrix ("L1-2", ..., "L5-6")
cor_col_56 <- paste0("L", target_layers[1], "-", target_layers[2])

# computed here for both splits rather than read from the stored corr_pat_spec,
# so that train and test come out of the same call. On the central frames
# patient_layer_corr() correlates the two profiles over the retained window
# only, which makes the restricted correlation a different feature from the
# full-domain one.
cor_full_train    <- patient_layer_corr(df_raw_train)
cor_full_test     <- patient_layer_corr(df_raw_test)
cor_central_train <- patient_layer_corr(df_central_train)
cor_central_test  <- patient_layer_corr(df_central_test)

# the marginal features are subset from the datasets built above, never
# recomputed; the correlation column is appended by row position, both following
# the patient order of df_raw (blocks of 6 consecutive rows, one per layer)
bind_cor <- function(x, corr_mat) {
  cbind(x[, cols_56, drop = FALSE],
        corr_mat[, cor_col_56, drop = FALSE])
}

train_datasets_layers56 <- list(
  layers56 = bind_cor(train_datasets[["means_logCV_ratios"]], cor_full_train)
)

test_datasets_layers56 <- list(
  layers56 = bind_cor(test_datasets[["means_logCV_ratios"]], cor_full_test)
)


# ------------------------------------------------------------------------------
# Central region AND layers 5-6: the two restrictions composed
# Marginal features from the central-region dataset, correlation from the
# central-region raw frame, so that no feature refers to the periphery

train_datasets_central56 <- list(
  central_layers56 = bind_cor(train_datasets_central[["means_logCV_ratios"]],
                              cor_central_train)
)

test_datasets_central56 <- list(
  central_layers56 = bind_cor(test_datasets_central[["means_logCV_ratios"]],
                              cor_central_test)
)

saveRDS(train_datasets_central,   "datasets/train_datasets_central.rds")
saveRDS(test_datasets_central,    "datasets/test_datasets_central.rds")
saveRDS(train_logCV_central,      "datasets/train_logCV_central.rds")
saveRDS(test_logCV_central,       "datasets/test_logCV_central.rds")
saveRDS(train_datasets_layers56,  "datasets/train_datasets_layers56.rds")
saveRDS(test_datasets_layers56,   "datasets/test_datasets_layers56.rds")
saveRDS(train_datasets_central56, "datasets/train_datasets_central_layers56.rds")
saveRDS(test_datasets_central56,  "datasets/test_datasets_central_layers56.rds")
saveRDS(cor_central_train,        "datasets/cor_pat_spec_central_train.rds")
saveRDS(cor_central_test,         "datasets/cor_pat_spec_central_test.rds")

# y_train / y_test are unchanged: same patients, same order


# ------------------------------------------------------------------------------
# The four configurations, in order of increasing restriction

configs <- list(
  full = list(
    train = train_datasets[["means_logCV_ratios"]],
    test  = test_datasets[["means_logCV_ratios"]]
  ),
  central = list(
    train = train_datasets_central[["means_logCV_ratios"]],
    test  = test_datasets_central[["means_logCV_ratios"]]
  ),
  layers56 = list(
    train = train_datasets_layers56[["layers56"]],
    test  = test_datasets_layers56[["layers56"]]
  ),
  central_layers56 = list(
    train = train_datasets_central56[["central_layers56"]],
    test  = test_datasets_central56[["central_layers56"]]
  )
)

config_levels <- names(configs)

# share of the raw measurements each configuration still sees: the window keeps
# the positions with (i-1)/(m-1) in [0.25, 0.75], the layer restriction keeps 2
# of the 6 layers, and the two multiply
meta_cols    <- c("names", "patient_id", "label")
spatial_cols <- setdiff(colnames(df_raw_train), meta_cols)

frac_domain <- sum(!is.na(df_central_train[, spatial_cols])) /
  sum(!is.na(df_raw_train[, spatial_cols]))
frac_layers <- length(target_layers) / 6

data_kept <- c(full             = 1,
               central          = frac_domain,
               layers56         = frac_layers,
               central_layers56 = frac_domain * frac_layers)


# ------------------------------------------------------------------------------
# Refit and evaluate
# Only the two model classes carried forward from the baseline chapter. The test
# matrices are never passed to a fit function: models[[mod_name]] is called on
# cfg$train only, and cfg$test appears exclusively inside evaluate_final_model(),
# after the fit is complete

models <- list(
  logreg = fit_logreg,
  rf     = fit_rf
)

test_results <- list()

for (cfg_name in config_levels) {
  cfg <- configs[[cfg_name]]

  for (mod_name in names(models)) {
    key <- paste(mod_name, cfg_name, sep = "_")
    message("Fitting ", key, " (", ncol(cfg$train), " features) ...")

    set.seed(26) # reproducibility: same seed as the baseline chapter
    fit <- models[[mod_name]](cfg$train, y_train, folds)

    # extract_final_model() lives in R/eval_utils.R: it knows how nestcv.glmnet
    # and nestcv.train each store their final fit
    entry <- extract_final_model(list(fit = fit, dataset = cfg_name, model = mod_name))

    test_results[[key]] <- list(
      dataset  = cfg_name,
      model    = mod_name,
      n_feat   = ncol(cfg$train),
      fit      = fit,
      final    = entry,
      cv_eval  = nestcv_modeval(fit, folds),   # training-set estimate, for context
      eval     = evaluate_final_model(entry, cfg$test, y_test)
    )
  }
}


# ------------------------------------------------------------------------------
# Results: one table per model class
# Tthe comparison is full-vs-restricted within a model class
# cv_macro_F1 is the nested-CV estimate on the training set, reported alongside
# the test metrics because it rests on 265 patients rather than 86

summary_df <- summary_test(test_results) %>%
  dplyr::mutate(
    n_feat      = sapply(test_results, function(r) r$n_feat),
    data_kept   = unname(data_kept[dataset]),
    cv_macro_F1 = sapply(test_results, function(r) r$cv_eval$metrics$macro_F1),
    dataset     = factor(dataset, levels = config_levels)
  ) %>%
  dplyr::select(dataset, model, n_feat, data_kept,
                cv_macro_F1, macro_F1, balanced_acc, accuracy, kappa) %>%
  dplyr::arrange(model, dataset)

for (mod_name in names(models)) {
  cat("\n\n== ", toupper(mod_name), " -- test set (", length(y_test),
      " patients), by increasing restriction ==\n", sep = "")
  print(summary_df %>%
          dplyr::filter(model == mod_name) %>%
          dplyr::select(-model),
        n = Inf, digits = 3)
}

# the tuned hyperparameters, one set per configuration
cat("\n-- hyperparameters retuned per configuration --\n")
for (key in names(test_results)) {
  entry <- test_results[[key]]$final
  if (entry$type == "glmnet") {
    cat(sprintf("%-25s alpha = %.2f, lambda = %.4f\n",
                key, entry$alpha, entry$lambda))
  } else {
    # sapply(as.character) rather than unlist(): bestTune stores splitrule as a
    # factor, which unlist() would turn into its integer code
    bt <- entry$final_fit$bestTune
    cat(sprintf("%-25s %s\n", key,
                paste(names(bt), sapply(bt, as.character),
                      sep = " = ", collapse = ", ")))
  }
}


# ------------------------------------------------------------------------------
# Per-class F1
# The headline metrics can stay flat while the composition shifts underneath, so
# the per-class breakdown is reported as well.
# (order of classes is based on distances of the covariance matrices from class
# healthy)

perclass_F1 <- perclass_metric(test_results, "per_class_F1") %>%
  dplyr::mutate(dataset = factor(dataset, levels = config_levels)) %>%
  dplyr::arrange(model, dataset)

for (mod_name in names(models)) {
  cat("\n-- ", toupper(mod_name), ": per-class F1 on the test set --\n", sep = "")
  print(perclass_F1 %>%
          dplyr::filter(model == mod_name) %>%
          dplyr::select(-model),
        n = Inf, digits = 3)
}

class_order <- readRDS("results/class_ord_from_healthy.rds")

perclass_long <- perclass_F1 %>%
  tidyr::pivot_longer(dplyr::all_of(class_order),
                      names_to = "class", values_to = "F1") %>%
  dplyr::mutate(
    class = factor(class, levels = class_order),
    model = factor(model, levels = names(models),
                   labels = c("Logistic regression", "Random forest"))
  )

perclass_heatmap <- ggplot(perclass_long, aes(x = class, y = dataset, fill = F1)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", F1)), size = 3.6) +
  facet_wrap(~ model, ncol = 1) +
  scale_fill_gradient2(
    low      = "white",
    mid      = "#FDDBC7",
    high     = "#B2182B",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "F1"
  ) +
  scale_y_discrete(limits = rev(config_levels)) +
  labs(
    title    = "Per-class F1 on the held-out test set",
    subtitle = "Rows: increasing restriction of the feature set (top: full, bottom: layers 5-6 on the central region)",
    x = "True class",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid    = element_blank(),
    strip.text    = element_text(face = "bold"),
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.text     = element_text(size = 11)
  )

print(perclass_heatmap)
ggsave("figures/classification/perclass_F1_restrictions.png", plot = perclass_heatmap,
       width = 9, height = 7, units = "in", dpi = 300)


# ------------------------------------------------------------------------------
saveRDS(test_results, "results/confirmation_test_results.rds")
saveRDS(summary_df,   "results/confirmation_test_summary.rds")

summary_df <- readRDS("results/confirmation_test_summary.rds")
write.csv(summary_df, "results/confirmation_test_summary.csv", row.names = FALSE)

