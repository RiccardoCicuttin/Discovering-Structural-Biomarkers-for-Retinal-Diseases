source("R/packages.R")

# Feature engineering on the raw (wide) thickness profiles

# compute features from original dataframe
compute_features <- function(df_split) {

  # dplyr::select is qualified here: MASS (attached via fda.usc) masks select()
  spatial_cols <- df_split %>%
    dplyr::select(where(is.numeric), -patient_id) %>%
    colnames()

  df_summary <- df_split %>%
    mutate(
      mean_val = rowMeans(across(all_of(spatial_cols)), na.rm = TRUE),
      sd_val   = apply(across(all_of(spatial_cols)), 1, sd, na.rm = TRUE)
    ) %>%
    dplyr::select(patient_id, label, names, mean_val, sd_val)

  # Layer order from the data
  layer_order <- df_split %>% pull(names) %>% unique()

  df_wide <- df_summary %>%
    pivot_wider(
      id_cols     = c(patient_id, label),
      names_from  = names,
      values_from = c(mean_val, sd_val)
    ) %>%
    # Rename from original names to layer1...layer6 / sd_layer1...6
    rename_with(~ paste0("layer",    1:6), all_of(paste0("mean_val_", layer_order))) %>%
    rename_with(~ paste0("sd_layer", 1:6), all_of(paste0("sd_val_",   layer_order)))

  df_wide %>%
    mutate(
      ratio_l1l2 = layer1 / layer2,
      ratio_l2l3 = layer2 / layer3,
      ratio_l3l4 = layer3 / layer4,
      ratio_l4l5 = layer4 / layer5,
      ratio_l5l6 = layer5 / layer6
    )
    # %>% relocate(patient_id, label, .after = last_col())
}

add_logCV <- function(df) {
  for (k in 1:6) {
    mn <- df[[paste0("layer",    k)]]
    sd <- df[[paste0("sd_layer", k)]]
    df[[paste0("logCV_layer", k)]] <- log(sd / mn)
  }
  df
}

# as_x selects the requested columns from a data frame and returns a plain
# numeric matrix
# glmnet and nestcv require x to be a numeric matrix
as_x <- function(df, cols) as.matrix(df[, cols])
