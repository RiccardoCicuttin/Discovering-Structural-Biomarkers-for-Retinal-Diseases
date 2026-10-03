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


# Restrict the spatial domain of a raw (wide) dataframe to `range`
# The features produced downstream are computed by the same code as
# the full-domain ones, the domain being the only difference
restrict_domain <- function(df_split, range = c(0, 1)) {

  spatial_cols <- df_split %>%
    dplyr::select(where(is.numeric), -patient_id) %>%
    colnames()

  vals <- as.matrix(df_split[, spatial_cols])

  for (i in seq_len(nrow(vals))) {
    obs    <- which(!is.na(vals[i, ]))   # this row's measurements
    m      <- length(obs)
    t_grid <- seq(0, 1, length.out = m)
    keep   <- t_grid >= range[1] & t_grid <= range[2]
    vals[i, obs[!keep]] <- NA_real_
  }

  df_split[spatial_cols] <- as.data.frame(vals)
  df_split
}


# For each patient (a block of n_layers consecutive rows), compute the
# inter-layer correlation matrix and flatten its lower triangle into one
# row of features (e.g. "L1-2", "L1-3", ...).
#
# na.omit() drops the positions that are unmeasured in ANY layer, so the
# correlation is computed on the positions the patient's layers have in common
# On a domain-restricted frame this is the retained window since
# restrict_domain() blanks the same normalized interval in every row and the
# number of measurements is constant across a patient's layers
patient_layer_corr <- function(df_raw, n_layers = 6, n_features = 880) {
  n_patients <- nrow(df_raw) / n_layers
  pairs <- combn(n_layers, 2)
  col_names <- paste0("L", pairs[1, ], "-", pairs[2, ])

  out <- t(sapply(seq_len(n_patients), function(p) {
    rows <- ((p - 1) * n_layers + 1):(p * n_layers)
    patient_block <- na.omit(t(df_raw[rows, 1:n_features]))
    cor_mat <- cor(patient_block)
    cor_mat[lower.tri(cor_mat)]
  }))

  colnames(out) <- col_names
  out
}
