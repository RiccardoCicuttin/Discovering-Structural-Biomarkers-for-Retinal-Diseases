source("R/packages.R")

smooth_patient_single_lam <- function(patient_lyr, k, lam, basis){
  
  # basis has to be chosen outside the function
  
  meta_cols <- c("names", "patient_id", "label")
  qual <- patient_lyr[, meta_cols]
  
  # the sampling grid is different among patients of different classes
  num_cols <- setdiff(colnames(patient_lyr), meta_cols)
  y <- as.numeric(patient_lyr[, num_cols])
  y <- y[!is.na(y)]
  
  m <- length(y)
  t_grid <- seq(0, 1, length.out = m)
  
  fdpar <- fdPar(basis, Lfdobj = 2, lambda = lam)
  smooth_res <- smooth.basis(t_grid, y, fdpar)
  
  list(
    smooth_res = smooth_res$fd,
    #lambda = lam,
    label = qual$label,
    id = qual$patient_id
    #layer = qual$names
  )
}

gcv_patient_one_lyr <- function(patient_lyr, k, lambdas, basis){
  
  meta_cols <- c("names", "patient_id", "label")
  qual <- patient_lyr[, meta_cols]
  
  num_cols <- setdiff(colnames(patient_lyr), meta_cols)
  y <- as.numeric(patient_lyr[, num_cols])
  y <- y[!is.na(y)]
  
  m <- length(y)
  t_grid <- seq(0, 1, length.out = m)
  
  # find the gcv associated to each smoothing parameter lambda 
  gcvs <- sapply(lambdas, function(lam){
    fdpar <- fdPar(basis, Lfdobj = 2, lambda = lam)
    smth <- smooth.basis(t_grid, y, fdpar)
    smth$gcv
  })
  
  gcvs
}

smooth_layer <- function(layer_df, k, lambdas, basis){
  
  gcvs_list <- lapply(seq_len(nrow(layer_df)), function(i) {
    gcv_patient_one_lyr(layer_df[i, , drop = FALSE], k, lambdas, basis)
  })
  
  # build the gcvs dataframe: 
  # each row is a patient
  # each column is the gcv associated to a certain lambda
  gcvs_df <- do.call(rbind, gcvs_list)
  gcvs_means = colMeans(gcvs_df)
  
  # lambda optimum is chosen as the one with the minimum mean
  lam_opt = lambdas[which.min(gcvs_means)]
  
  # then smooth with this chosen lambda
  results <- lapply(seq_len(nrow(layer_df)), function(i) {
    smooth_patient_single_lam(layer_df[i, , drop = FALSE], k, lam_opt, basis)
  })
  
  list(smoothed = results, lambda_opt = lam_opt)
}


# plot single layer profile
plot_patient_profile <- function(df, class, id, k, basis, lambda,
                                 show_points = TRUE) {
  
  patient_long <- df %>%
    filter(label == class, patient_id == id)
  
  if (nrow(patient_long) == 0)
    warning("No rows found for class = '", class, "', patient_id = ", id)
  
  patient_long <- patient_long %>%
    mutate(names = fct_inorder(names)) %>%
    dplyr::select(-label, -patient_id) %>%
    pivot_longer(
      cols      = where(is.numeric),
      names_to  = "Raw_Column",
      values_to = "Thickness"
    ) %>%
    filter(!is.na(Thickness)) %>%
    group_by(names) %>%
    mutate(
      Location = row_number(),
      t = seq(0, 1, length.out = n())
    ) %>%
    ungroup() %>%
    dplyr::select(-Raw_Column)
  
  smooth_one <- function(y) {
    m      <- length(y)
    t_grid <- seq(0, 1, length.out = m)
    fdpar  <- fdPar(basis, Lfdobj = 2, lambda = lambda)
    smth   <- smooth.basis(t_grid, y, fdpar)
    as.numeric(eval.fd(t_grid, smth$fd))
  }
  
  patient_long <- patient_long %>%
    group_by(names) %>%
    mutate(Thickness_smooth = smooth_one(Thickness)) %>%
    ungroup()
  
  patient_long <- patient_long %>%
    group_by(t) %>%
    arrange(names, .by_group = TRUE) %>%
    mutate(Cumulative_raw = cumsum(Thickness)) %>%
    ungroup()
  
  n_layers       <- nlevels(patient_long$names)
  fill_palette   <- viridisLite::turbo(n_layers)
  point_palette  <- colorspace::darken(rev(fill_palette), amount = 0.4)
  names(point_palette) <- levels(patient_long$names)
  
  p <- ggplot(patient_long, aes(x = t)) +
    geom_area(
      aes(y = Thickness_smooth, fill = names),
      stat = "identity",
      position = "stack",
      alpha = 0.85,
      color = "white",
      linewidth = 0.2
    ) +
    scale_fill_viridis_d(option = "turbo") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.50", "0.75", "1"),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x    = "Normalised spatial position",
      y    = "Cumulative Thickness (pixel)",
      fill = "Retinal Layer"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position  = "right",
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.minor = element_blank()
    )
  
  if (show_points) {
    p <- p +
      geom_point(
        aes(y = Cumulative_raw, color = names),
        size = 0.8,
        alpha = 0.7,
        shape = 16,
        show.legend = FALSE
      ) +
      scale_color_manual(values = point_palette)
  }
  
  p
}













