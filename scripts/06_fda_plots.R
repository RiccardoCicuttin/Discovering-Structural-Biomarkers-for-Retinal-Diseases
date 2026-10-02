source("R/packages.R")
source("R/fda_utils.R")

my_cols = c(
  "MH"     = "#F28E2B", 
  "DR"     = "mediumpurple", 
  "CSR"    = "#E15759", 
  "AMD"    = "seagreen3",
  "Healthy" = "dodgerblue"
)


# ------------------------------------------------------------------------------
# Representation of each layer colored by class
# For each layer, for each patient, evaluate the fd on a common grid
# and store (layer, patient_id, label, t, value) for ggplot.
t_eval <- seq(0, 1, length.out = 200)   # common evaluation grid for plotting

curves_long <- do.call(rbind, lapply(1:6, function(lyr_idx) {
  
  layer_entry <- df_smooth_by_lyr[[lyr_idx]]
  
  do.call(rbind, lapply(layer_entry$smoothed, function(obs) {
    
    vals <- as.numeric(eval.fd(t_eval, obs$smooth_res))   
    
    tibble::tibble(
      layer = paste0("layer", lyr_idx),
      id    = obs$id,
      label = obs$label,
      t     = t_eval,
      value = vals
    )
  }))
}))

layers_smooth_plot <- curves_long %>%
  ggplot(aes(x = t, y = value, color = label, group = interaction(label, id))) +
  geom_line(alpha = 0.4, linewidth = 0.6) +    # was 0.3
  facet_wrap(~ layer, ncol = 3, scales = "free_y") +
  scale_color_manual(values = my_cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1))) +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    panel.grid.major = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey85", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text       = element_text(face = "bold")
  ) +
  labs(
    title = "Smoothed retinal layer profiles by class",
    x     = "Normalised spatial position",
    y     = "Thickness (pixels)",
    color = "Diagnosis"
  )
layers_smooth_plot

ggsave("figures/functional_data_analysis/layers_smooth_plot.png", plot = layers_smooth_plot, width = 8, height = 6, 
       units = "in", dpi = 300)

# ------------------------------------------------------------------------------
 # Plot single layer profile
df_raw <- readRDS("datasets/df_raw.rds")
h5_plot <- plot_patient_profile(df_raw, "Healthy", 5, 100, basis, 1e-04, show_points = F)
h5_plot
ggsave("figures/functional_data_analysis/smooth_healthy5.png", plot = h5_plot, width = 9, height = 4.5, 
       units = "in", dpi = 300)
amd12_plot <- plot_patient_profile(df_raw, "AMD", 12, 100, basis, 1e-04, show_points = F)
amd12_plot
ggsave("figures/functional_data_analysis/smooth_AMD12.png", plot = amd12_plot, width = 9, height = 4.5, 
       units = "in", dpi = 300)


# ------------------------------------------------------------------------------
# Plot the eigenfunctions
n_pc  <- 6
# cum_ve[n_pc] is the share of the retained (M = 10) variance these panels
# carry; print(round(cum_ve, 3)) in 05_eda_functional.R and quote it in the
# report rather than hard-coding it here.
n_lyr  <- 6

eig_long <- do.call(rbind, lapply(1:n_pc, function(k) {
  do.call(rbind, lapply(1:n_lyr, function(lyr) {
    tibble::tibble(
      psi_idx = paste0("psi[", k, "]"),
      layer   = paste0("layer", lyr),
      t       = grid,
      value   = mfpca_eigfuns[[k]]@X[lyr, ] # k-th eigenfunction, lyr-th component
    )
  }))
}))

eig_long$psi_idx <- factor(eig_long$psi_idx, levels = paste0("psi[", 1:n_pc, "]"))
eig_long$layer   <- factor(eig_long$layer,   levels = paste0("layer", 1:n_lyr))

mfpca_fun_plot <- ggplot(eig_long, aes(x = t, y = value, color = layer)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ psi_idx, ncol = 3, labeller = label_parsed) +
  scale_color_viridis_d(option = "turbo") +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  labs(
    #title    = "First 8 multivariate eigenfunctions from MFPCA",
    #subtitle = expression("Each panel shows " * psi[k](t) * " = (" * psi[k]^(1)(t) * ", ..., " * psi[k]^(6) (t) * ")"),
    x        = "Normalised spatial position",
    y        = "Eigenfunction value",
    color    = "Retinal layer"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )
mfpca_fun_plot

ggsave("figures/functional_data_analysis/mfpca_fun_plot.png", plot = mfpca_fun_plot, width = 8, height = 6, 
       units = "in", dpi = 300)

# ------------------------------------------------------------------------------
# Boxplot of first 8 mfpcas scores by class
n_pc <- 6

scores_long <- mfpca_features %>%
  dplyr::select(label, dplyr::all_of(paste0("MFPC", 1:n_pc))) %>%
  tidyr::pivot_longer(
    cols      = dplyr::starts_with("MFPC"),
    names_to  = "pc",
    values_to = "score"
  ) %>%
  dplyr::mutate(pc = factor(pc, levels = paste0("MFPC", 1:n_pc)))

mfpca_scores_boxplot <- ggplot(scores_long, aes(x = label, y = score, fill = label)) +
  geom_boxplot(color = "black", outlier.shape = NULL, outlier.fill = "white",
               outlier.size = 1.5, median.linewidth = 0.6) +
  facet_wrap(~ pc, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = my_cols) +
  theme_minimal() +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.border     = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  ) +
  labs(
    title = paste0("Distribution of the first ", n_pc, " MFPCA scores by class"),
    x     = "Diagnosis",
    y     = "Score"
  )
mfpca_scores_boxplot

ggsave("figures/functional_data_analysis/mfpca_scores_boxplot.png", plot = mfpca_scores_boxplot, width = 8, height = 6, 
       units = "in", dpi = 300)


# ------------------------------------------------------------------------------
# Plot mean and variance functions per class
classes <- levels(factor(patient_meta$label))
t_grid  <- fd_layers[[1]]@argvals[[1]]

# Build long tibble for mean functions, one row per (layer, class, grid point)
mean_long <- do.call(rbind, lapply(1:6, function(lyr) {
  do.call(rbind, lapply(classes, function(cls) {
    tibble::tibble(
      layer = paste0("layer", lyr),
      class = cls,
      t     = t_grid,
      value = mean_funs_byclass[[lyr]][[cls]]@X[1, ]
    )
  }))
}))

# Same structure for variance functions
var_long <- do.call(rbind, lapply(1:6, function(lyr) {
  do.call(rbind, lapply(classes, function(cls) {
    tibble::tibble(
      layer = paste0("layer", lyr),
      class = cls,
      t     = t_grid,
      value = var_funs_byclass[[lyr]][[cls]]@X[1, ]
    )
  }))
}))

# Plot class-mean functions
p_means <- ggplot(mean_long, aes(x = t, y = value, color = class)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ layer, ncol = 3, scales = "free_y") +
  scale_color_manual(values = my_cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  labs(
    #title = "Class-specific mean thickness profiles",
    x     = "Normalised spatial position",
    y     = "Thickness (pixel)",
    color = "Diagnosis"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# Plot class-variance functions
p_vars <- ggplot(var_long, aes(x = t, y = value, color = class)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ layer, ncol = 3, scales = "free_y") +
  scale_color_manual(values = my_cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  labs(
    #title = "Class-specific pointwise variance profiles",
    x     = "Normalised spatial position",
    y     = "Variance",
    color = "Diagnosis"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p_means
p_vars

ggsave("figures/functional_data_analysis/mean_funs.png", plot = p_means, width = 8, height = 6, 
       units = "in", dpi = 300)
ggsave("figures/functional_data_analysis/var_funs.png", plot = p_vars, width = 8, height = 6, 
       units = "in", dpi = 300)


# ------------------------------------------------------------------------------
# Perturbation: what shape change does component k encode
# the six-layer mean is perturbed by a fixed multiple of the
# component's own standard deviation sqrt(lambda_k)
# Every component is scaled by its own SD, so the panels are directly
# comparable, and the +/- pair shows the mode symmetrically in both directions.
n_pc <- 4

# Population mean function per layer (across all training patients)
# 6 × 750 matrix
mu_mat <- t(sapply(fd_layers, function(fl) meanFunction(fl)@X[1, ]))

c_sd     <- 2                  # perturbation size, in SDs of the component (1 is also conventional)
lambda_k <- mfpca_fit$values   # eigenvalues: score variance on each component

perturb_sd_long <- do.call(rbind, lapply(1:n_pc, function(k) {
  do.call(rbind, lapply(1:6, function(lyr) {
    psi_k_lyr <- mfpca_eigfuns[[k]]@X[lyr, ]
    delta     <- c_sd * sqrt(lambda_k[k]) * psi_k_lyr
    mu_lyr    <- mu_mat[lyr, ]

    tibble::tibble(
      pc    = paste0("psi[", k, "]"),
      layer = paste0("layer", lyr),
      curve = rep(c("mean", "plus", "minus"), each = length(t_grid)),
      t     = rep(t_grid, 3),
      value = c(mu_lyr, mu_lyr + delta, mu_lyr - delta)
    )
  }))
}))

perturb_sd_long$pc    <- factor(perturb_sd_long$pc,    levels = paste0("psi[", 1:n_pc, "]"))
perturb_sd_long$layer <- factor(perturb_sd_long$layer, levels = paste0("layer", 1:6))
perturb_sd_long$curve <- factor(perturb_sd_long$curve, levels = c("mean", "plus", "minus"))

mfpca_lyr_perturbations_sd <- ggplot(perturb_sd_long,
                                     aes(x = t, y = value,
                                         color = curve, linetype = curve)) +
  geom_line(linewidth = 0.8) +
  facet_grid(pc ~ layer, scales = "free_y", labeller = labeller(pc = label_parsed)) +
  scale_color_manual(
    values = c(mean = "grey30", plus = "#D62728", minus = "#1F77B4"),
    labels = c(mean = "mean", plus = "mean + c*sd", minus = "mean - c*sd")
  ) +
  scale_linetype_manual(
    values = c(mean = "solid", plus = "dashed", minus = "dashed"),
    labels = c(mean = "mean", plus = "mean + c*sd", minus = "mean - c*sd")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  labs(
    x        = "Normalised spatial position",
    y        = "Perturbed thickness",
    color    = paste0("Perturbation (c = ", c_sd, ")"),
    linetype = paste0("Perturbation (c = ", c_sd, ")")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1)
  )

mfpca_lyr_perturbations_sd

ggsave("figures/functional_data_analysis/mfpca_lyr_perturbations_sd.png", plot = mfpca_lyr_perturbations_sd, width = 8, height = 6,
       units = "in", dpi = 300)















