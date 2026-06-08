source("R/packages.R")
source("R/fda_utils.R")

my_cols = c(
  "MH"     = "#F28E2B", 
  "DR"     = "mediumpurple", 
  "CSR"    = "#E15759", 
  "AMD"    = "seagreen3",
  "Healthy" = "dodgerblue"
)

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

curves_long %>%
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


# Plot single layer profile
plot_patient_profile(df_raw_train, "CSR", 50, 100, basis, 1e-04, show_points = T)
# csr 50


# Plot the eigenfunctions
n_pc  <- 8
#cum_ve[8] : 0.87
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

ggplot(eig_long, aes(x = t, y = value, color = layer)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ psi_idx, ncol = 4, labeller = label_parsed) +
  scale_color_viridis_d(option = "turbo") +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  labs(
    title    = "First 8 multivariate eigenfunctions from MFPCA",
    subtitle = expression("Each panel shows " * psi[k](t) * " = (" * psi[k]^(1)(t) * ", ..., " * psi[k]^(6) (t) * ")"),
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


# Boxplot of first 8 mfpcas scores by class
n_pc <- 8

scores_long <- mfpca_features %>%
  dplyr::select(label, dplyr::all_of(paste0("MFPC", 1:n_pc))) %>%
  tidyr::pivot_longer(
    cols      = dplyr::starts_with("MFPC"),
    names_to  = "pc",
    values_to = "score"
  ) %>%
  dplyr::mutate(pc = factor(pc, levels = paste0("MFPC", 1:n_pc)))

ggplot(scores_long, aes(x = label, y = score, fill = label)) +
  geom_boxplot(color = "black", outlier.shape = NULL, outlier.fill = "white",
               outlier.size = 1.5, median.linewidth = 0.6) +
  facet_wrap(~ pc, ncol = 4, scales = "free_y") +
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
    title = "Distribution of the first 8 MFPCA scores by class",
    x     = "Diagnosis",
    y     = "Score"
  )


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
    title = "Class-specific mean thickness profiles",
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
    title = "Class-specific pointwise variance profiles",
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
    panel.grid.minor = element_blank()
  )

p_means
p_vars


# Perturbation
n_pc <- 4
# if I move the population mean in the direction of
# psi_k by the average amount class j moves

# Population mean function per layer (across all training patients)
mu_mat <- t(sapply(fd_layers, function(fl) meanFunction(fl)@X[1, ]))   # 6 × 750

# Class-mean MFPCA scores: average score per (class, PC)
class_mean_scores <- mfpca_features %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(dplyr::across(dplyr::starts_with("MFPC"), mean),
                   .groups = "drop")
# class_mean_scores: rows = classes, columns = label + MFPC1..MFPC20

# Build long tibble: one row per (PC, layer, class, grid point) for the perturbed curve
perturb_long <- do.call(rbind, lapply(1:n_pc, function(k) {
  do.call(rbind, lapply(1:6, function(lyr) {
    # Layer-lyr component of psi_k
    psi_k_lyr <- mfpca_eigfuns[[k]]@X[lyr, ]
    
    do.call(rbind, lapply(classes, function(cls) {
      # Average score of class cls on PC k
      xi_bar_jk <- class_mean_scores[[paste0("MFPC", k)]][class_mean_scores$label == cls]
      
      tibble::tibble(
        pc    = paste0("psi[", k, "]"),
        layer = paste0("layer", lyr),
        class = cls,
        t     = t_grid,
        value = mu_mat[lyr, ] + xi_bar_jk * psi_k_lyr
      )
    }))
  }))
}))

perturb_long$pc    <- factor(perturb_long$pc, levels = paste0("psi[", 1:n_pc, "]"))
perturb_long$layer <- factor(perturb_long$layer, levels = paste0("layer", 1:6))


ggplot(perturb_long, aes(x = t, y = value, color = class)) +
  geom_line(linewidth = 0.8) +
  facet_grid(pc ~ layer, scales = "free_y", labeller = labeller(pc = label_parsed)) +
  scale_color_manual(values = my_cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  labs(
    title    = "Class-specific perturbations along each MFPCA component",
    subtitle = expression("Per panel: " * hat(mu)^(l)(t) + bar(xi)[list(j,k)] %.% psi[k]^(l)(t)),
    x        = "Normalised spatial position",
    y        = "Perturbed thickness",
    color    = "Diagnosis"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

















