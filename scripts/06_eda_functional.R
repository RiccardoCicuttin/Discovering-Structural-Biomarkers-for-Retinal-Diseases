source("R/packages.R")
source("scripts/01_preprocessing.R")

# Single patient retina plot
plot_patient_profile <- function(df, class, id) {
  
  patient_long <- df %>%
    filter(label == class, patient_id == id)
  
  if (nrow(patient_long) == 0)
    warning("No rows found for class = '", class, "', patient_id = ", id,
            ". Check spelling and that patient_id was correctly parsed.")
  
  patient_long <- patient_long %>%
    mutate(names = fct_inorder(names)) %>%
    select(-label, -patient_id) %>%
    pivot_longer(
      cols      = where(is.numeric),
      names_to  = "Raw_Column",
      values_to = "Thickness"
    ) %>%
    filter(!is.na(Thickness)) %>%
    group_by(names) %>%
    mutate(Location = row_number()) %>%
    ungroup() %>%
    select(-Raw_Column)
  
  # ── Smoothing block: comment out the next 8 lines to plot raw values ──────
  patient_long <- patient_long %>%
    group_by(names) %>%
    mutate(
      Thickness = tryCatch(
        predict(loess(Thickness ~ Location, span = 0.25)),
        error = function(e) {
          warning("loess failed for layer '", unique(names), "', using raw values. ", e$message)
          Thickness
        }
      )
    ) %>%
    ungroup()
  # ── End smoothing block ───────────────────────────────────────────────────
  
  ggplot(patient_long, aes(x = Location, y = Thickness, fill = names)) +
    geom_area(alpha = 0.85, color = "white", linewidth = 0.2) +
    scale_fill_viridis_d(option = "turbo") +
    labs(
      title    = paste0("Topographical profile: ", class, ", patient ", id),
      subtitle = "Individual patient, aggregated across sequential spatial locations",
      x        = "Spatial Location (Sequential Index)",
      y        = "Cumulative Thickness (pixel)",
      fill     = "Retinal Layer"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position  = "right",
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.minor = element_blank()
    )
}

# Usage
plot_patient_profile(df_orig, "Healthy", 14)











