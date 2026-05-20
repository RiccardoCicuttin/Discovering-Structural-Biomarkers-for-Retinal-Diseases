library(ggplot2)
library(GGally)
library(patchwork)
library(dplyr)
library(corrplot)
library(MVN)
library(car)
library(mvtnorm)
library(ggridges)
library(forcats)

# Personalized colors
my_cols = c(
  "MH"     = "#F28E2B", 
  "DR"     = "mediumpurple", 
  "CSR"    = "#E15759", 
  "AMD"    = "seagreen3",
  "Healthy" = "dodgerblue"
)

retina_df_means = retina_df[, c(1:6, 18)]

# Correlation matrices
png("figures/exploratory_data_analysis/corrplot.png", width = 2000, height = 1400, res = 200)
par(mfrow = c(2,3))

for (m in unique(retina_df_means$label)) {
  
  df = subset(retina_df_means, label == m)[, 1:6]
  cor_mat = cor(df, use = "pairwise.complete.obs")
  
  corrplot(cor_mat, method = "color",addCoef.col = "black",number.cex = 1,tl.col = "#2c3e50",tl.srt = 45 )
  title(m, line = 3, cex.main = 1.5)
}

df=subset(retina_df_means)[,1:6]
cor_mat =cor(df, use = "pairwise.complete.obs")
corrplot(cor_mat, method = "color",addCoef.col = "black",number.cex = 1,tl.col = "#2c3e50",tl.srt = 45 )
title("Global", line = 3, cex.main = 1.5)

# Analysis of correlation matrices based on spd-matrix distances
library(CovTools)
library(abind)
library(ggdendro)
df_list_by_class = split(retina_df[, 1:17], retina_df$label)
corr_m_list <- lapply(df_list_by_class, function(df) {
  as.matrix(cor(df, use = "pairwise.complete.obs"))
})

mats <- list(Healthy = corr_m_list$Healthy,
             AMD     = corr_m_list$AMD,
             CSR     = corr_m_list$CSR,
             DR      = corr_m_list$DR,
             MH      = corr_m_list$MH)

D_AI <- CovDist(abind(mats, along = 3), method = "AIRM")   # affine-invariant Riemannian metric
D_LE <- CovDist(abind(mats, along = 3), method = "LERM")   # log-Euclidean

rownames(D_AI) <- colnames(D_AI) <- names(mats)
rownames(D_LE) <- colnames(D_LE) <- names(mats)

d_from_healthy_LE <- D_LE["Healthy", ]
d_from_healthy_LE <- d_from_healthy_LE[names(d_from_healthy_LE) != "Healthy"]
ranking_LE <- sort(d_from_healthy_LE)
ranking_LE

d_from_healthy_AI <- D_AI["Healthy", ]
d_from_healthy_AI <- d_from_healthy_AI[names(d_from_healthy_AI) != "Healthy"]
ranking_AI <- sort(d_from_healthy_AI)
ranking_AI

hc <- hclust(as.dist(D_AI), method = "average")
ggdendrogram(hc, rotate = FALSE) + labs(
  title = "Hierarchical clustering of classes by covariance structure",
  x = "", 
  y = "AIM distance for SPD matrices") + theme_minimal(base_size = 14) + 
  theme(
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 13),
    axis.text.x  = element_text(size = 13),   # class names
    axis.text.y  = element_text(size = 11)
  )



# Modified format for plotting purposes
data_long <- retina_df_means %>%
  select(label, starts_with("layer")) %>%
  pivot_longer(cols = starts_with("layer"), 
               names_to = "layer_name", 
               values_to = "thickness") %>%
  mutate(label = factor(label, levels = names(my_cols)))


# Ridgeplots
ridgeplots = ggplot(data_long, aes(x = thickness, y = label, fill = label)) +
  geom_density_ridges(
    scale = 1.5, 
    alpha = 1, 
    color = "black", 
    rel_min_height = 0.005,
    quantile_lines = TRUE, 
    quantiles = 2,
    vline_color = "grey40",
    vline_linetype = "dashed"
  ) +
  facet_wrap(~ layer_name, scales = "free_x") + 
  #scale_fill_brewer(palette = "my_cols") +
  scale_fill_manual(values = my_cols) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey85", linewidth = 0.25),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(vjust = 0)
  ) + 
  labs(
    title = "Density estimates for the mean width of retinal layers",
    x = "Thickness measured in pixels",
    y = "Diagnosis" 
  )

ggsave("figures/exploratory_data_analysis/ridgeplots.png", plot = ridgeplots, width = 8, height = 6, 
      units = "in", dpi = 300)


# Boxplots
boxplots = ggplot(data_long, aes(x = label, y = thickness, fill = label)) +
  geom_boxplot(color = "black", outlier.shape = NULL, outlier.fill = "white", 
               outlier.size = 1.5, median.linewidth = 0.6) +
  facet_wrap(~ layer_name, ncol = 3, scales = "free_y") +
  #scale_fill_brewer(palette = "mycols") +
  scale_fill_manual(values = my_cols) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.border = element_rect(color = "black", fill = NA), 
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  ) +
  labs(title = "Retinal thickness per layer",
       x = "Diagnosis",
       y = "Thickness (pixel)")

ggsave("figures/exploratory_data_analysis/boxplots.png", plot = boxplots, width = 8, height = 6, 
       units = "in", dpi = 300)


# Testing normality

# Univariate normality (qqplots)
qqplots = ggplot(data_long, aes(sample = thickness, color = label)) +
  # sample quantiles vs theoretical ones
  stat_qq(alpha = 0.6, size = 1.2) +
  
  stat_qq_line(color = "black", linetype = "dashed", linewidth = 0.5) +
  
  facet_grid(layer_name ~ label, scales = "free_y") +
  
  scale_color_manual(values = my_cols) +
  theme_minimal() +
  theme(
    legend.position = "none",
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
    
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text.x = element_text(face = "bold", size = 10), 
    strip.text.y = element_text(face = "bold", size = 10, angle = 270), 
    
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 8, color = "black")
  ) +
  labs(
    title = "Univariate normality: layers vs diagnosis",
    x = "Theoretical normal quantiles",
    y = "Sample quantiles (thickness in pixel)"
  )

ggsave("figures/exploratory_data_analysis/qqplots.png", plot = qqplots, width = 8, height = 6, 
       units = "in", dpi = 300)


# Multivariate normality
labels=levels(as.factor(retina_df_means$label)) # there are five diagnosis (4 diseases + healthy)
for(i in 1:5){ 
  print(mvn(data = retina_df_means %>% filter(label == labels[i]) %>% select(-label))$multivariate_normality)
}
# no class can be considered normally distributed

# Box-Cox transform
lambda = powerTransform(retina_df_means[,-7])

retina_df_means_transf=bind_cols(
  bcPower(retina_df_means[1], lambda$lambda[1]),
  bcPower(retina_df_means[2], lambda$lambda[2]),
  bcPower(retina_df_means[3], lambda$lambda[3]),
  bcPower(retina_df_means[4], lambda$lambda[4]),
  bcPower(retina_df_means[5], lambda$lambda[5]),
  bcPower(retina_df_means[6], lambda$lambda[6]),
  retina_df_means[,7]
)

labels_tr=levels(as.factor(retina_df_means_transf$label))
for(i in 1:5){
  print(mvn(data = retina_df_means_transf %>% filter(label == labels[i]) %>% select(-label))$multivariate_normality)
}
# even after the suggested Box-Cox transform, there is no evidence for each class
# that data are normally distributed 


# PERMANOVA
library(vegan)
adonis2(as.matrix(retina_df[,-18])~retina_df$label, permutations = 9999)
# class membership explains 32.9% of the total multivariate variance

# Outlier analysis

# Find outliers based on Mahalanobis distance 
x_bar = colMeans(retina_df_means[,-7])
S = cov(retina_df_means[,-7])
d2 = mahalanobis(retina_df_means[,-7], center = x_bar, cov = S)
retina_df_means_outliers = retina_df_means[which(d2 > qchisq(0.95, df = 6)),]

retina_df_means_no_outs = retina_df_means[which(d2 <= qchisq(0.95, df = 6)),]
for(i in 1:5){ 
  print(mvn(data = retina_df_means_no_outs %>% filter(label == labels[i]) %>% select(-label))$multivariate_normality)
}
# even without outliers, there is no evidence for each class 
# that data are normally distributed 


# Plots of average retina profile per class
retina_plots_list =  vector(mode='list', length=5)
labels = c("AMD", "CSR", "DR", "MH", "Healthy")
for(i in 1:5){
  long_profiles <- read_excel(path, sheet = i) %>%
    select(-image) %>% 
    mutate(names = fct_inorder(names)) %>%
    group_by(names) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
    mutate(label = "normal") %>%
    select(-label) %>% 
    pivot_longer(
      cols = -names,                 
      names_to = "Raw_Column_Name",         
      values_to = "Mean_Thickness"   
    ) %>%
    #mutate(names = factor(names, levels = original_names, labels = new_names)) %>%
    # Group by the layer so our location index restarts at 1 for each layer
    group_by(names) %>%
    # row_number() generates the strict sequential sequence (1 to M) 
    # based on the original column order, ignoring the actual column names completely
    mutate(Location = row_number()) %>%
    ungroup() %>%
    # Clean the dataset by removing the raw character names
    select(-Raw_Column_Name)
  
  retina_plots_list[[i]] = ggplot(long_profiles, aes(x = Location, y = Mean_Thickness, fill = names)) +
    geom_area(alpha = 0.85, color = "white", linewidth = 0.2) +
    scale_fill_viridis_d(option = "turbo") + 
    labs(
      title = paste("Average Topographical Profile:", labels[i]),
      subtitle = "Aggregated across sequential spatial locations",
      x = "Spatial Location (Sequential Index)",
      y = "Cumulative Thickness (pixel)",
      fill = "Retinal Layer"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.minor = element_blank()
    )
  
}

#retina_plots_list[[4]]


# ingle patient retina plot
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

# Check skewness of sd features
sd_cols <- datasets[["means_sd_ratios"]][, grepl("^sd_", colnames(datasets[["means_sd_ratios"]]))]
apply(sd_cols, 2, function(x) moments::skewness(x, na.rm = TRUE))










`




