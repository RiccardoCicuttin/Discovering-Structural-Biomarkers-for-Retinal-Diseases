source("R/packages.R")

train_features = readRDS("datasets/train_features.rds")

df_list_by_class = split(train_features[, 3:19], train_features$label)
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

class_ord_from_healthy <- c("Healthy", names(ranking_AI))
saveRDS(class_ord_from_healthy, "results/class_ord_from_healthy.rds")
