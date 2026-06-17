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
  y = "AIRM distance for SPD matrices") + theme_minimal(base_size = 14) + 
  theme(
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 13),
    axis.text.x  = element_text(size = 13),   # class names
    axis.text.y  = element_text(size = 11)
  )

class_ord_from_healthy <- c("Healthy", names(ranking_AI))
saveRDS(class_ord_from_healthy, "results/class_ord_from_healthy.rds")


# ------------------------------------------------------------------------------
# Patient specific correlation between layers
y_train<-readRDS("datasets/y_train.rds")

# build train set
df_raw_train<-readRDS("datasets/df_raw_train.rds")
matr=matrix(0,nrow=265,ncol=15)
j=1
for(i in seq(1,1590,by=6)){
  df = subset(df_raw_train[i:(i+5),1:880])
  df = t(df)
  df_num=na.omit(df)
  cor_mat = cor(df_num)
  val_sup=t(cor_mat)[lower.tri(t(cor_mat),diag=FALSE)]
  matr[j,]=val_sup
  j=j+1
}

v=c("L1-2","L1-3","L1-4","L1-5","L1-6","L2-3","L2-4","L2-5","L2-6","L3-4","L3-5","L3-6","L4-5","L4-6","L5-6")
cor_pat_spec_train=as.data.frame(matr)

colnames(cor_pat_spec_train)=v
train_datasets=readRDS("datasets/train_datasets.rds")
as_x <- function(df, cols) as.matrix(df[, cols])
train_datasets[["corr_pat_spec"]] <- as_x(cor_pat_spec_train)
train_datasets <- saveRDS(train_datasets , "datasets/train_datasets.rds")

# build test set
df_raw_test<-readRDS("datasets/df_raw_test.rds")
matr=matrix(0,nrow=86,ncol=15)
j=1
for(i in seq(1,516,by=6)){
  df = subset(df_raw_test[i:(i+5),1:880])
  df = t(df)
  df_num=na.omit(df)
  cor_mat = cor(df_num)
  val_sup=t(cor_mat)[lower.tri(t(cor_mat),diag=FALSE)]
  matr[j,]=val_sup
  j=j+1
}

v=c("L1-2","L1-3","L1-4","L1-5","L1-6","L2-3","L2-4","L2-5","L2-6","L3-4","L3-5","L3-6","L4-5","L4-6","L5-6")
cor_pat_spec_test=as.data.frame(matr)

colnames(cor_pat_spec_test)=v
test_datasets=readRDS("datasets/test_datasets.rds")
as_x <- function(df, cols) as.matrix(df[, cols])
test_datasets[["corr_pat_spec"]] <- as_x(cor_pat_spec_test)
test_datasets <- saveRDS(test_datasets , "datasets/test_datasets.rds")

boxplot(cor_pat_spec)


# ------------------------------------------------------------------------------
# Clustering considering the correlation as features
hclust.s=hclust(dist(cor_pat_spec),method="single")
hclust.a=hclust(dist(cor_pat_spec),method="average")
hclust.c=hclust(dist(cor_pat_spec,method="manhattan"),method="complete")

plot(hclust.c,hang=-0.1,xlab='',labels=F,cex=0.6,sub='')

cluster.c=cutree(hclust.c,k=3)
plot(cor_pat_spec,col=cluster.c + 1L, pch=19)


# ------------------------------------------------------------------------------
# Clustering considering the distances between patient correlation matrices
# we use the AIM distance for spd matrices

y_train      <- readRDS("datasets/y_train.rds")
df_raw_train <- readRDS("datasets/df_raw_train.rds")

# Build patient-specific SPD matrices and label alignment
# Build patient_meta from df_raw_train
# The data is structured as 6 consecutive rows per patient (one per layer). 
# Group by (patient_id, label) preserves the patient ordering from df_raw_train 
# label is kept for later evaluations
patient_meta_train <- df_raw_train %>%
  distinct(patient_id, label) %>%
  mutate(row_idx = row_number())

patient_meta_train == train_patients

n_patients <- nrow(patient_meta_train)

# Build the 6 × 6 SPD array, one slice per patient
# we need to compute the full 6 x 6 matrix to work in the SPD geometry
spd_array <- array(NA, dim = c(6, 6, n_patients))

for (i in seq_len(n_patients)) {
  # 6 consecutive rows in df_raw_train correspond to one patient's 6 layers
  rows <- ((i - 1) * 6 + 1):(i * 6)
  patient_block <- df_raw_train[rows, 1:880] %>% t() %>% as.matrix()
  patient_block <- na.omit(patient_block)
  
  spd_array[, , i] <- cor(patient_block, use = "complete.obs")
}

# Compute pairwise distances on SPD manifold
D_AIRM   <- as.matrix(CovDist(spd_array, method = "AIRM")) # affine invariant Riemaniann metric
D_LogEuc <- as.matrix(CovDist(spd_array, method = "LERM"))   # log-Euclidean

# Hierarchical clustering (Ward linkage) 
hc_AIRM   <- hclust(as.dist(D_AIRM),   method = "ward.D")
# hc_LogEuc <- hclust(as.dist(D_LogEuc), method = "average")

# Cut at k=5
clust_AIRM   <- cutree(hc_AIRM,   k = 5)
clust_LogEuc <- cutree(hc_LogEuc, k = 5)

# Confusion matrix 
confusion_AIRM <- table(Cluster = clust_AIRM, Label = y_train)
print(confusion_AIRM)

# Per-cluster purity: for each cluster, what fraction belongs to its
# dominant class? 
cluster_purity <- function(clusters, labels) {
  cm <- table(clusters, labels)
  sum(apply(cm, 1, max)) / sum(cm)
}
cat("Cluster purity (AIRM):  ", round(cluster_purity(clust_AIRM,   y_train), 3), "\n")
cat("Cluster purity (LogEuc):", round(cluster_purity(clust_LogEuc, y_train), 3), "\n")

# Plot confusion matrix as a heatmap
pheatmap(
  as.matrix(confusion_AIRM),
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  scale           = "none",
  color           = colorRampPalette(c("white", "#2166AC"))(50),
  display_numbers = TRUE,
  number_format   = "%d",
  number_color    = "grey20",
  main            = "Clusters (AIRM) × True disease labels",
  fontsize_row    = 11,
  fontsize_col    = 11,
  border_color    = "grey85",
  cellwidth       = 50,
  cellheight      = 40
)

# Visualize the dendrogram colored by true labels
dend <- as.dendrogram(hc_AIRM)
labels_in_dend_order <- y_train[hc_AIRM$order]

# Apply class colors to leaves
my_cols = c(
  "MH"     = "#F28E2B", 
  "DR"     = "mediumpurple", 
  "CSR"    = "#E15759", 
  "AMD"    = "seagreen3",
  "Healthy" = "dodgerblue"
)
leaf_cols <- my_cols[as.character(labels_in_dend_order)]

dend <- dend %>%
  set("leaves_col", leaf_cols) %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.8) %>%
  set("labels", "")   # hide individual patient labels (too crowded)

plot(dend, main = "Hierarchical clustering (AIRM, Ward linkage)")
legend("topright", legend = names(my_cols), col = my_cols, pch = 19, cex = 0.9, bty = "n")


# ------------------------------------------------------------------------------
# Plot the correlation matrix of a single patient
plot_patient_cor <- function(df_raw_train, patient_id, label, my_cols = NULL) {
  
  # Find the rows belonging to this patient
  patient_rows <- df_raw_train %>%
    mutate(row_idx = row_number()) %>%
    filter(patient_id == !!patient_id, label == !!label) %>%
    pull(row_idx)
  
  if (length(patient_rows) == 0)
    stop("No rows found for patient_id = ", patient_id, ", label = ", label)
  
  # Extract the 6 × 880 measurement matrix and transpose to 880 × 6 for cor()
  patient_block <- df_raw_train[patient_rows, 1:880] %>% t() %>% as.matrix()
  patient_block <- na.omit(patient_block)
  
  cor_mat <- cor(patient_block, use = "complete.obs")
  rownames(cor_mat) <- colnames(cor_mat) <- paste0("layer", 1:6)
  
  pheatmap(
    cor_mat,
    cluster_rows    = FALSE,
    cluster_cols    = FALSE,
    scale           = "none",
    color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks          = seq(-1, 1, length.out = 101),
    border_color    = "grey85",
    cellwidth       = 50,
    cellheight      = 50,
    fontsize        = 12,
    angle_col       = 0,
    display_numbers = TRUE,
    number_format   = "%.2f",
    number_color    = "grey20",
    main            = paste0("Inter-layer correlations: ", label, ", patient ", patient_id)
  )
}

# Usage
plot_patient_cor(df_raw_train, patient_id = 103, label = "Healthy")
plot_patient_cor(df_raw_train, patient_id = 11, label = "AMD")
plot_patient_cor(df_raw_train, patient_id = 51, label = "MH")
plot_patient_cor(df_raw_train, patient_id = 41, label = "DR")

























