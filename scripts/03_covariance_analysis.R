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


#patient specific correlation between layers
y_train<-readRDS("datasets/y_train.rds")
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
class(cor_pat_spec_train)
train_datasets <- saveRDS(train_datasets , "datasets/train_datasets.rds")

boxplot(cor_pat_spec)

hclust.s=hclust(dist(cor_pat_spec),method="single")
hclust.a=hclust(dist(cor_pat_spec),method="average")
hclust.c=hclust(dist(cor_pat_spec,method="manhattan"),method="complete")

plot(hclust.c,hang=-0.1,xlab='',labels=F,cex=0.6,sub='')

cluster.c=cutree(hclust.c,k=3)
plot(cor_pat_spec,col=cluster.c + 1L, pch=19)
