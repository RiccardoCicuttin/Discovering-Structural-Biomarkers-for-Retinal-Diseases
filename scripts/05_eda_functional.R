source("R/packages.R")
source("R/fda_utils.R")

df_raw_train <- readRDS("datasets/df_raw_train.rds")
layer_name_map <- readRDS("results/layer_name_map.rds")

# create a list containing data of a single layer
raw_train_by_layer <- list()
for (l in unique(df_raw_train$names)) {
  raw_train_by_layer[[l]] = df_raw_train %>% filter(names == l)
  raw_train_by_layer[[l]] = as.data.frame(raw_train_by_layer[[l]])
}

saveRDS(raw_train_by_layer, "datasets/raw_train_by_layer.rds")

# check how many patients have negative recordings per layer, divided by class
neg_counts <- do.call(rbind, lapply(1:6, function(lyr_idx) {
  
  df <- raw_train_by_layer[[lyr_idx]]
  meta_cols    <- c("names", "patient_id", "label")
  spatial_cols <- setdiff(colnames(df), meta_cols)
  
  # For each patient (row), TRUE if any spatial value is negative
  has_negative <- apply(df[, spatial_cols], 1, function(row) {
    any(row < 0, na.rm = TRUE)
  })
  
  # Count patients with any negative, by class
  tapply(has_negative, df$label, sum) |>
    as.list() |>
    tibble::as_tibble() |>
    dplyr::mutate(layer = paste0("layer", lyr_idx), .before = 1)
}))

print(neg_counts)


# ------------------------------------------------------------------------------
# Smoothing
# create common basis for each
basis <- create.bspline.basis(c(0,1), 60, norder=4)
# number of basis functions
K = 60
# possible rule of thumb: K = min(n/4, 40)  
# search for the optimal smoothing lambda over a logarithmic grid
lambda_grid = 10^seq(-4, 4, by = 0.5)
# results of the smoothing
df_smooth_by_lyr <- lapply(c(1:6), function(i){
  smooth_layer(raw_train_by_layer[[i]], K, lambda_grid, basis)
})
# example first patient of layer 1:
# df_smooth_by_lyr[[1]]$smoothed[[1]]$smooth_res 

saveRDS(df_smooth_by_lyr, "results/df_smooth_by_lyr.rds")
df_smooth_by_lyr <- readRDS("results/df_smooth_by_lyr.rds")

# ------------------------------------------------------------------------------
# MFPCA
# To use MFPCA you need a multiFunData object 
# which is a list of univariate funData objects

# To get a funData object you need:
# the domain (x-values) 
# the values of the different observations (y-values)

# domain (evaluation grid)
grid <- seq(0, 1, length.out = 750) 

# Use the previous smoothing results to build the funData
build_funData_layer <- function(layer_entry, t_grid) {
  
  # Evaluate every patient's fd on the grid
  # stack them into a matrix
  # Each row is one patient's curve computed at the 750 grid points
  mat <- t(sapply(layer_entry$smoothed, function(obs) {
    as.numeric(eval.fd(t_grid, obs$smooth_res))
  }))
  # mat: n_patients × 750 
  
  # funData constructor: argvals = grid, X = matrix 
  funData(argvals = t_grid, X = mat)
}

# fd_layers is a list of 6 funData objects, one per layer
fd_layers <- lapply(df_smooth_by_lyr, build_funData_layer, t_grid = grid)

mfd <- multiFunData(fd_layers)

# we need to save the couple (id, label) of each patient
# they can't go into the multiFunData object
patient_meta <- tibble::tibble(
  id    = sapply(df_smooth_by_lyr[[1]]$smoothed, function(obs) obs$id),
  label = sapply(df_smooth_by_lyr[[1]]$smoothed, function(obs) obs$label)
)

# The MFPCA algorithm implemented reduces the multivariate problem to a 
# univariate FPCA on each component, followed by a combination step. 
# uniExpansions specifies the basis used to represent each component's curves 
# and eigenfunctions during the univariate FPCA stage
# k is the dimension of the basis used to represent each component's eigenfunctions 
# during the univariate FPCA. It should not exceed the dimension of the smoothing
# basis used previously, since FPCA cannot recover spatial detail 
# that has already been smoothed away
uni_exp <- lapply(1:6, function(i) {
  list(
    type = "splines1D", # type of basis functions
    k = 30)
})

mfpca_fit <- MFPCA(
  mFData = mfd,
  M = 20, # number of pcs
  uniExpansions = uni_exp,
  fit = TRUE
)

# save the scores: n x M matrix
# Dataset of scores
scores_df <- as.data.frame(mfpca_fit$scores)
colnames(scores_df) <- paste0("MFPC", 1:ncol(scores_df))
# reattach IDs and labels (row order is preserved)
mfpca_features <- dplyr::bind_cols(patient_meta, scores_df)

# group means of scores
mean_scores <- mfpca_features %>% dplyr::select(-'id') %>% group_by(label) %>% 
  summarise_all(list(Mean_score = mean))

# Cumulative variance explained
ve <- mfpca_fit$values / sum(mfpca_fit$values)
cum_ve <- cumsum(ve)
print(round(cum_ve, 3))

# eigenfunctions (vectors of six functions)
# mfpca_fit$functions is a  multiFunData: the M eigenfunctions psi_1,...,psi_M
# save separately each eigenfunction
mfpca_eigfuns <- sapply(1:20, function(i){
  
  mat = t(sapply(1:6, function(j){
    mfpca_fit$functions[[j]]@X[i, ]
  }))
  
  funData(argvals = grid, X = mat)
  
})
# mfpca_fit$functions[[1]]@X[20,] == mfpca_eigfuns[[20]]@X[1,]

#mfpca_fit$scores       # n x M matrix: the scores xi_{ik}
#mfpca_fit$vectors      # coefficient matrix for the eigenfunction expansion
#mfpca_fit$meanFunction # multiFunData: the cross-layer mean


# ------------------------------------------------------------------------------
classes <- levels(factor(patient_meta$label))

# mean functions
mean_funs <- lapply(fd_layers, meanFunction)
mean_funs_byclass <- lapply(seq_along(fd_layers), function(i) {
  fd_layer <- fd_layers[[i]]
  lapply(setNames(classes, classes), function(cls) {
    rows_cls <- which(patient_meta$label == cls)
    meanFunction(fd_layer[rows_cls])
  })
})
# Access:
# mean_funs_byclass[[layer_idx]][[class_name]] 
# is a funData with 1 curve
  
# variance functions
var_funs <- lapply(fd_layers, function(fl) {
  v <- matrixStats::colVars(fl@X)
  funData(argvals = fl@argvals, X = matrix(v, nrow = 1))
})

var_funs_byclass <- lapply(seq_along(fd_layers), function(i) {
  fd_layer <- fd_layers[[i]]
  lapply(setNames(classes, classes), function(cls) {
    rows_cls <- which(patient_meta$label == cls)
    X_sub    <- fd_layer@X[rows_cls, , drop = FALSE]
    v        <- matrixStats::colVars(X_sub)
    funData(argvals = fd_layer@argvals, X = matrix(v, nrow = 1))
  })
})


# ------------------------------------------------------------------------------
# Smoothed curves evaluated on the grid, organized as an
# n x L x M array (patients x layers x positions) 
# The funData objects in fd_layers already hold each layer's 
# n x M matrix in their @X slot, so we
# just stack the 6 layers along a new middle dimension
# fd_layers[[l]]@X  is  n x M   (patients x positions) for layer l


grid   <- seq(0, 1, length.out = 750)
labels <- factor(patient_meta$label) 
n <- length(labels)
L <- length(fd_layers)                  
M  <- length(grid)                      

# Build the n x L x M array
layer_array <- array(NA_real_, dim = c(n, L, M))
for (ell in seq_len(L)) {
  layer_array[, ell, ] <- fd_layers[[ell]]@X   # n x M slice for this layer
}

# B / (B + sum_k W_k) is a trace-based descriptive ratio: the fraction
# of total functional dispersion (summed across layers and positions)
# explained by class membership
# it's computed by vegan::adonis2

# need the pairwise L2 distance matrix in the product space H = (L^2)^6
# to calculate it

# Convert each funData layer to fdata (fda.usc class)
# fd_layers[[ell]]@X is the n x M matrix already evaluated on the grid
fdata_layers <- lapply(fd_layers, function(fd_l) {
  fdata(fd_l@X, argvals = grid)
})

# metric.lp computes the n x n matrix of pairwise L2 distances for one layer
# default is L2 norm
D_sq_layers <- lapply(fdata_layers, function(fd) {
  as.matrix(metric.lp(fd, lp = 2))^2    # square to get squared distances
})

# Total squared multivariate distance: sum across layers
D_sq_total <- Reduce("+", D_sq_layers)

# Convert to a dist object (taking the square root back to actual distances)
D <- as.dist(sqrt(D_sq_total))

# PERMANOVA
set.seed(2026)
permanova_res <- adonis2(D ~ labels, permutations = 9999)
print(permanova_res)

eta2 <- permanova_res$R2[1]                   # the eta^2 effect size
cat(sprintf("\neta^2 (fraction of functional variability explained by class): %.4f\n",
            eta2))

# Pointwise wilk's lambda
# At each position t, treat the L=6 layer values as a multivariate response
# and run a one-way MANOVA across classes. 

# We plot 1 - Lambda(t) so that large = strong

compute_one_minus_wilks <- function(arr, class_labels) {
  M <- dim(arr)[3]
  out <- numeric(M)
  for (t_idx in seq_len(M)) {
    Y <- arr[, , t_idx]                      
    # manova requires a matrix response
    # class_labels is the grouping factor
    fit  <- manova(Y ~ class_labels)
    lam  <- summary(fit, test = "Wilks")$stats[1, "Wilks"]
    out[t_idx] <- 1 - lam
  }
  out
}

sep_curve <- compute_one_minus_wilks(layer_array, labels)

# Tidy data frame for plotting
sep_df <- tibble(
  location  = grid,
  separation = sep_curve
)

sep_curve <- ggplot(sep_df, aes(x = location, y = separation)) +
  geom_point(size = 1.6, color = "#1565C0", alpha = 0.7) +
  labs(
    #title    = "Pointwise multivariate class separation along the retina",
    #subtitle = expression("1 - Wilks' " * Lambda * "(t)"),
    x = "Retinal location (normalized)",
    y = expression("1 - " * Lambda * "(t)")
  ) +
  ylim(0, max(sep_curve) * 1.05) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(sep_curve)

separation_results <- list(
  eta2          = eta2,
  permanova     = permanova_res,
)
saveRDS(separation_results, "results/separation_results.rds")

ggsave("figures/functional_data_analysis/sep_curve.png", plot = sep_curve, width = 8, height = 6, 
       units = "in", dpi = 300)




























