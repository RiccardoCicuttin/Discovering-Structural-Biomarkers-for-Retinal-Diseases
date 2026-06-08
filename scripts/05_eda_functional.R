source("R/packages.R")
source("R/fda_utils.R")

df_raw_train <- readRDS("datasets/df_raw_train.rds")
layer_name_map <- readRDS("results/layer_name_map.rds")

# create a list containing data of a single layer
raw_train_by_layer = list()
for (l in unique(df_raw_train$names)) {
  raw_train_by_layer[[l]] = df_raw_train %>% filter(names == l)
  raw_train_by_layer[[l]] = as.data.frame(raw_train_by_layer[[l]])
}

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


# Smoothing
# create common basis for each
basis <- create.bspline.basis(c(0,1), 60, norder=4)
# number of basis functions
K = 100
# search for the optimal smoothing lambda over a logarithmic grid
lambda_grid = 10^seq(-4, 4, by = 0.5)
# results of the smoothing
df_smooth_by_lyr <- lapply(c(1:6), function(i){
  smooth_layer(raw_train_by_layer[[i]], K, lambda_grid, basis)
})
# example first patient of layer 1:
# df_smooth_by_lyr[[1]]$smoothed[[1]]$smooth_res 
saveRDS(df_smooth_by_lyr, "results/df_smooth_by_lyr.rds")
  

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
# reattach IDs and labels — row order is preserved
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
# Access: mean_funs_byclass[[layer_idx]][[class_name]] → funData with 1 curve
  
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


mfpca_fit$scores       # n x M matrix: the scores xi_{ik}
mfpca_fit$vectors      # coefficient matrix for the eigenfunction expansion
mfpca_fit$meanFunction # multiFunData: the cross-layer mean
















































