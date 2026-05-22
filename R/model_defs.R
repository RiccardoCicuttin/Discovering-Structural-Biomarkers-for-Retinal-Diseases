source("R/packages.R")

# Each model will have a specific function that performs nested cv
# Fit functions: 

# 1. Multiclass penalized logistic regression
fit_logreg <- function(x, y, folds) {
  
  # hyperparameter to be tuned in the elastic net
  alphas = c(0, 0.1, 0.25, 0.5, 0.75, 1)
  
  nestcv.glmnet(
    y = y, x = x,
    family = "multinomial",
    alphaSet = alphas,
    standardize = TRUE, 
    outer_folds = folds,
    n_outer_folds = 10,
    n_inner_folds = 10,
    pass_outer_folds = TRUE, 
    cv.cores = parallel::detectCores(logical = FALSE))
  # pass_outer_folds uses the outer folds computed externally to fit the final model
  # instead of new ones, useful for reproducibility.
  # final CV (on the entire dataset, runs once at the end) mirrors the inner CVs 
  # you needs n_outer_folds = n_inner_folds otherwise the metrics you compute from 
  # the nested procedure do not refer to the actual method you used to fit the model
  # (nested cv is just a performance assessment)
}

# Setting object to pass to the caret::train function
# which is used inside the inner loop of nested cv
inner_ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE, # required for logLoss
  # the model must give P(Y=k|x) for each class
  summaryFunction = mnLogLoss, # multinomial logLoss
  savePredictions = "final",
  allowParallel = FALSE
)

# 2. Support vector machines
fit_svm <- function(x, y, folds) {
  
  # this is a commonly used heuristic to tune svm hyperparameters
  tg <- expand.grid(
    C = 2^seq(-3, 7, by = 2),
    sigma = 2^seq(-5, 1, by = 2)
    # sigma since caret "svmRadial" uses
    # K(x, x') = exp( -sigma * ||x - x'||^2 )
  )
  
  nestcv.train(
    y = y, x = x,
    method = "svmRadial",
    tuneGrid = tg,
    trControl = inner_ctrl,
    metric = "logLoss",
    preProcess = c("center", "scale"),
    outer_folds = folds,
    n_outer_folds = 10,
    n_inner_folds = 10,
    pass_outer_folds = TRUE, 
    cv.cores = parallel::detectCores(logical = FALSE)
  )
}

# 3. KNN
fit_knn = function(x, y, folds){
  
  nestcv.train(
    y = y, 
    x = x,
    method = "knn", 
    metric = "logLoss",             
    preProcess = c("center", "scale"), 
    tuneLength = 3,                     
    trControl = inner_ctrl,
    outer_folds = folds,
    n_outer_folds = 10,
    n_inner_folds = 10,
    tuneGrid = data.frame(k = seq(5, 21, by=2)),
    pass_outer_folds = TRUE, 
    cv.cores = parallel::detectCores(logical = FALSE)
  )
}

# 4. Random forest
fit_rf <- function(x, y, folds) {
  
  tg <- expand.grid(
    # number of randomly selected covariates for each tree
    mtry = sort(unique(c(2, 4, 6, max(2, floor(sqrt(ncol(x))))))), # sqrt(p) usually the default
    # splits are chosen to maximize the decrease in Gini impurity
    splitrule = "gini",
    # minimum number of observations a node must contain to be eligible for further splitting
    min.node.size = c(1, 5) # 1 is the default
  )
  
  nestcv.train(
    y = y, x = x,
    # ranger is a faster alternative to "rf" 
    method = "ranger", 
    tuneGrid = tg,
    trControl = inner_ctrl,
    metric = "logLoss",
    outer_folds  = folds,
    n_outer_folds = 10,
    num.trees = 500,
    importance = "impurity",
    cv.cores = parallel::detectCores(logical = FALSE)
  )
}