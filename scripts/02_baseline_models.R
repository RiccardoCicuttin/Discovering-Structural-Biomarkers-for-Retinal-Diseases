library(caret)
library(nestedcv)
library(MLmetrics)
library(glmnet)

# 3 different sets of features
# saved them in a list for practicality
datasets = list(
  # nest cv functions require a matrix in input
  means = as.matrix(retina_df_means[, -ncol(retina_df_means)]),
  means_sd = as.matrix(retina_df[, -c(13:18)]),
  means_sd_ratios = as.matrix(retina_df[, -ncol(retina_df)])
)

y = retina_df$label

set.seed(2026)

# the k (external) folds need to be common between models
folds = createFolds(y, k = 10, returnTrain = FALSE)

# function to obtain evaluation metrics from results of a nested cv
# this gives an estimate of the performance of the procedure used to 
# choose hyperparaeters and fit the model
nestcv_modeval <- function(obj) {
  
  # Outer CV predictions (unbiased)
  preds = as.factor(obj$output$predy)
  # After the loop, every patient has exactly one prediction, 
  # made by the model whose training fold did not include them. 
  # obj$output is a synthetic dataset of size n in which every
  # prediction is out of sample for the model that produced it
  truth = obj$output$testy
  # those are the actual labels, just permutated due to the procedure
  
  cm = caret::confusionMatrix(preds, truth)
  
  # those measures are not of a single model that actually exists
  # they're estimated by pooling 10 procedure instances 
  list(
    predictions = preds,
    truth = truth,
    confusion_mat = cm,
    metrics = list(
      macro_F1 = mean(cm$byClass[, "F1"], na.rm = TRUE),
      per_class_F1 = cm$byClass[, "F1"],
      balanced_acc = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE)
    )
  )
}

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
    # ( nested cv is just a performance assessment)
}

#
inner_ctrl <- trainControl(
  method          = "cv",
  number          = 10,
  classProbs      = TRUE,
  summaryFunction = mnLogLoss, #multinomial logLoss
  savePredictions = "final",
  allowParallel = FALSE
)

#
#

# if one wants to add a model, he just has to write the fit_ function and add the model to this list
models <- list(
  logreg = fit_logreg,
  rf = fit_rf,
  svm = fit_svm,
  knn = fit_knn
)

models = list(logreg = fit_logreg, rf = fit_rf)

# list where to save the results of each model fitted on a specific dataset
results = list()

# actual loop for each model on each dataset
for (ds_name in names(datasets)) {
  for (mod_name in names(models)) {
    key <- paste(mod_name, ds_name, sep = "_") # key to lookup for results, ex: logreg_means
    #message("Fitting ", key, " ...")
    set.seed(2026) # reproducibility
    fit <- models[[mod_name]](datasets[[ds_name]], y, folds) # models[[mod_name]] selects the correct fit function
    results[[key]] <- list(dataset = ds_name,
                           model   = mod_name,
                           fit = fit, 
                           eval = nestcv_modeval(fit))
  }
}

# summary of the metrics used for the comparison
# rbind the resuts of the second argument 
summary_tbl <- do.call(rbind, lapply(results, function(r) {
  m <- r$eval$metrics
  tibble::tibble(
    dataset      = r$dataset,
    model        = r$model,
    macro_F1     = m$macro_F1,
    balanced_acc = m$balanced_acc
  )
}))

print(summary_tbl)

summary_tbl[which.max(summary_tbl$macro_F1) , ]



# NESTED CROSS VALIDATION
# (abbiamo pochi dati e iperparametri da trovare)

# da capire:
# nei fold va usata una ripartizione proporzionale? come fare? possibile con questa funzione?
# problema: alcuni fold potrebbero avere quasi solo sani
# alcuni fold potrebbero non avere pazienti di una qualche classe

# da capire
inner_ctrl = trainControl(
  method = "cv", 
  number = 5, 
  classProbs = TRUE,                  # Absolutely required for multiclass metrics
  summaryFunction = multiClassSummary # THIS IS THE KEY!
)

# da capire!
nested_multiclass_rf <- nestcv.train(
  y = as.data.frame(retina_df_means)$label, 
  x = as.data.frame(retina_df_means[,-7]),
  method = "knn",                      # Using Random Forest here
  #Select your model here (e.g., "rf", "knn", "svmRadial")
  
  # Tell the algorithm to use the Mean F1 to select the best hyperparameters
  metric = "Mean_F1", # IMPORTANTE USARE QUESTA METRICA              
  
  tuneLength = 3,                     
  outer_cv = 5,                       
  trControl = inner_ctrl
)

# View the results
#print(nested_multiclass_rf)
print(nested_multiclass_rf$summary)

labels = as.data.frame(retina_df)$label

#KNN
library(class)

retina_std_means=scale(retina_df_means[,-7], center= FALSE, scale= TRUE)

k_values <- 1:10
wcss <- sapply(k_values, function(k) {
  kmeans(retina_std_means, centers = k, nstart = 20)$tot.withinss
})

plot(k_values, wcss,  #error plot
     type = "b",
     pch = 19,              
     col = "blue",          
     xlab = "Cluster (K)", 
     ylab = "Error (WCSS)",
     main = "elbow method")

#we select k=5 for the KNN


par(mfrow = c(2,3))
#layer1-2
x <- seq(min(retina_std_means[,1]), max(retina_std_means[,1]), length=200)
y <- seq(min(retina_std_means[,2]), max(retina_std_means[,2]), length=200)
xy<-expand.grid(layer1=x,layer2=y)
data.knn5 <- knn(train = retina_std_means[,1:2], test = xy, cl = retina_df_means$label, k = 5)
z <- as.numeric(data.knn5)
cl <- as.factor(retina_df_means$label)
plot(retina_std_means[,1:2], main="k-NN with k = 5", xlab='layer1', ylab='layer2', 
     pch=20, col=my_cols[as.numeric(cl)],
     cex.main=1.2)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5,3.5, 4.5), 
        drawlabels=FALSE, add=TRUE, lwd=2, col="black")

#layer2-3
x <- seq(min(retina_std_means[,2]), max(retina_std_means[,2]), length=200)
y <- seq(min(retina_std_means[,3]), max(retina_std_means[,3]), length=200)
xy<-expand.grid(layer2=x,layer3=y)
data.knn5 <- knn(train = retina_std_means[,2:3], test = xy, cl = retina_df_means$label, k = 5)
z <- as.numeric(data.knn5)
cl <- as.factor(retina_df_means$label)
plot(retina_std_means[,2:3], main="k-NN with k = 5", xlab='layer2', ylab='layer3', 
     pch=20, col=my_cols[as.numeric(cl)],
     cex.main=1.2)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5,3.5, 4.5), 
        drawlabels=FALSE, add=TRUE, lwd=2, col="black")

#layer3-4
x <- seq(min(retina_std_means[,3]), max(retina_std_means[,3]), length=200)
y <- seq(min(retina_std_means[,4]), max(retina_std_means[,4]), length=200)
xy<-expand.grid(layer3=x,layer4=y)
data.knn5 <- knn(train = retina_std_means[,3:4], test = xy, cl = retina_df_means$label, k = 5)
z <- as.numeric(data.knn5)
cl <- as.factor(retina_df_means$label)
plot(retina_std_means[,3:4], main="k-NN with k = 5", xlab='layer3', ylab='layer4', 
     pch=20, col=my_cols[as.numeric(cl)],
     cex.main=1.2)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5,3.5, 4.5), 
        drawlabels=FALSE, add=TRUE, lwd=2, col="black")

#layer4-5
x <- seq(min(retina_std_means[,4]), max(retina_std_means[,4]), length=200)
y <- seq(min(retina_std_means[,5]), max(retina_std_means[,5]), length=200)
xy<-expand.grid(layer4=x,layer5=y)
data.knn5 <- knn(train = retina_std_means[,4:5], test = xy, cl = retina_df_means$label, k = 5)
z <- as.numeric(data.knn5)
cl <- as.factor(retina_df_means$label)
plot(retina_std_means[,4:5], main="k-NN with k = 5", xlab='layer4', ylab='layer5', 
     pch=20, col=my_cols[as.numeric(cl)],
     cex.main=1.2)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5,3.5, 4.5), 
        drawlabels=FALSE, add=TRUE, lwd=2, col="black")
#layer5-6
x <- seq(min(retina_std_means[,5]), max(retina_std_means[,5]), length=200)
y <- seq(min(retina_std_means[,6]), max(retina_std_means[,6]), length=200)
xy<-expand.grid(layer5=x,layer6=y)
data.knn5 <- knn(train = retina_std_means[,5:6], test = xy, cl = retina_df_means$label, k = 5)
z <- as.numeric(data.knn5)
cl <- as.factor(retina_df_means$label)
plot(retina_std_means[,5:6], main="k-NN with k = 5", xlab='layer5', ylab='layer6', 
     pch=20, col=my_cols[as.numeric(cl)],
     cex.main=1.2)
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5,3.5, 4.5), 
        drawlabels=FALSE, add=TRUE, lwd=2, col="black")













