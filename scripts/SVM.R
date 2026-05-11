## Support Vector Machines

library(e1071)
library(caret)
library(nestedcv)
library(MLmetrics)

retina_std_means = scale(retina_df_means[,-7], center = FALSE, scale = TRUE)


inner_ctrl = trainControl(
  method = "cv", 
  number = 5, 
  classProbs = TRUE,                  # Absolutely required for multiclass metrics
  summaryFunction = multiClassSummary # THIS IS THE KEY!
)


# da capire!
nested_multiclass_rf <- nestcv.train(
  y = as.data.frame(retina_df_means)$label, 
  x = as.data.frame(retina_std_means[,-7]),
  method = "svmRadial",                      # Using Random Forest here
  #Select your model here (e.g., "rf", "knn", "svmRadial")
  
  # Tell the algorithm to use the logLoss to select the best hyperparameters
  metric = "logLoss", # IMPORTANTE USARE QUESTA METRICA              
  
  tuneLength = 3,                     
  outer_cv = 5,                       
  trControl = inner_ctrl
)

print(nested_multiclass_rf$summary)
nested_multiclass_rf$bestTunes
nested_multiclass_rf$outer_result

pred = nested_multiclass_rf$output$predy
truth = nested_multiclass_rf$output$testy

cm = caret::confusionMatrix(pred, truth)

macro_F1 = mean(cm$byClass[, "F1"], na.rm = TRUE)

## Healthy vs AMD
first_two = retina_std_means[,c(1,6)]
lab = factor(retina_df_means$label)

indices_sep = which(lab == 'Healthy' | lab == 'AMD')
dat_sep = data.frame(x = first_two[indices_sep, 1:2],
                     y = factor(lab[indices_sep]))

head(dat_sep)

# Jitter ? non so se serve
set.seed(2026)
dat_sep[,1:2] <- dat_sep[,1:2] + cbind(rnorm(nrow(dat_sep), sd=0.025))

plot(dat_sep[,1:2], main = 'first two layers (Healthy vs AMD)', xlab = 'layer1',
     ylab = 'layer2', pch = 19, col = my_cols[as.numeric(dat_sep$y)])

svm = svm(y ~., data = dat_sep, kernel = 'radial', cost = 10, gamma = 0.8, scale = FALSE)

plot(svm, dat_sep, 
     col = my_cols,
     #main = 'Linearly separable Healthy vs AMD', 
     pch = 19)


## Healthy vs CSR
first_two = retina_std_means[,c(1,6)]
lab = factor(retina_df_means$label)

indices_sep = which(lab == 'Healthy' | lab == 'CSR')
dat_sep = data.frame(x = first_two[indices_sep, 1:2],
                     y = factor(lab[indices_sep]))

head(dat_sep)

# Jitter ? non so se serve
set.seed(2026)
dat_sep[,1:2] <- dat_sep[,1:2] + cbind(rnorm(nrow(dat_sep), sd=0.025))

plot(dat_sep[,1:2], main = 'first two layers (Healthy vs CSR)', xlab = 'layer1',
     ylab = 'layer2', pch = 19, col = c('hotpink', 'royalblue')[as.numeric(dat_sep$y)])

svm = svm(y ~., data = dat_sep, kernel = 'radial', cost = 10, scale = FALSE)

plot(svm, dat_sep, 
     col = c('hotpink', 'royalblue'),
     main = 'Linearly separable Healthy vs CSR', 
     pch = 19)


## Healthy vs DR
first_two = retina_std_means[,c(1,6)]
lab = factor(retina_df_means$label)

indices_sep = which(lab == 'Healthy' | lab == 'DR')
dat_sep = data.frame(x = first_two[indices_sep, 1:2],
                     y = factor(lab[indices_sep]))

head(dat_sep)

# Jitter ? non so se serve
set.seed(2026)
dat_sep[,1:2] <- dat_sep[,1:2] + cbind(rnorm(nrow(dat_sep), sd=0.025))

plot(dat_sep[,1:2], main = 'first two layers (Healthy vs DR)', xlab = 'layer1',
     ylab = 'layer2', pch = 19, col = c('hotpink', 'royalblue')[as.numeric(dat_sep$y)])

svm = svm(y ~., data = dat_sep, kernel = 'radial', cost = 50, scale = FALSE)

plot(svm, dat_sep, 
     col = c('hotpink', 'royalblue'),
     main = 'Linearly separable Healthy vs AMD', 
     pch = 19)

## Healthy vs MH
first_two = retina_std_means[,c(1,6)]
lab = factor(retina_df_means$label)

indices_sep = which(lab == 'Healthy' | lab == 'MH')
dat_sep = data.frame(x = first_two[indices_sep, 1:2],
                     y = factor(lab[indices_sep]))

head(dat_sep)

# Jitter ? non so se serve
set.seed(2026)
dat_sep[,1:2] <- dat_sep[,1:2] + cbind(rnorm(nrow(dat_sep), sd=0.025))

plot(dat_sep[,1:2], main = 'first two layers (Healthy vs MH)', xlab = 'layer1',
     ylab = 'layer6', pch = 19, col = c('hotpink', 'royalblue')[as.numeric(dat_sep$y)])

svm = svm(y ~., data = dat_sep, kernel = 'radial', cost = 100, scale = FALSE)

plot(svm, dat_sep, 
     col = c('hotpink', 'royalblue'),
     main = 'Linearly separable Healthy vs AMD', 
     pch = 19)





