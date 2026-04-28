# 3 different sets of features
# retina_df_means
retina_df_mns_sd = retina_df[, c(1:12, 18)]
# retina_df

library(caret)
library(nestedcv)
library(MLmetrics)

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





