source("R/packages.R", echo=FALSE)
source("R/feature_utils.R")

path <- "data/PatientsData.xlsx"

# read raw data
read_sheet_raw <- function(path, sheet, label) {
  read_excel(path, sheet = sheet) %>%
    mutate(
      label      = label,
      patient_id = as.integer(str_extract(image, "(?<=[A-Za-z_])\\d+(?=\\.jpe?g$)"))
    ) %>%
    dplyr::select(-image) %>%
    relocate(patient_id, label, .after = last_col())
}

df_raw <- map2_dfr(
  list(1,       2,     3,    4,     5),
  c("AMD",   "CSR",  "DR", "MH", "Healthy"),
  ~ read_sheet_raw(path, .x, .y)
)

# Save the mapping from clean names (layer1...layer6) to original names
layer_order    <- df_raw %>% filter(label == "AMD") %>% pull(names) %>% unique()
layer_name_map <- setNames(layer_order, paste0("layer", 1:6))
# layer_name_map["layer1"] 
saveRDS(layer_name_map, "results/layer_name_map.rds")

# recode the "names" column using the inverted map
name_recode <- setNames(paste0("layer", 1:6), layer_order)
df_raw <- df_raw %>%
  mutate(names = name_recode[names]) 

# each patient is uniquely identified by (label, id)
patient_key <- df_raw %>% distinct(label, patient_id) %>% mutate(row_idx = row_number())

# create partition train - test based on the keys of patients
set.seed(2026)
train_idx      <- createDataPartition(patient_key$label, p = 0.75, list = FALSE)
train_patients <- patient_key[ train_idx, ]
test_patients  <- patient_key[-train_idx, ]

df_raw_train <- df_raw %>% semi_join(train_patients, by = c("label", "patient_id"))
df_raw_test  <- df_raw %>% semi_join(test_patients,  by = c("label", "patient_id"))

saveRDS(train_patients, "datasets/train_patients.rds")
saveRDS(test_patients,  "datasets/test_patients.rds")
saveRDS(df_raw_train,   "datasets/df_raw_train.rds")
saveRDS(df_raw_test,    "datasets/df_raw_test.rds")
saveRDS(df_raw,         "datasets/df_raw.rds")

# compute_features lives in R/feature_utils.R
train_features <- compute_features(df_raw_train)
test_features  <- compute_features(df_raw_test)

# sd and mean of each layer are correlated
# coefficient of variation shares less information
het_comparison <- map_dfr(1:6, function(k) {
  mn  <- train_features[[paste0("layer",    k)]]
  sd  <- train_features[[paste0("sd_layer", k)]]
  cv  <- sd / mn
  lcv <- log(cv)
  
  tibble(
    layer          = paste0("layer", k),
    cor_mean_sd    = cor(mn, sd,  use = "complete.obs"),  
    cor_mean_CV    = cor(mn, cv,  use = "complete.obs"),  
    cor_mean_logCV = cor(mn, lcv, use = "complete.obs"),  
    skew_CV        = skewness(cv,  na.rm = TRUE),         
    skew_logCV     = skewness(lcv, na.rm = TRUE)      
  )
})

print(het_comparison, digits = 3)

train_logCV <- add_logCV(train_features)
test_logCV  <- add_logCV(test_features)

mean_cols  <- paste0("layer",       1:6)
sd_cols    <- paste0("sd_layer",    1:6)
ratio_cols <- paste0("ratio_l",     c("1l2","2l3","3l4","4l5","5l6"))
logCV_cols <- paste0("logCV_layer", 1:6)

train_datasets <- list(
  means              = as_x(train_features, mean_cols),
  means_sd           = as_x(train_features, c(mean_cols, sd_cols)),
  means_sd_ratios    = as_x(train_features, c(mean_cols, sd_cols, ratio_cols)),
  means_logCV_ratios = as_x(train_logCV,    c(mean_cols, logCV_cols, ratio_cols))
)

test_datasets <- list(
  means              = as_x(test_features, mean_cols),
  means_sd           = as_x(test_features, c(mean_cols, sd_cols)),
  means_sd_ratios    = as_x(test_features, c(mean_cols, sd_cols, ratio_cols)),
  means_logCV_ratios = as_x(test_logCV,    c(mean_cols, logCV_cols, ratio_cols))
)

y_train <- factor(train_features$label)
y_test  <- factor(test_features$label)


# datasets/ is local only (inserted in .gitignore) 
saveRDS(train_datasets,      "datasets/train_datasets.rds")
saveRDS(test_datasets, "datasets/test_datasets.rds")
saveRDS(y_train,       "datasets/y_train.rds")
saveRDS(y_test,        "datasets/y_test.rds")
saveRDS(train_features, "datasets/train_features.rds")
saveRDS(test_features,  "datasets/test_features.rds")
saveRDS(train_logCV,    "datasets/train_logCV.rds")
saveRDS(test_logCV,     "datasets/test_logCV.rds")



