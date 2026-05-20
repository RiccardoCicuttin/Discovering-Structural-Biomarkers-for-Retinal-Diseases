library(readxl)
library(dplyr) 
library(tidyr)
library(purrr)
library(stringr)
library(caret)

path = "data/PatientsData.xlsx"

data_AMD = read_excel(path, sheet = 1) %>% 
  # mean of every row
  mutate(patient_id = as.integer(str_extract(image, "(?<=[A-Za-z_])\\d+(?=\\.jpe?g$)")),
         mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd),) %>%
  # select just the mean and the sd value, the name of the layer, and the image ID
  select(mean_value, sd_val, names, image, patient_id) %>%
  # use as observation ID the image ID to change the format
  pivot_wider(names_from = names, 
              values_from = c(mean_value, sd_val), 
              id_cols = c(image, patient_id)) %>%
  # erase the image ID and add the disease label column
  select(-image) %>%
  mutate(
    ratio_l1l2 = .[[1]] / .[[2]],
    ratio_l2l3 = .[[2]] / .[[3]],
    ratio_l3l4 = .[[3]] / .[[4]],
    ratio_l4l5 = .[[4]] / .[[5]],
    ratio_l5l6 = .[[5]] / .[[6]]
  ) %>% mutate(label = "AMD")

data_CSR = read_excel(path, sheet = 2) %>% 
  mutate(patient_id = as.integer(str_extract(image, "(?<=[A-Za-z_])\\d+(?=\\.jpe?g$)")), 
         mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  select(mean_value, sd_val, names, image, patient_id) %>% 
  pivot_wider(names_from = names, values_from = c(mean_value, sd_val), id_cols = c(image, patient_id)) %>%
  select(-image) %>% 
  mutate(
    ratio_l1l2 = .[[1]] / .[[2]],
    ratio_l2l3 = .[[2]] / .[[3]],
    ratio_l3l4 = .[[3]] / .[[4]],
    ratio_l4l5 = .[[4]] / .[[5]],
    ratio_l5l6 = .[[5]] / .[[6]]
  ) %>% mutate(label = "CSR") 

data_DR = read_excel(path, sheet = 3) %>% 
  mutate(patient_id = as.integer(str_extract(image, "(?<=[A-Za-z_])\\d+(?=\\.jpe?g$)")),
         mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  select(mean_value, sd_val, names, image, patient_id) %>% 
  pivot_wider(names_from = names, values_from = c(mean_value, sd_val), id_cols = c(image, patient_id)) %>% 
  select(-image) %>%
  mutate(
    ratio_l1l2 = .[[1]] / .[[2]],
    ratio_l2l3 = .[[2]] / .[[3]],
    ratio_l3l4 = .[[3]] / .[[4]],
    ratio_l4l5 = .[[4]] / .[[5]],
    ratio_l5l6 = .[[5]] / .[[6]]
  ) %>% mutate(label = "DR")

data_MH = read_excel(path, sheet = 4) %>% 
  mutate(patient_id = as.integer(str_extract(image, "(?<=[A-Za-z_])\\d+(?=\\.jpe?g$)")),
         mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  select(mean_value, sd_val, names, image, patient_id) %>% 
  pivot_wider(names_from = names, values_from = c(mean_value, sd_val), id_cols = c(image, patient_id)) %>%
  select(-image) %>% 
  mutate(
    ratio_l1l2 = .[[1]] / .[[2]],
    ratio_l2l3 = .[[2]] / .[[3]],
    ratio_l3l4 = .[[3]] / .[[4]],
    ratio_l4l5 = .[[4]] / .[[5]],
    ratio_l5l6 = .[[5]] / .[[6]]
  ) %>% mutate(label = "MH")

data_healthy = read_excel(path, sheet = 5) %>% 
  mutate(patient_id = as.integer(str_extract(image, "(?<=[A-Za-z_])\\d+(?=\\.jpe?g$)")),
         mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  select(mean_value, sd_val, names, image, patient_id) %>% 
  pivot_wider(names_from = names, values_from = c(mean_value, sd_val), id_cols = c(image, patient_id)) %>% 
  select(-image) %>%
  mutate(
    ratio_l1l2 = .[[1]] / .[[2]],
    ratio_l2l3 = .[[2]] / .[[3]],
    ratio_l3l4 = .[[3]] / .[[4]],
    ratio_l4l5 = .[[4]] / .[[5]],
    ratio_l5l6 = .[[5]] / .[[6]]
  ) %>% mutate(label = "Healthy")

# merge into one data set the data sets obtained
retina_df = bind_rows(data_healthy, data_AMD, data_CSR, data_DR, data_MH)

# Original names are not self-explicative and presenten spaces in them
# names(data)[1:6]
# We saved them using a map, with keys the new names
# New names are "layer_j" with j=1,..,6. 
# layer1 is the most internal one
original_names = colnames(retina_df[,-c(1, 13:18)])
new_names = c(paste0("layer", 1:6), paste0("sd_layer", 1:6))
retina_names_map = setNames(original_names, new_names)
#retina_names_map['layer1']

names(retina_df)[na.omit(1:12)] = new_names

# NAs are not present
sum(is.na(retina_df))==0 #TRUE


# Data in original format, with patient ID 
df_orig = bind_rows(read_excel(path, sheet = 1),
                    read_excel(path, sheet = 2),
                    read_excel(path, sheet = 3),
                    read_excel(path, sheet = 4),
                    read_excel(path, sheet = 5))
sheet_indices <- 1:5
labels <- c("Healthy", "MH", "DR", "CSR", "AMD")  

df_orig <- map2_dfr(sheet_indices, labels, function(idx, lbl) {
  read_excel(path, sheet = idx) %>%
    mutate(
      label      = lbl,
      patient_id = as.integer(str_extract(image, "(?<=[A-Za-z_])\\d+(?=\\.jpe?g$)"))
    ) %>%
    select(-image)
})

# --- Build a patient-level key table ----
# One row per patient: just the (label, patient_id) pair.
# This is the unit we stratify on.
patient_key <- retina_df %>%
  distinct(label, patient_id)

# --- Stratified split on the patient keys ----
# createDataPartition stratifies on 'label', so class proportions
# are preserved in both splits.  p = 0.75 means ~75 % of patients
# from each class go to training.
set.seed(123)
train_idx <- createDataPartition(patient_key$label, p = 0.75, list = FALSE)

train_patients <- patient_key[ train_idx, ]
test_patients  <- patient_key[-train_idx, ]

# --- Subset the TRANSFORMED (wide) dataset ----
# semi_join keeps only rows whose (label, patient_id) appears in the
# right-hand table.  It's an "existence filter": no duplication, no
# column additions — just a clean subset.
train_set <- retina_df %>% semi_join(train_patients, by = c("label", "patient_id"))
test_set  <- retina_df %>% semi_join(test_patients,  by = c("label", "patient_id"))

# --- Subset the ORIGINAL (long) dataset ----
# Exactly the same join keys.  Each patient has 6 rows in df_orig
# (one per layer), and semi_join keeps all 6 for every patient that
# landed in the corresponding split.
train_orig <- df_orig %>% semi_join(train_patients, by = c("label", "patient_id"))
test_orig  <- df_orig %>% semi_join(test_patients,  by = c("label", "patient_id"))




retina_df_logCV <- retina_df %>%
  mutate(across(
    all_of(paste0("sd_layer", 1:6)),
    ~ log(.x / retina_df[[paste0("layer", str_extract(cur_column(), "\\d+"))]]),
    .names = "logCV_layer{str_extract(.col, '\\\\d+')}"
  )) %>%
  select(all_of(paste0("layer", 1:6)),
         starts_with("logCV_layer"),
         starts_with("ratio_"),
         label)







