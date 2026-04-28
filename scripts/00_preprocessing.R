library(readxl)
library(dplyr) 
library(tidyr)
library(purrr)
library(stringr)

path = "data/PatientsData.xlsx"

data_AMD = read_excel(path, sheet = 1) %>% 
  # mean of every row
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  # select just the mean and the sd value, the name of the layer, and the image ID
  select(mean_value, sd_val, names, image) %>%
  # use as observation ID the image ID to change the format
  pivot_wider(names_from = names, 
              values_from = c(mean_value, sd_val), 
              id_cols = last_col()) %>%
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
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  select(mean_value, sd_val, names, image) %>% 
  pivot_wider(names_from = names, values_from = c(mean_value, sd_val), id_cols = last_col()) %>%
  select(-image) %>% 
  mutate(
    ratio_l1l2 = .[[1]] / .[[2]],
    ratio_l2l3 = .[[2]] / .[[3]],
    ratio_l3l4 = .[[3]] / .[[4]],
    ratio_l4l5 = .[[4]] / .[[5]],
    ratio_l5l6 = .[[5]] / .[[6]]
  ) %>% mutate(label = "CSR") 

data_DR = read_excel(path, sheet = 3) %>% 
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  select(mean_value, sd_val, names, image) %>% 
  pivot_wider(names_from = names, values_from = c(mean_value, sd_val), id_cols = last_col()) %>% 
  select(-image) %>%
  mutate(
    ratio_l1l2 = .[[1]] / .[[2]],
    ratio_l2l3 = .[[2]] / .[[3]],
    ratio_l3l4 = .[[3]] / .[[4]],
    ratio_l4l5 = .[[4]] / .[[5]],
    ratio_l5l6 = .[[5]] / .[[6]]
  ) %>% mutate(label = "DR")

data_MH = read_excel(path, sheet = 4) %>% 
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  select(mean_value, sd_val, names, image) %>% 
  pivot_wider(names_from = names, values_from = c(mean_value, sd_val), id_cols = last_col()) %>%
  select(-image) %>% 
  mutate(
    ratio_l1l2 = .[[1]] / .[[2]],
    ratio_l2l3 = .[[2]] / .[[3]],
    ratio_l3l4 = .[[3]] / .[[4]],
    ratio_l4l5 = .[[4]] / .[[5]],
    ratio_l5l6 = .[[5]] / .[[6]]
  ) %>% mutate(label = "MH")

data_healthy = read_excel(path, sheet = 5) %>% 
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE), 
         sd_val = apply(select(., where(is.numeric)), na.rm=TRUE,1,sd)) %>%
  select(mean_value, sd_val, names, image) %>% 
  pivot_wider(names_from = names, values_from = c(mean_value, sd_val), id_cols = last_col()) %>% 
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
original_names = colnames(retina_df[,-c(13:18)])
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
labels <- c("normal", "MH", "DR", "CSR", "AMD")  

df_orig <- map2_dfr(sheet_indices, labels, function(idx, lbl) {
  read_excel(path, sheet = idx) %>%
    mutate(
      label      = lbl,
      patient_id = as.integer(str_extract(image, "(?<=[A-Za-z_])\\d+(?=\\.jpe?g$)"))
    ) %>%
    select(-image)
})









