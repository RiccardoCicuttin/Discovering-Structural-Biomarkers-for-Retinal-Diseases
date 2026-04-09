library(readxl)
library(dplyr) 
library(tidyr)

path = "data/PatientsData.xlsx"

data_AMD = read_excel(path, sheet = 1) %>% 
  # mean of every row
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE)) %>%
  # select just the mean value, the name of the layer, and the image ID
  select(mean_value, names, image) %>% 
  # use as observation ID the image ID to change the format
  pivot_wider(names_from = names, values_from = mean_value, id_cols = last_col()) %>%
  # erase the image ID and add the disease label column
  select(-image) %>% mutate(label = "AMD")

data_CSR = read_excel(path, sheet = 2) %>% 
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE)) %>%
  select(mean_value, names, image) %>% 
  pivot_wider(names_from = names, values_from = mean_value, id_cols = last_col()) %>%
  select(-image) %>% mutate(label = "CSR")

data_DR = read_excel(path, sheet = 3) %>% 
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE)) %>%
  select(mean_value, names, image) %>% 
  pivot_wider(names_from = names, values_from = mean_value, id_cols = last_col()) %>% 
  select(-image) %>% mutate(label = "DR")

data_MH = read_excel(path, sheet = 4) %>% 
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE)) %>%
  select(mean_value, names, image) %>% 
  pivot_wider(names_from = names, values_from = mean_value, id_cols = last_col()) %>%
  select(-image) %>% mutate(label = "MH")

data_healthy = read_excel(path, sheet = 5) %>% 
  mutate(mean_value=rowMeans(select(., where(is.numeric)), na.rm=TRUE)) %>%
  select(mean_value, names, image) %>% 
  pivot_wider(names_from = names, values_from = mean_value, id_cols = last_col()) %>% 
  select(-image) %>% mutate(label = "Healthy")

# merge into one data set the data sets obtained
retina_df = bind_rows(data_healthy, data_AMD, data_CSR, data_DR, data_MH)

# Original names are not self-explicative and presenten spaces in them
# names(data)[1:6]
# We saved them using a map, with keys the new names
# New names are "layer_j" with j=1,..,6. 
# layer1 is the most internal one
original_names = colnames(retina_df[,-7])
new_names = paste0("layer", 1:6);
retina_names_map = setNames(original_names, new_names)
#retina_names_map['layer1']

names(retina_df)[na.omit(1:6)] = new_names

# NAs are not present
sum(is.na(retina_df))==0 #TRUE













