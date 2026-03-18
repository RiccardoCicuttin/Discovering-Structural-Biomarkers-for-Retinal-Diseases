# 1. Caricamento delle librerie necessarie
library(readxl)
library(dplyr) # per manipolare i dati

# 2. Definizione del percorso relativo
path <- "data/PatientsData.xlsx"

# 3. Lettura dei fogli Excel

data_AMD <- read_excel(path, sheet = 1) 

data_CSR <- read_excel(path, sheet = 2)

data_DR <- read_excel(path, sheet = 3)

data_MH <- read_excel(path, sheet = 4)

data_normal <- read_excel(path, sheet = 5)

# Verifica visiva dei dati caricati
head(data_AMD)

