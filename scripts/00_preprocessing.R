# 1. Caricamento delle librerie necessarie
library(readxl)
library(dplyr) # per manipolare i dati

# 2. Definizione del percorso relativo
path <- "data/PatientsData.xlsx"

# 3. Lettura dei fogli Excel
# Se il primo foglio contiene il Layer 1 e il secondo il Layer 2, farai così:
data_AMD <- read_excel(path, sheet = 1) 
# oppure usando il nome del foglio: sheet = "NomeDelFoglio"

data_CSR <- read_excel(path, sheet = 2)

# Verifica visiva dei dati caricati
head(data_AMD)

