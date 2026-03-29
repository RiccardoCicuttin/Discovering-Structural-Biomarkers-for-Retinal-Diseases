install.packages("GGally")
library(ggplot2)
library(GGally)



ggpairs(data, 
        columns = 1:6, 
        aes(color = label, alpha = 0.5),
        # 1. Definisce i nomi personalizzati
        columnLabels = c("ILM-NFL", "NFL-IPL", "IPL-INL", "INL-OPL", "OPL-ONL", "ISOS-RPE"),
        # 2. Sposta le etichette sulla diagonale
        axisLabels = "internal",
        # 3. Rimuove le correlazioni in alto
        upper = list(continuous = "blank"),
        # 4. Rimuove i grafici a creste (densità) dalla diagonale
        diag = list(continuous = "blankDiag")
) + 
  theme_minimal() +
  # Opzionale: pulizia estetica per le etichette diagonali
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(), # Rimuove le etichette esterne ormai duplicate
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  ) +
  labs(title = "Matrice di correlazione degli spessori retinici",
       x = "Spessore Misurato (pixel)",
       y = "Spessore Misurato (pixel)")


  
