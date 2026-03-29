install.packages("GGally")
library(ggplot2)
library(GGally)


ggpairs(retina_df, 
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



ggpairs(retina_df, 
        columns = 1:6, 
        aes(color = label, alpha = 0.5), # Colore attivo per i punti in basso
        
        columnLabels = c("ILM-NFL", "NFL-IPL", "IPL-INL", "INL-OPL", "OPL-ONL", "ISOS-RPE"),
        axisLabels = "internal",
        
        # PARTE ALTA: Un solo valore numerico (correlazione globale)
        upper = list(
          continuous = wrap("cor", mapping = aes(color = NULL), size = 5)
        ),
        
        # DIAGONALE: Solo i nomi dei layer
        diag = list(continuous = "blankDiag"),
        
        # PARTE BASSA: Grafici a punti colorati per patologia
        lower = list(continuous = "points")
) + 
  theme_minimal() +
  theme(
    strip.background = element_blank(), # Rimuove i bordi delle etichette esterne
    strip.text = element_blank(),       # Rimuove il testo delle etichette esterne
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  ) +
  labs(x = "Spessore (pixel)", y = "Spessore (pixel)")
