library(ggplot2)
library(GGally)
library(patchwork)
library(dplyr)
library(corrplot)
# Scatteplots
# Assicurati che il tuo vettore my_cols sia caricato nell'ambiente
# my_cols <- c("normal" = "#4E79A7", "MH" = "#F28E2B", "DR" = "#E15759", "CSR" = "#76B7B2", "AMD" = "#59A14F")

ggpairs(
  retina_df, 
  columns = 1:6, 
  # Rimuoviamo alpha da aes() per evitare che venga trattato come variabile
  mapping = aes(color = label), 
  
  columnLabels = c("ILM-NFL", "NFL-IPL", "IPL-INL", "INL-OPL", "OPL-ONL", "ISOS-RPE"),
  axisLabels = "internal", # Forza i nomi sulla diagonale
  
  # PARTE ALTA: Un solo valore numerico (correlazione globale in nero)
  upper = list(
    continuous = wrap("cor", mapping = aes(color = NULL), color = "black", size = 5, fontface = "bold")
  ),
  
  # DIAGONALE: Lasciamo vuota la geometria così restano solo i nomi di axisLabels
  diag = list(continuous = "blankDiag"),
  
  # PARTE BASSA: Punti colorati per patologia con alpha applicato fisicamente
  lower = list(
    continuous = wrap("points", alpha = 0.5, size = 1.2)
  )
) + 
  # --- COLORI PERSONALIZZATI ---
  scale_color_manual(values = my_cols) +
  scale_fill_manual(values = my_cols) +
  
  # --- TEMA E RIQUADRI (Stile Boxplot precedente) ---
  theme_minimal() +
  theme(
    # Riquadro nero spesso attorno a ogni facet
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    # Griglie per valutare la dispersione
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey85", linewidth = 0.25),
    
    # Rimuove bordi e testi esterni ridondanti
    strip.background = element_blank(), 
    strip.text = element_blank(),       
    
    # Font coerenti
    axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  labs(x = "Spessore (pixel)", y = "Spessore (pixel)")


# 1. Definiamo le variabili e le etichette per la diagonale
vars <- c("layer1", "layer2", "layer3", "layer4", "layer5", "layer6")
nomi_diag <- c("ILM-NFL", "NFL-IPL", "IPL-INL", "INL-OPL", "OPL-ONL", "ISOS-RPE")


# Modified format for plotting purposes
data_long <- retina_df %>%
  select(label, starts_with("layer")) %>%
  pivot_longer(cols = starts_with("layer"), 
               names_to = "layer_name", 
               values_to = "thickness")

# personalize colours
my_cols <- c(
  "normal" = "dodgerblue", 
  "MH"     = "#F28E2B", 
  "DR"     = "mediumpurple", 
  "CSR"    = "#E15759", 
  "AMD"    = "seagreen3"
)

# Ridgeplots
ggplot(data_long, aes(x = thickness, y = label, fill = label)) +
  geom_density_ridges(
    scale = 1.5, 
    alpha = 1, 
    color = "black", 
    rel_min_height = 0.005,
    quantile_lines = TRUE, 
    quantiles = 2,
    vline_color = "grey40",
    vline_linetype = "dashed"
  ) +
  # Applichiamo il labeller personalizzato al facet
  facet_wrap(~ layer_name, scales = "free_x") + 
  #scale_fill_brewer(palette = "my_cols") +
  scale_fill_manual(values = my_cols) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    # Griglie interne
    panel.grid.major = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey85", linewidth = 0.25),
    # Riquadro nero rigido attorno a ogni facet (supporto geometrico)
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    # Riquadro e testo del titolo del facet identici al boxplot
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    # Aggiustamento specifico per i ridgeplot per allineare il testo Y alle linee
    axis.text.y = element_text(vjust = 0)
  ) + 
  labs(
    title = "Density estimates for the mean width of retinal layers",
    x = "Width measured in pixels",
    y = "Diagnosis" 
  )


# Boxplots
ggplot(data_long, aes(x = label, y = thickness, fill = label)) +
  geom_boxplot(color = "black", outlier.shape = NULL, outlier.fill = "white", 
               outlier.size = 1.5, fatten = 2) +
  facet_wrap(~ layer_name, ncol = 3, scales = "free_y") +
  #scale_fill_brewer(palette = "mycols") +
  scale_fill_manual(values = my_cols) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.border = element_rect(color = "black", fill = NA), 
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  ) +
  labs(title = "Retinal width per layer",
       x = "Diagnosis",
       y = "Width (pixel)")

#corrplot

par(mfrow = c(2,3))

for (m in unique(retina_df$label)) {
  
  df <- subset(retina_df, label == m)[, 1:6]
  cor_mat <- cor(df, use = "pairwise.complete.obs")
  
  corrplot(cor_mat, method = "color",addCoef.col = "black",number.cex = 1,tl.col = "#2c3e50" )
  title(m, line = 1.5, cex = 1)
}

