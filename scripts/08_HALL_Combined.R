# ==============================================================================
# Script Name: 08_HALL_Combined.R
# Description: Harmonizes and visualizes GSEA Hallmark pathway enrichment shifts 
#              across distinct biological conditions (Aging: Old vs Young, 
#              Intervention: Exercise/Sedentary vs Old, and Pharmacological: 
#              GSK1016790A vs Vehicle). Generates high-dimensional dotplot matrices
#              displaying normalized enrichment scores (NES) and statistical 
#              significance (-log10 FDR) across osteogenic or neural tissues.
# Input Data:  MSigDB GSEA export tables (*.csv or data frames):
#              - 'gsea_report_for_na_pos_*' / 'gsea_report_for_na_neg_*' 
#                corresponding to pairwise comparison outputs (O vs Y, S vs O, GSK vs Veh).
# ==============================================================================
 
# 0. Load Required Libraries 
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
 
# 1. Data Ingestion and Multi-Cohort Integration 
# Bind directional GSEA summary reports across comparison cohorts
# Group annotations:
#   - "O vs Y": Old (16M) vs Young (3M) physiological aging baseline
#   - "S vs O": Treadmill Running (Exercise) vs Sedentary Old control
#   - "GSK vs Veh": TRPV4 agonist intervention vs Vehicle-treated control

df_oy_up    <- gsea_report_for_na_pos_1720664519716 %>% mutate(Comparison = "O vs Y")
df_oy_down  <- gsea_report_for_na_neg_1720664519716 %>% mutate(Comparison = "O vs Y")
df_so_up    <- gsea_report_for_na_pos_1720665223327 %>% mutate(Comparison = "S vs O")
df_so_down  <- gsea_report_for_na_neg_1720665223327 %>% mutate(Comparison = "S vs O")
df_gsk_up   <- gsea_report_for_na_pos_1782893335594 %>% mutate(Comparison = "GSK vs Veh")
df_gsk_down <- gsea_report_for_na_neg_1782893335594 %>% mutate(Comparison = "GSK vs Veh")

# Concatenate all directional comparison subsets into a unified master frame
gsea_combined <- bind_rows(df_oy_up, df_oy_down, df_so_up, df_so_down, df_gsk_up, df_gsk_down)
 
# 2. Data Curation and Factor Level Serialization 
plot_data <- gsea_combined %>%
  # Strip redundant prefixes to enhance gene set readability on Y-axis
  mutate(Pathway = str_replace(NAME, "HALLMARK_", "")) %>%
# Calculate -log10 FDR q-value for geometric point size mapping (pseudo-count 1e-5 avoids Inf)
mutate(Neg_Log_FDR = -log10(FDR.q.val + 1e-5)) %>%
# Enforce biological and intervention progression across X-axis
mutate(Comparison = factor(Comparison, levels = c("O vs Y", "S vs O", "GSK vs Veh")))

# Dynamically rank pathways along Y-axis based on their NES in physiological aging (O vs Y)
pathway_order <- plot_data %>% 
  filter(Comparison == "O vs Y") %>% 
  arrange(NES) %>% 
  pull(Pathway)

plot_data <- plot_data %>%
  mutate(Pathway = factor(Pathway, levels = pathway_order))
 
# 3. High-Dimensional GSEA Bubble Matrix Visualization 
p_combined <- ggplot(plot_data, aes(x = Comparison, y = Pathway)) +
  # Map geometric point size to statistical confidence and fill color to directional enrichment (NES)
  geom_point(aes(size = Neg_Log_FDR, color = NES)) +
# Classic diverging palette: Red (Activated/Upregulated), Blue (Suppressed/Downregulated)
scale_color_gradient2(
  low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0,
  name = "NES"
) +
scale_size_continuous(
  name = expression(-log[10](FDR)),
  range = c(2, 8)
) +
theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
labs(title = "HALLMARK pathways shift across interventions")[cite: 3, 5]

# Render plot to graphic device
print(p_combined)

# Save publication-grade vector PDF
# ggsave("Fig_Bone_GSEA_Combined_Shift.pdf", plot = p_combined, width = 7, height = 10, dpi = 300)
 
# 4. Source Data Export (Nature Portfolio Compliance) 
target_dir <- "/users/kenny/Desktop/2025/运动/sub/re-sub/Table"

df_source_gsea <- plot_data %>%
  dplyr::select(Pathway, Comparison, NES, FDR.q.val, Neg_Log_FDR)

write.csv(df_source_gsea, 
          file = file.path(target_dir, "Source_Data_GSEA_Hallmark_Shifts.csv"), 
          row.names = FALSE)