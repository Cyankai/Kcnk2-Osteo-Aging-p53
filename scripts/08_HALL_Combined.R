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
#   - bo : Bone
#   - br : Brain
#   - "O vs Y": Old (16M) vs Young (3M) physiological aging baseline
#   - "S vs O": Treadmill Running (Exercise) vs Sedentary Old control
#   - "GSK vs Veh": TRPV4 agonist intervention vs Vehicle-treated control

df_br_oy_up    <- gsea_report_for_na_pos_1720664519716 %>% mutate(Comparison = "O vs Y")
df_br_oy_down  <- gsea_report_for_na_neg_1720664519716 %>% mutate(Comparison = "O vs Y")
df_br_so_up    <- gsea_report_for_na_pos_1720665223327 %>% mutate(Comparison = "S vs O")
df_br_so_down  <- gsea_report_for_na_neg_1720665223327 %>% mutate(Comparison = "S vs O")
df_br_gsk_up   <- gsea_report_for_na_pos_1782893335594 %>% mutate(Comparison = "GSK vs Veh")
df_br_gsk_down <- gsea_report_for_na_neg_1782893335594 %>% mutate(Comparison = "GSK vs Veh")

df_bo_oy_up    <- gsea_report_for_na_pos_1716970278396 %>% mutate(Comparison = "O vs Y")
df_bo_oy_down  <- gsea_report_for_na_neg_1716970278396 %>% mutate(Comparison = "O vs Y")
df_bo_so_up    <- gsea_report_for_na_pos_1721372323907 %>% mutate(Comparison = "S vs O")
df_bo_so_down  <- gsea_report_for_na_neg_1721372323907 %>% mutate(Comparison = "S vs O")
df_bo_gsk_up   <- gsea_report_for_na_pos_1782893508589 %>% mutate(Comparison = "GSK vs Veh")
df_bo_gsk_down <- gsea_report_for_na_neg_1782893508589 %>% mutate(Comparison = "GSK vs Veh")

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(fmsb) # 必须加载 fmsb 包绘制雷达图

# ==============================================================================
# 1. 数据清洗与合并模块 (区分骨与脑)
# ==============================================================================

# 合并脑组织数据
plot_data_brain <- bind_rows(
  df_br_oy_up, df_br_oy_down, 
  df_br_so_up, df_br_so_down, 
  df_br_gsk_up, df_br_gsk_down
) %>%
  mutate(Pathway = str_replace(NAME, "HALLMARK_", "")) %>%
  mutate(Comparison = factor(Comparison, levels = c("O vs Y", "S vs O", "GSK vs Veh")))

# 合并骨组织数据
plot_data_bone <- bind_rows(
  df_bo_oy_up, df_bo_oy_down, 
  df_bo_so_up, df_bo_so_down, 
  df_bo_gsk_up, df_bo_gsk_down
) %>%
  mutate(Pathway = str_replace(NAME, "HALLMARK_", "")) %>%
  mutate(Comparison = factor(Comparison, levels = c("O vs Y", "S vs O", "GSK vs Veh")))


# ==============================================================================
# 2. 定义核心通路与通用绘图函数
# ==============================================================================

# 14条核心通路及美化标签
core_hallmarks <- c(
  "OXIDATIVE_PHOSPHORYLATION", "FATTY_ACID_METABOLISM", "GLYCOLYSIS",
  "EPITHELIAL_MESENCHYMAL_TRANSITION", "WNT_BETA_CATENIN_SIGNALING", "ANGIOGENESIS",
  "IL6_JAK_STAT3_SIGNALING", "INFLAMMATORY_RESPONSE", "COMPLEMENT", "INTERFERON_GAMMA_RESPONSE",
  "P53_PATHWAY", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "E2F_TARGETS", "G2M_CHECKPOINT"
)
short_labels <- c(
  "OxPhos", "FA Metab", "Glycolysis",
  "EMT", "Wnt/β-Cat", "Angiogenesis",
  "IL-6/STAT3", "Inflammation", "Complement", "IFN-γ",
  "p53 Axis", "ROS", "E2F Targets", "G2M Check"
)

# 封装自动化绘图函数
generate_tissue_plots <- function(data, tissue_name, color_theme) {
  
  # --- 通用数据宽格式准备 ---
  df_wide <- data %>%
    dplyr::select(Pathway, Comparison, NES) %>%
    tidyr::pivot_wider(names_from = Comparison, values_from = NES) %>%
    dplyr::rename(Exercise = `S vs O`, GSK = `GSK vs Veh`, Aging = `O vs Y`)
  
  # ----------------------------------------------------------------------------
  # 图 A：14 条核心通路雷达图 (Radar Chart)
  # ----------------------------------------------------------------------------
  # 提取雷达图所需的格式
  radar_prep <- data %>%
    filter(Pathway %in% core_hallmarks) %>%
    filter(Comparison %in% c("S vs O", "GSK vs Veh")) %>%
    dplyr::select(Pathway, Comparison, NES) %>%
    tidyr::pivot_wider(names_from = Pathway, values_from = NES) %>%
    as.data.frame()
  
  rownames(radar_prep) <- radar_prep$Comparison
  radar_matrix <- radar_prep[, core_hallmarks] # 确保列顺序对应
  colnames(radar_matrix) <- short_labels
  
  # 设定极值界限
  max_nes <- 2.5
  min_nes <- -2.5
  radar_data <- rbind(rep(max_nes, ncol(radar_matrix)), rep(min_nes, ncol(radar_matrix)), radar_matrix)
  
  # 颜色配置
  colors_border <- c("#3498DB", color_theme)
  colors_fill   <- c(scales::alpha("#3498DB", 0.20), scales::alpha(color_theme, 0.20))
  
  # 绘制雷达图 (直接输出到绘图设备)
  par(mar = c(3, 3, 3, 3), xpd = TRUE)
  radarchart(
    radar_data,
    axistype = 1, pcol = colors_border, pfcol = colors_fill,
    plwd = 2.5, plty = 1, cglcol = "grey75", cglty = 2, axislabcol = "grey40",
    caxislabels = seq(min_nes, max_nes, length.out = 5), cglwd = 1.0, vlcex = 0.85,
    title = paste0(tissue_name, ": Shift of Key Aging Hallmarks")
  )
  legend(
    x = "topright", # 将图例移至右上角，避开底部密集的通路标签
    inset = c(-0.1, 0), # 微调水平和垂直位置
    legend = c("Exercise (S vs O)", "TRPV4 Agonist (GSK vs Veh)"),
    bty = "n", pch = 19, col = colors_border, text.col = "black",
    cex = 0.9, pt.cex = 1.5, horiz = FALSE # 改为垂直排列
  )
  
  # ----------------------------------------------------------------------------
  # 图 B：保真度拟合散点图 (Exercise vs GSK)
  # ----------------------------------------------------------------------------
  p_fidelity <- ggscatter(
    df_wide %>% drop_na(Exercise, GSK), 
    x = "Exercise", y = "GSK", size = 2.5, color = color_theme, alpha = 0.7,
    add = "reg.line", add.params = list(color = "black", linetype = "dashed"),
    conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 5,
    xlab = "Exercise Effect (NES)", ylab = "Pharmacological Effect (NES)",
    title = paste0(tissue_name, ": Transcriptomic Fidelity")
  ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", inherit.aes = FALSE) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # ----------------------------------------------------------------------------
  # 图 C：衰老逆转拟合散点图 (Aging vs GSK)
  # ----------------------------------------------------------------------------
  p_rescue <- ggscatter(
    df_wide %>% drop_na(Aging, GSK), 
    x = "Aging", y = "GSK", size = 2.5, color = "#2CA02C", alpha = 0.7,
    add = "reg.line", add.params = list(color = "black", linetype = "dashed"),
    conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 5,
    xlab = "Aging Effect (NES: O vs Y)", ylab = "Pharmacological Effect (NES)",
    title = paste0(tissue_name, ": Reversal of Aging")
  ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", inherit.aes = FALSE) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 打印 ggplot 图表
  print(p_fidelity)
  print(p_rescue)
}


# ==============================================================================
# 3. 运行并生成图表
# ==============================================================================

# 为防止雷达图画幅过小被截断，如果你在 RStudio 中运行，建议先稍微拉大 Plots 窗口

# 生成骨组织 (Bone) 的 3 张图 (使用红色主题代表 GSK)
cat("Generating plots for Bone Tissue...\n")
generate_tissue_plots(plot_data_bone, tissue_name = "Bone Tissue", color_theme = "#E74C3C")

# 生成脑组织 (Brain) 的 3 张图 (使用紫色主题代表 GSK，以示区分)
cat("Generating plots for Brain Tissue...\n")
generate_tissue_plots(plot_data_brain, tissue_name = "Brain Tissue", color_theme = "#9B59B6")
