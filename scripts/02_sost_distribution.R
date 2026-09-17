# ==============================================================================
# Script Name: 02_sost_distribution.R
# Description: Performs spatial transcriptomic/proteomic mapping of Sclerostin 
#              expression across the cortical bone. It generates combined 
#              scatter and raincloud plots to quantify the endosteal accumulation 
#              (Inner zone) of Sclerostin in aged mice and its reversal by exercise.
# Input Data:  'clean_data' (Dataframe containing Group, Depth, and Fluo intensity).
# Output:      Combined spatial distribution plots (Scatter + Raincloud).
# ============================================================================== 

# 1. 加载必要的包
library(ggplot2)
library(dplyr)
library(gghalves)
library(ggpubr)
library(patchwork)

# 全局经典配色 
group_colors <- c("Young" = "#4DBBD5", "Old" = "#E64B35", "Exercise" = "#00A087")

# ==============================================================================
# 2. 统一数据源与动态划分区域 (基于 clean_data 自动生成所有指标)
# ==============================================================================
# 假设您的 clean_data 包含: Group, Depth, Fluo

filtered_data <- clean_data %>%
  # 剔除无效深度
  filter(Depth >= 0 & Depth <= 100) %>%
  # 核心：根据 Depth 自动划分 Spatial Zone
  mutate(
    Zone = case_when(
      Depth <= 100/3 ~ "Inner",
      Depth > 100/3 & Depth <= 200/3 ~ "Mid",
      Depth > 200/3 ~ "Outer"
    ),
    # 锁定因子顺序，确保两张图的组别和分面顺序完全统一
    Group = factor(Group, levels = c("Young", "Old", "Exercise")),
    Zone = factor(Zone, levels = c("Inner", "Mid", "Outer")) # 保留您要求的雨云图从深到浅的排版
  )

# ==============================================================================
# 3. 动态计算散点图顶部的百分比标注
# ==============================================================================
anno_df <- filtered_data %>%
  group_by(Group, Zone) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    Total = sum(Count),
    Percentage = Count / Total,
    Label = paste0(round(Percentage * 100, 1), "%"),
    # 将文字精准锚定在各区域 X 轴的中点
    Depth = case_when(
      Zone == "Inner" ~ 16.5,
      Zone == "Mid" ~ 50,
      Zone == "Outer" ~ 83.5
    ),
    Fluo = Inf # 让文字永远悬浮在最顶部
  )

# ==============================================================================
# 4. 绘制左侧：自带门控与动态定量标注的空间散点图
# ==============================================================================
p_scatter_final <- ggplot(filtered_data, aes(x = Depth, y = Fluo)) +
  
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "gray30", linewidth = 0.8) +
  geom_vline(xintercept = c(33, 66), linetype = "dashed", color = "gray60", linewidth = 0.8) +
  
  geom_point(aes(color = Group), alpha = 0.7, size = 1.8, shape = 16) +
  
  # 自动调用算好的百分比
  geom_text(data = anno_df, aes(label = Label), y = Inf, vjust = 1.8, 
            size = 4.5, fontface = "bold", color = "black") +
  
  facet_wrap(~ Group, ncol = 3) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = c(0, 33, 66, 100)) +
  
  # 统一纵坐标名称，去掉 A.U. 使图面更干净
  labs(x = "Relative Depth (%)", 
       y = "Relative Sclerostin intensity") +
  
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.x = element_text(color = "black", face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(color = "black", margin = margin(r = 10)),
    axis.line = element_line(linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# ==============================================================================
# 5. 绘制右侧：基于同一数据的空间分布雨云图
# ==============================================================================
# 动态计算每个提琴图底部的 n 值
cell_counts <- filtered_data %>%
  group_by(Zone, Group) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(y_pos = min(filtered_data$Fluo) - 0.15)

my_comparisons <- list(c("Young", "Old"), c("Old", "Exercise"))

p_rain <- ggplot(filtered_data, aes(x = Group, y = Fluo, fill = Group, color = Group)) +
  
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "gray30", linewidth = 0.8) +
  
  geom_half_violin(side = "r", position = position_nudge(x = 0.15, y = 0),
                   alpha = 0.8, color = NA, trim = FALSE) +
  geom_half_point(side = "l", position = position_nudge(x = -0.1, y = 0),
                  size = 0.4, alpha = 0.6, transformation = position_jitter(width = 0.05)) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black",
               position = position_nudge(x = 0.05, y = 0)) +
  
  facet_wrap(~ Zone) +
  
  geom_text(data = cell_counts, aes(x = Group, y = y_pos, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 4, fontface = "bold", color = "black") +
  
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     p.adjust.method = "BH", label = "p.signif", vjust = 0.5, size = 5) +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  labs(y = "Relative Sclerostin intensity", x = NULL) +
  
  theme_classic(base_size = 15) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(color = "black", face = "bold", size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", margin = margin(r = 10)),
    axis.line = element_line(linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  coord_cartesian(clip = "off", ylim = c(min(filtered_data$Fluo)-0.2, max(filtered_data$Fluo)+0.4)) 

# ==============================================================================
# 6. 终极无缝拼图渲染
# ==============================================================================
final_combined_plot <- p_scatter_final + p_rain + 
  plot_layout(widths = c(1, 1.1)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20, face = "bold"))

print(final_combined_plot)
