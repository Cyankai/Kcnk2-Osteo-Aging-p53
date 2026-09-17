# ==============================================================================
# Script Name: 08_HALL_Combined.R
# Description: Integrates cross-species transcriptomic data (human and murine 
#              exercise cohorts). Identifies evolutionarily conserved translational 
#              targets by intersecting differentially expressed genes (DEGs) 
#              and visualizing the concordance of bone-derived factors.
# Input Data:  'data_GSK', 'data_S' (Murine DEGs), and 'data_H' (Human DEGs).
# Output:      Four-quadrant scatter plot highlighting conserved targeted pathways.
# ==============================================================================
## 四象限图

# 加载必要的包
library(dplyr)
library(tibble)

# 1. 处理小鼠 GSK 数据 (提取行名，转大写，重命名列)
df_GSK <- data_GSK %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene_Symbol") %>%
  mutate(Gene_Symbol = toupper(Gene_Symbol)) %>% # 转换为全大写
  select(Gene_Symbol, 
         Log2FC_Mouse_GSK = logFC, 
         Pval_Mouse_GSK = PValue)

# 2. 处理小鼠运动数据 data_S (提取行名，转大写，重命名列)
df_S <- data_S %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene_Symbol") %>%
  mutate(Gene_Symbol = toupper(Gene_Symbol)) %>% 
  select(Gene_Symbol, 
         Log2FC_Mouse_Ex = logFC, 
         Pval_Mouse_Ex = PValue)

# 3. 处理人类运动数据 data_H (列 X 改名，确保大写，重命名列)
df_H <- data_H %>%
  as.data.frame() %>%
  rename(Gene_Symbol = X) %>%
  mutate(Gene_Symbol = toupper(Gene_Symbol)) %>%
  select(Gene_Symbol, 
         Log2FC_Human_Ex = logFC, 
         Pval_Human_Ex = PValue)

# 4. 将三个数据框完美合并 (以 Gene_Symbol 为桥梁)
# 推荐使用 inner_join (只保留在三个数据集中共同检测到的基因)
# 如果你想保留所有基因，就把 inner_join 改成 full_join
merged_data <- df_GSK %>%
  inner_join(df_S, by = "Gene_Symbol") %>%
  inner_join(df_H, by = "Gene_Symbol")

# 查看合并后的结果
head(merged_data)

library(ggplot2)
library(ggrepel) # 用于给基因打标签防重叠

# 设定显著性阈值，标记哪些是重要的靶点 (例如 P < 0.05 且 都有下降)
plot_df <- merged_data %>%
  mutate(
    # 找出在小鼠和人类运动后都显著下调的核心靶点 (第三象限)
    Target_Class = case_when(
      Log2FC_Mouse_Ex < 0 & Pval_Mouse_Ex < 0.05 & Log2FC_Human_Ex < 0 & Pval_Human_Ex < 0.05 ~ "Conserved Down",
      Log2FC_Mouse_Ex > 0 & Pval_Mouse_Ex < 0.05 & Log2FC_Human_Ex > 0 & Pval_Human_Ex < 0.05 ~ "Conserved Up",
      TRUE ~ "Not Significant"
    ),
    # 只给核心靶点加上文字标签，比如 DKK1, SOST 等
    Label = ifelse(Target_Class != "Not Significant", Gene_Symbol, NA)
  )

# 绘制四象限图
ggplot(plot_df, aes(x = Log2FC_Mouse_Ex, y = Log2FC_Human_Ex, color = Target_Class)) +
  geom_point(alpha = 0.7, size = 2) +
  # 添加十字参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # 添加基因名称
  geom_text_repel(aes(label = Label), size = 4, fontface = "bold.italic", max.overlaps = 20) +
  scale_color_manual(values = c("Conserved Down" = "#1F77B4", "Conserved Up" = "#D62728", "Not Significant" = "grey80")) +
  theme_classic() +
  labs(x = "Log2 Fold Change (Mouse Exercise)",
       y = "Log2 Fold Change (Human Exercise)",
       title = "Translational Concordance of Bone-Derived Factors")
