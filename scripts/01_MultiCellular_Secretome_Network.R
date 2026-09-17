# ==============================================================================
# Script Name: 01_MultiCellular_Secretome_Network.R
# Description: Constructs a multi-cellular secretome coordination network based 
#              on DIA proteomics data. This script evaluates topological shifts 
#              in circulating osteokines and visualizes the exercise-induced 
#              systemic clearance of bone-derived Wnt antagonists.
# Input Data:  'pro_data' (Proteomics expression matrix), 'cor_matrix' (Correlation),
#              and 'p_value_matrix' (Significance matrix).
# Output:      A topological network plot highlighting spatial and functional niches.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 加载必备包
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(dplyr)
  library(ggforce)
  library(ggnewscale)
})

# ------------------------------------------------------------------------------
# 2. 从 pro_data 中动态计算“运动引起的增减比例 (Log2FC)”
# ------------------------------------------------------------------------------
# 假设 pro_data 已在环境中，且列名 'A' 开头代表基线，'B' 开头代表运动后
cols_A <- grep("^A", colnames(pro_data), value = TRUE)
cols_B <- grep("^B", colnames(pro_data), value = TRUE)

# 计算两组的平均表达量
mean_A <- rowMeans(pro_data[, cols_A, drop = FALSE], na.rm = TRUE)
mean_B <- rowMeans(pro_data[, cols_B, drop = FALSE], na.rm = TRUE)

# 计算对数变化倍数 (Log2 Fold Change)
# 正值代表运动后升高，负值代表运动后下降
gene_log2fc <- log2(mean_A / mean_B)

# 将计算结果存入数据框备用
val_df <- data.frame(
  name = rownames(pro_data),
  Log2FC = gene_log2fc,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 3. 动态提取相关性连线 (Edges)
# ------------------------------------------------------------------------------
# 提取 25 个核心靶点列表
core_genes <- c("SOST", "DKK1", "SERPINE2", "CXCL12", "ACP5", "CST3", 
                "COL1A1", "COMP", "CLEC11A", "STC2", "COL1A2", "OMD", 
                "ECM2", "FMOD", "MYOC", "FBLN7", "PRG4", "F13A1", "TNC", 
                "MIA", "FRZB", "PRG2", "SMOC2", "MGP", "COL11A1")

# 从全局矩阵中子集化核心基因
sub_cor <- cor_matrix[core_genes, core_genes]
sub_pval <- p_value_matrix[core_genes, core_genes]

# 严格阈值提取：P < 0.01 且绝对相关系数 >= 0.5 (只取上三角防止重复连线)
edge_indices <- which(sub_pval < 0.01 & abs(sub_cor) >= 0.6 & upper.tri(sub_cor), arr.ind = TRUE)

edges_df <- data.frame(
  from = rownames(sub_cor)[edge_indices[, 1]],
  to   = colnames(sub_cor)[edge_indices[, 2]],
  r_value = sub_cor[edge_indices]
) %>%
  mutate(
    weight = abs(r_value),
    Direction = ifelse(r_value > 0, "Positive", "Negative")
  )

# ------------------------------------------------------------------------------
# 4. 组装节点信息 (Nodes) 与单细胞拓扑领地
# ------------------------------------------------------------------------------
nodes_df <- data.frame(name = core_genes) %>%
  left_join(val_df, by = "name") %>%
  mutate(
    # 处理矩阵中缺失的基因，默认FC设为0 (无颜色变化)
    Log2FC = ifelse(is.na(Log2FC), 0, Log2FC),
    # 划分单细胞谱系领地
    CellType = case_when(
      name %in% c("SOST", "DKK1") ~ "Osteocyte",
      name %in% c("SERPINE2", "CXCL12") ~ "BMSC",
      name %in% c("ACP5", "CST3") ~ "Osteoclast",
      name %in% c("COL1A1", "COMP", "CLEC11A", "STC2", "COL1A2", "OMD", "ECM2", "FMOD", "MYOC", "FBLN7") ~ "Osteoblast & Matrix",
      TRUE ~ "Chondrocyte"
    ),
    # 生成带有数值标签的完美文本格式 (例如：SOST\n-3.12)
    label_text = sprintf("%s\n%.2f", name, Log2FC)
  )

# ------------------------------------------------------------------------------
# 5. 构建网络与精准坐标映射 (保证排版永不重叠)
# ------------------------------------------------------------------------------
g <- tbl_graph(nodes = nodes_df, edges = edges_df, directed = FALSE) %>%
  mutate(Degree = centrality_degree())

coords_df <- data.frame(
  name = c("SOST", "DKK1", "STC2", "FMOD", "COL1A1", "COMP", "CLEC11A", "COL1A2", "MYOC", "FBLN7", "ECM2", "OMD",
           "CXCL12", "SERPINE2", "ACP5", "CST3", "PRG4", "F13A1", "TNC", "MIA", "FRZB", "PRG2", "SMOC2", "MGP", "COL11A1"),
  custom_x = c(0.25, -0.65, -0.70, 0.25, -1.40, -2.10, -2.35, -1.30, 0.95, 1.45, 1.90, 2.30, 
               2.45, 2.75, 0.85, 0.15, -1.75, -2.35, -2.95, -2.15, -2.85, -0.55, 0.35, 1.15, 1.30),
  custom_y = c(0.05, 0.05, 1.35, 1.35, 2.00, 1.45, 2.30, 2.90, 1.90, 2.65, 1.40, 0.65, 
               2.10, 1.25, -0.55, -2.55, -0.75, -1.20, -0.65, -2.05, -2.55, -1.60, -2.05, -2.50, -1.45)
)

layout_data <- create_layout(g, layout = "fr")
layout_data$x <- coords_df$custom_x[match(layout_data$name, coords_df$name)]
layout_data$y <- coords_df$custom_y[match(layout_data$name, coords_df$name)]

distinct_cell_colors <- c("Osteocyte" = "#D62728", "Osteoblast & Matrix" = "#2CA02C", 
                          "Chondrocyte" = "#1F77B4", "BMSC" = "#9467BD", "Osteoclast" = "#FF7F0E")

# ------------------------------------------------------------------------------
# 6. 图表渲染 (渐变色映射 + 真实数值双行标签)
# ------------------------------------------------------------------------------
p_final <- ggraph(layout_data) +
  
  # 图层 1：单细胞来源领地 (背景)
  geom_mark_hull(aes(x = x, y = y, fill = CellType, color = CellType, group = CellType),
                 concavity = 1.8, expand = unit(7, "mm"), radius = unit(6, "mm"),
                 alpha = 0.12, linewidth = 0.5, show.legend = TRUE) +
  scale_fill_manual(values = distinct_cell_colors, name = "Cellular Origin\n(Niche Territory)") +
  scale_color_manual(values = distinct_cell_colors, guide = "none") +
  
  new_scale_fill() +  # 释放 Fill 给节点的渐变色使用
  
  # 图层 2：跨细胞相关性网络连线
  geom_edge_link(aes(width = weight, color = Direction), alpha = 0.85) +
  scale_edge_width(range = c(0.8, 2.6), name = "|Correlation|") +
  scale_edge_color_manual(values = c("Positive" = "gray50", "Negative" = "#00A087")) +
  
  # 图层 3：节点中心色映射真实表达量变化 (Log2FC)
  geom_node_point(aes(size = Degree, fill = Log2FC), shape = 21, color = "black", stroke = 1.2) +
  scale_fill_gradient2(low = "#313695", mid = "white", high = "#E64B35", midpoint = 0, 
                       name = "Exercise Response\n(Log2 Fold Change)") +
  scale_size_continuous(range = c(6, 13), name = "Network Degree") +
  
  # 图层 4：双行基因标签 (基因名 + 真实 Log2FC 值)
  geom_node_text(aes(label = label_text), repel = TRUE, size = 3.6, fontface = "bold.italic", 
                 color = "black", lineheight = 0.8, bg.color = "white", bg.r = 0.15) +
  
  # 注释标签
  # annotate("text", x = -0.20, y = 0.55, label = "★ Primary Responsive Hub", 
  #          color = "#B30000", fontface = "bold", size = 4.2, hjust = 0.5) +
  
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5))) +
  labs(title = "Exercise-Induced Multi-Cellular Secretome Coordination")

# 输出高清成果图
print(p_final)