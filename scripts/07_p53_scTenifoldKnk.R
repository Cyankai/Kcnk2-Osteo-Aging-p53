# ==============================================================================
# Script Name: 07_p53_scTenifoldKnk.R
# Description: Conducts an in silico gene perturbation (virtual knockout) of Trp53 
#              (p53) using the scTenifoldKnk machine learning framework. Evaluates 
#              the structural gene regulatory network (GRN) shifts and the rescue 
#              of SASP and Wnt antagonists within the aged Kcnk2+ IMC subpopulation.
# Input Data:  Sub-setted Seurat object of 16M Kcnk2+ cells.
# Output:      'p53_ko_results.rds' (GRN adjacency matrices, manifold alignment).
# ==============================================================================
# =================================================================
# 第一步：加载包 
# =================================================================
library(Seurat)
library(scTenifoldKnk) 
library(dplyr)
library(Matrix)

# =================================================================
# 第二步：安全的 Subsetting 操作与原始 Count 提取
# =================================================================
# 提取 16M 的 Kcnk2+ 衰老亚群
kcnk2_aged_cells <- subset(pbmc, subset = orig.ident == "16M" & cell_type == "Kcnk2+Cell")

# 强制从原始 RNA Assay 提取未标准化的 Counts (Seurat v5 写法)
count_matrix <- LayerData(kcnk2_aged_cells, assay = "RNA", layer = "counts")

# (如果是 Seurat v4，请换用下面这句：)
# count_matrix <- GetAssayData(kcnk2_aged_cells, assay = "RNA", slot = "counts")

# =================================================================
# 第三步：加入“白名单”的基因过滤 (极为关键)
# =================================================================
# 1. 明确我们需要重点关注的靶点 (敲除靶点 + 观察靶点)
core_targets <- c("Trp53", "Dkk1", "Sost", "Limch1", "Kcnk2")

# 2. 基础过滤逻辑：在至少 5% 的目标细胞中表达量大于 0
expressed_genes <- rowSums(count_matrix > 0) >= (ncol(count_matrix) * 0.05)

# 3. 🚨 白名单强制保留：不论表达多低，只要这些靶点存在于矩阵中，就强制设为 TRUE
present_targets <- core_targets[core_targets %in% rownames(count_matrix)]
expressed_genes[present_targets] <- TRUE

# 4. 执行过滤
count_matrix_filtered <- count_matrix[expressed_genes, ]

print(paste("过滤后剩余的基因数：", nrow(count_matrix_filtered)))
print(paste("参与构建网络的细胞数：", ncol(count_matrix_filtered)))
print(paste("成功强制保留的核心靶点数：", length(present_targets)))

# =================================================================
# 第四步：执行 p53 (Trp53) 虚拟敲除 (scTenifoldKnk)
# =================================================================
set.seed(2026)
p53_ko_results <- scTenifoldKnk(countMatrix = count_matrix_filtered, 
                                gKO = "Trp53",     # 小鼠的 p53 基因名
                                nc_nNet = 16)      # 线程数请根据你的服务器配置调整 (如 10 到 80)

# =================================================================
# 第五步：追踪目标标志物的因果逆转
# =================================================================
# 智能提取差异表达表格
if(is.data.frame(p53_ko_results)) {
  diff_reg <- p53_ko_results
} else {
  df_index <- which(sapply(p53_ko_results, is.data.frame))[1]
  diff_reg <- p53_ko_results[[df_index]]
}

# 设定你要验证的衰老与结构靶点
validation_genes <- c("Dkk1", "Sost", "Limch1", "Kcnk2")

# 提取并按表达变化 (FC) 排序
validation_targets <- diff_reg %>% 
  filter(gene %in% validation_genes) %>%
  select(gene, distance, FC, p.value, p.adj) %>%
  arrange(FC)

cat("\n=== 虚拟敲除 Trp53 (p53) 后核心靶基因的变化 ===\n")
print(validation_targets)

save(p53_ko_results, file='p53_ko_results.Rdata' )

# 1. 检查结果对象是否存在且非空
exists("p53_ko_results")
str(p53_ko_results)

# 2. 检查输出数据框的前几行与维度
head(p53_ko_results)
dim(p53_ko_results)

# 3. 检查关键下游靶点 (如 Kcnk2, Cdkn1a 等) 的扰动排序与显著性
p53_ko_results[p53_ko_results$gene %in% c('Trp53',"Kcnk2", "Cdkn1a", "Enpp1", "Dkk1", "Sost"), ]

# 提取核心差异调控数据框
diff_reg <- p53_ko_results$diffRegulation

# 1. 查看受 Trp53 扰动最显著的 Top 20 基因
head(diff_reg, 20)

# 2. 检查关键通路基因 (Kcnk2, Cdkn1a, Enpp1, Dkk1, Sost, Limch1)
key_targets <- c("Kcnk2", "Cdkn1a", "Enpp1", "Dkk1", "Sost", "Limch1", "Il6", "Serpine1")
diff_reg[diff_reg$gene %in% key_targets, ]

# 筛选显著受扰动基因 (根据标准可设 p.adj < 0.05 或 Z > 1.96)
sig_perturbed_genes <- diff_reg[diff_reg$p.adj < 0.05 & diff_reg$Z > 1.96, ]
save(p53_ko_results, file='p53_ko_results.Rdata' )