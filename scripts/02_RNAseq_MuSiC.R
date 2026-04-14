######  MuSiC 反卷积算法  #####
# 1. 加载必要的包
library(Seurat)
library(MuSiC)
library(Biobase)
library(SingleCellExperiment)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)

# 步骤 1：处理 Bulk RNA-seq 数据 (TPM_value)
print("正在构建 Bulk 数据的 ExpressionSet 对象...")
# 确保 TPM_value 是矩阵格式，并转换为 MuSiC 需要的 ExpressionSet
bulk_matrix <- as.matrix(TPM_value)
bulk_eset <- Biobase::ExpressionSet(assayData = bulk_matrix)

# 步骤 2：处理单细胞数据 (pbmc)
print("正在准备单细胞参考数据集...")
# ⚠️ 注意：MuSiC 算法强制要求单细胞数据中必须包含“样本来源(Sample)”的信息
if (!"SampleID" %in% colnames(pbmc@meta.data)) {
  pbmc$SampleID <- pbmc$orig.ident 
}

# 提取单细胞的原始 count 矩阵进行运算 (反卷积必须用 counts 而不是 scale_data)
sce_ref <- as.SingleCellExperiment(pbmc, assay = "RNA")

# 步骤 3：正式运行 MuSiC 算法
print("🚀 正在运行 MuSiC 反卷积算法，请稍候 (根据细胞数量可能需要几分钟)...")
music_results <- music_prop(
  bulk.mtx = Biobase::exprs(bulk_eset),  # Bulk 表达矩阵
  sc.sce = sce_ref,                      # 单细胞参考集
  clusters = 'cell_type',                # 你指定的细胞亚群列名
  samples = 'SampleID',                  # 样本来源列名
  verbose = TRUE
)

# 步骤 4：提取结果并保存
estimated_proportions <- data.frame(music_results$Est.prop.weighted)  

# 1.  decon_res 变量
decon_res <- music_results 

# 2. 从反卷积结果生成 df_plot
df_plot <- estimated_proportions
df_plot$Sample <- rownames(df_plot)

# 3. 提取分组信息 (Group)。
df_plot$Group <- gsub(".TPM*", "", df_plot$Sample) 



#######  IMC_Only_Exact_Pvalue ######

# 1. 专门提取 IMC 数据
df_imc <- df_plot %>% select(Sample, Group, Proportion = IMC)
df_imc$Group <- factor(df_imc$Group, levels = c("Y", "O", "S"))

# 2. 绘制单独的 IMC 比例图
p_imc <- ggplot(df_imc, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = c("Y" = "#2ca02c", "O" = "#d62728", "S" = "#1f77b4")) +
  theme_classic(base_size = 15) +
  labs(title = "Inferred Proportion of Kcnk2+ IMCs",
       x = "", y = "Cell Proportion") +
  scale_x_discrete(labels=c("Y" = "Young", "O" = "Aged\nSed.", "S" = "Aged\nExe.")) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  stat_compare_means(comparisons = list(c("Y", "O"), c("O", "S")),
                     method = "t.test",
                     label = "p.format")

ggsave("IMC_Only_Exact_Pvalue.pdf", plot = p_imc, width = 4, height = 5)
print(p_imc)


####### Ref的选择 #####

DefaultAssay(pbmc) <- "RNA" # 切换到原始表达 Assay

# 1. 绘制多基因小提琴图 (Violin Plot)
p_vln <- VlnPlot(pbmc,
                 features = c("Col1a1", "Sp7"),
                 group.by = "cell_type", 
                 pt.size = 0,
                 ncol = 1) + 
  theme(axis.title.x = element_blank()) + 
  plot_annotation(title = "Single-cell Expression Comparison")

# 2. 绘制 UMAP 表达分布图 (Feature Plot)
p_feat <- FeaturePlot(pbmc,
                      features = c("Col1a1", "Sp7"),
                      ncol = 2,
                      order = TRUE) + 
  plot_annotation(title = "UMAP Feature Distribution")

# 3. 整合输出并保存
final_validation_plot <- p_vln | p_feat
ggsave("Col1a1_Reference_Validation.pdf", plot = final_validation_plot, width = 12, height = 8)
print(final_validation_plot)


######### Lineage-Normalized Pathogenic Burden #########

print("🌟 步骤 1: 选取靶标基因与内参基因")
target_genes <- c("Dkk1", "Sost", "Kcnk2", "Enpp1")
reference_gene <- "Col1a1"  

valid_genes <- intersect(target_genes, rownames(TPM_value))
if (!(reference_gene %in% rownames(TPM_value))) {
  stop(paste("错误：内参基因", reference_gene, "不在 TPM 矩阵中，请更换为 Col1a1 或 Runx2"))
}

print("🌟 步骤 2: 计算谱系内相对分泌指数 (Relative Burden Index)")
ref_expr <- as.numeric(TPM_value[reference_gene, ]) + 0.01

norm_expr <- data.frame(Sample = colnames(TPM_value))
for (gene in valid_genes) {
  norm_expr[[gene]] <- as.numeric(TPM_value[gene, ]) / ref_expr
}

print("🌟 步骤 3: 与 IMC 比例合并并作图")
corr_df_norm <- merge(df_plot[, c("Sample", "Group", "IMC")], norm_expr, by = "Sample")

corr_long_norm <- corr_df_norm %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Norm_Expression")

corr_long_norm$Gene <- factor(corr_long_norm$Gene, levels = valid_genes)
corr_long_norm$Group <- factor(corr_long_norm$Group, levels = c("Y", "O", "S"))

print("🌟 步骤 4: 绘制归一化后的完美回归图")
p_corr_norm <- ggplot(corr_long_norm, aes(x = IMC, y = Norm_Expression)) +
  geom_smooth(method = "lm", color = "#b22222", fill = "grey80", alpha = 0.4) + 
  geom_point(aes(fill = Group), shape = 21, size = 4, color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("Y" = "#2ca02c", "O" = "#d62728", "S" = "#1f77b4")) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  theme_classic(base_size = 14) +
  labs(title = "Lineage-Normalized Pathogenic Burden",
       x = "Inferred Proportion of Kcnk2+ IMCs",
       y = paste0("Relative Expression (Normalized to ", reference_gene, ")")) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold.italic", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Normalized_IMC_Correlation.pdf", plot = p_corr_norm, width = 10, height = 4)
print(p_corr_norm)


######  Z-score 靶基因集打分法量化 p53 活性联合分析  #####

print("🌟 步骤 1: 提取参考基因 (Col1a1) 并进行相对定量校准")
col1a1_expr <- as.numeric(bulk_matrix["Col1a1", ]) + 1e-6
p53_targets <- c("Cdkn1a", "Mdm2", "Bax", "Phlda3", "Ccng1", "Gadd45a", "Zmat3", "Bbc3")
valid_targets <- intersect(p53_targets, rownames(bulk_matrix))
target_expr <- bulk_matrix[valid_targets, ]

target_expr_norm <- sweep(target_expr, 2, col1a1_expr, FUN = "/")

z_expr_norm <- t(scale(t(target_expr_norm)))
p53_activity_norm <- colMeans(z_expr_norm, na.rm = TRUE)

print("🌟 步骤 3: 重新构建数据框 (使用归一化后的数据)")
imc_proportion <- decon_res$Est.prop.weighted[, "IMC"] 

df_analysis_norm <- data.frame(
  Sample = colnames(bulk_matrix),
  Group = metadata$Group, 
  IMC_Proportion = as.numeric(imc_proportion),
  p53_Activity_Norm = as.numeric(p53_activity_norm)
)

df_analysis_norm$Group <- factor(df_analysis_norm$Group, levels = c("Y", "O", "S"))
my_comparisons <- list(c("Y", "O"), c("O", "S"))

plot_box_norm <- function(y_var, y_label, title) {
  ggboxplot(df_analysis_norm, x = "Group", y = y_var,
            fill = "Group",
            palette = c("#4DBBD5", "#E64B35", "#00A087"), 
            width = 0.5, size = 0.4, outlier.shape = NA, add = "jitter", 
            add.params = list(size = 1.2, color = "#4d4d4d", alpha = 0.8, jitter = 0.15)) +
    stat_compare_means(comparisons = my_comparisons, 
                       method = "t.test", label = "p.value", 
                       step.increase = 0.12, bracket.size = 0.3, 
                       tip.length = 0.02, size = 3.5) +       
    labs(title = title, x = "", y = y_label) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(linewidth = 0.4), 
      axis.ticks = element_line(linewidth = 0.4) 
    )
}

b1_norm <- plot_box_norm("p53_Activity_Norm", "Col1a1-Normalized Score", "p53 Pathway Activity")
b2_norm <- plot_box_norm("IMC_Proportion", "Estimated Proportion", "IMC (Kcnk2+) Accumulation")
box_plot_final <- (b1_norm | b2_norm) 
print(box_plot_final)

ggsave("Bulk_Validation_Boxplots_Normalized.pdf", box_plot_final, width = 10, height = 8)
 
######  Tp53 逆转象限图  ###### 

# 1. 构建 Bulk 的 Metadata（给 p53 Z-score 那块用）
metadata <- data.frame(
  Sample = df_plot$Sample,
  Group = df_plot$Group
)
rownames(metadata) <- metadata$Sample

# 2. 计算各个组别的平均 TF 活性
mean_Y <- rowMeans(tf_activity[, metadata$Sample[metadata$Group == "Y"], drop=FALSE])
mean_O <- rowMeans(tf_activity[, metadata$Sample[metadata$Group == "O"], drop=FALSE])
mean_S <- rowMeans(tf_activity[, metadata$Sample[metadata$Group == "S"], drop=FALSE])

# 3. 计算差异效应量
act_diff_z <- data.frame(
  TF = rownames(tf_activity),
  Delta_Aging = mean_O - mean_Y,      # 老化效应：老 - 年轻
  Delta_Exercise = mean_S - mean_O    # 运动干预效应：运动 - 老
)

print("🌟 步骤 1: 抓取 Top 10 拯救因子与 Tp53")
act_diff_z_clean <- act_diff_z %>%
  filter(!grepl("^A0A", TF))

top_rescued_tfs_clean <- act_diff_z_clean %>%
  filter(Quadrant == "Rescued by Exercise") %>%
  arrange(desc(Distance)) %>%
  head(10) %>%
  pull(TF)

tfs_to_label_clean <- unique(c("Tp53", top_rescued_tfs_clean))

act_diff_z_clean <- act_diff_z_clean %>%
  mutate(Label = ifelse(TF %in% tfs_to_label_clean, TF, ""))

print("🌟 步骤 2: 绘制宽幅舒适版全基因组逆转图")
# 动态调整极限值，防止模拟数据范围不够 30 导致报错，使用真实数据时您可以解除限制
max_x <- max(30, max(act_diff_z_clean$Delta_Aging) * 1.5)

final_plot_spaced <- ggplot(act_diff_z_clean, aes(x = Delta_Aging, y = Delta_Exercise)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(aes(color = Quadrant, size = Quadrant), alpha = 0.85) +
  scale_color_manual(values = c("Rescued by Exercise" = "#E64B35", 
                                "Reversed Up" = "#4DBBD5", 
                                "Not Rescued/Unchanged" = "grey80")) +
  scale_size_manual(values = c("Rescued by Exercise" = 3.5, 
                               "Reversed Up" = 3.5, 
                               "Not Rescued/Unchanged" = 1.5)) +
  scale_x_continuous(limits = c(min(act_diff_z_clean$Delta_Aging) - 2, max_x)) +
  geom_text_repel(aes(label = Label),
                  size = 5,                  
                  fontface = "bold.italic",   
                  color = "black",
                  box.padding = 1.2,          
                  point.padding = 0.5,        
                  nudge_x = 2,                
                  segment.color = "grey30",
                  segment.linewidth = 0.5,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  seed = 42) +                
  labs(title = "Global Transcriptional Network Reversal",
       subtitle = "Bottom-right: Transcription factors activated by aging and suppressed by exercise",
       x = "Aging Effect (\u0394 Activity, Old - Young)",
       y = "Exercise Effect (\u0394 Activity, Exe - Old)") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10) 
  )

print(final_plot_spaced)
ggsave("Tp53_Reversal_Quadrant_Plot.pdf", final_plot_spaced, width = 10, height = 8)