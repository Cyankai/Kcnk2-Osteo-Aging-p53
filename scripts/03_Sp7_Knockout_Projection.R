#####  锚点投射 ####
library(Seurat)
library(ggplot2)

# 1. 切换到 RNA Assay 
DefaultAssay(pbmc_yours) <- "RNA"
DefaultAssay(pbmc_public) <- "RNA"

# 2. 找到两个数据集共有的基因名
common_features <- intersect(rownames(pbmc_yours), rownames(pbmc_public))

# 3. 在共有基因中，寻找参考集 (pbmc_yours) 的高变基因
pbmc_yours <- FindVariableFeatures(pbmc_yours, nfeatures = 2000)
transfer_features <- intersect(VariableFeatures(pbmc_yours), common_features)

# 4. 重新寻找锚点 (显式提供 features 参数)
transfer_anchors <- FindTransferAnchors(
  reference = pbmc_yours, 
  query = pbmc_public, 
  dims = 1:30,
  features = transfer_features,       
  reference.reduction = "pca",
  normalization.method = "LogNormalize"
)

# 5. 执行标签映射
predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = pbmc_yours$cell_type,    # 确保这是你定义 IMC 的那一列
  dims = 1:30
)

# 6. 将结果写回并检查
pbmc_public <- AddMetaData(pbmc_public, metadata = predictions)

# 检查原作者的分类里，到底有多少细胞被识别成了你的 IMC
table(Original = pbmc_public$cell_type, Predicted = pbmc_public$predicted.id)

# 7. 可视化：看看这些被“抢出来”的 IMC 细胞分布在哪里
DimPlot(pbmc_public, group.by = "predicted.id", label = TRUE) + 
  ggtitle("Predicted Identities in Sp7 KO Data")


##### Lineage Deviation: Massive Expansion of IMCs in Sp7 KO #####
library(dplyr)

meta_public <- pbmc_public@meta.data

prop_data <- meta_public %>%
  group_by(orig.ident, predicted.id) %>% 
  summarise(Count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(Proportion = Count / sum(Count))

ggplot(prop_data, aes(x = orig.ident, y = Proportion, fill = predicted.id)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.6) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Lineage Deviation: Massive Expansion of IMCs in Sp7 KO",
       x = "Group", y = "Cell Proportion", fill = "Predicted Cell Type") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


######  Ectopic Antagonist Surge Upon Sp7 Ablation ####
library(tidyr)

genes_to_plot <- c("Dkk1", "Sost") 
expr_sp7 <- FetchData(pbmc_public, vars = c(genes_to_plot, "orig.ident", "predicted.id"))
colnames(expr_sp7)[3] <- "cell_type" 
target_cells <- c("LMP", "IMC", "Osteoblast", "Osteocyte")
expr_sp7 <- expr_sp7 %>% filter(cell_type %in% target_cells)

summary_sp7 <- expr_sp7 %>%
  group_by(orig.ident, cell_type) %>%
  summarise(across(all_of(genes_to_plot), sum), .groups = "drop") %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression")

delta_labels_sp7 <- summary_sp7 %>%
  pivot_wider(names_from = orig.ident, values_from = Expression) %>%
  mutate(Delta = KO - WT) %>%
  filter(Delta > 0) %>% 
  mutate(Label = sprintf("Δ+%.1f", Delta),
         Y_pos = KO,
         orig.ident = "KO") 

summary_sp7$cell_type <- factor(summary_sp7$cell_type, levels = target_cells)

ggplot(summary_sp7, aes(x = cell_type, y = Expression, fill = orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7, color = "black", linewidth = 0.5) +
  facet_wrap(~ Gene, scales = "free_y") + 
  geom_text(data = delta_labels_sp7,
            aes(x = cell_type, y = Y_pos, label = Label),
            position = position_dodge(width = 0.8),
            vjust = -0.5, fontface = "bold", color = "#C0392B", size = 3.5) +
  scale_fill_manual(values = c("WT" = "#3498DB", "KO" = "#E74C3C")) +
  theme_classic() +
  labs(title = "Ectopic Antagonist Surge Upon Sp7 Ablation",
       subtitle = "Floating annotations (Δ) represent the increment induced by maturation arrest",
       x = "Osteogenic Lineage", y = "Total Transcript Abundance") +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top"
  )


#########  Cell Cycle Arrest Induced by Sp7 Ablation  #####
library(stringr)

s.genes.mouse <- str_to_title(cc.genes$s.genes)
g2m.genes.mouse <- str_to_title(cc.genes$g2m.genes)

pbmc_public <- CellCycleScoring(
  object = pbmc_public,
  s.features = s.genes.mouse,
  g2m.features = g2m.genes.mouse,
  set.ident = FALSE
)

meta_sp7 <- pbmc_public@meta.data 
meta_sp7$cell_type <- meta_sp7$predicted.id 

meta_sp7 <- meta_sp7 %>% filter(cell_type %in% target_cells)
meta_sp7$cell_type <- factor(meta_sp7$cell_type, levels = target_cells)
meta_sp7$Phase <- factor(meta_sp7$Phase, levels = c("G1", "S", "G2M"))

cc_sp7_summary <- meta_sp7 %>%
  group_by(orig.ident, cell_type, Phase) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(orig.ident, cell_type) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

ggplot(cc_sp7_summary, aes(x = cell_type, y = Proportion, fill = Phase)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.4) +
  facet_wrap(~ orig.ident) + 
  theme_classic() +
  scale_fill_manual(values = c("G1" = "#E74C3C", "S" = "#F1C40F", "G2M" = "#3498DB")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cell Cycle Arrest Induced by Sp7 Ablation",
       subtitle = "Comparing cell cycle phase distribution across osteolineage",
       x = "Osteogenic Trajectory", y = "Proportion of Cells", fill = "Cell Cycle Phase") +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.8),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 12, hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.spacing = unit(1.5, "lines") 
  )


###### Sp7 GO #####  
library(clusterProfiler)
library(org.Mm.eg.db)

Idents(pbmc_public) <- "predicted.id"
imc_public <- subset(pbmc_public, idents = "IMC")
Idents(imc_public) <- "orig.ident"

imc_deg <- FindMarkers(imc_public, 
                       ident.1 = "KO", 
                       ident.2 = "WT", 
                       logfc.threshold = 0, 
                       min.pct = 0.05)      

imc_deg$gene <- rownames(imc_deg)
ranked_df <- imc_deg %>% arrange(desc(avg_log2FC))
gene_list <- ranked_df$avg_log2FC
names(gene_list) <- ranked_df$gene

gsea_go <- gseGO(geneList     = gene_list,
                 OrgDb        = org.Mm.eg.db,
                 keyType      = "SYMBOL",
                 ont          = "BP",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.25,
                 eps          = 0,            
                 nPermSimple  = 10000,        
                 pAdjustMethod= "BH",
                 verbose      = FALSE)

gsea_df <- as.data.frame(gsea_go)
 
mito_terms <- gsea_df %>% 
  filter(NES < 0) %>% 
  filter(grepl("mitochondr|mito", Description, ignore.case = TRUE)) %>% 
  head(10)

mito_terms$Description <- paste0("GOBP_", toupper(gsub(" ", "_", mito_terms$Description)))
mito_terms$minus_NES <- abs(mito_terms$NES) 
mito_terms$minus_log10_FDR <- -log10(mito_terms$p.adjust)

mito_terms$Description <- factor(mito_terms$Description, 
                                 levels = mito_terms$Description[order(mito_terms$minus_NES)])

ggplot(mito_terms, aes(x = minus_NES, y = Description)) +
  geom_point(aes(size = setSize, color = minus_log10_FDR)) +
  scale_color_gradient(low = "blue", high = "magenta", name = bquote(-log[10]~"(FDR)")) +
  scale_size_continuous(name = "Counts", range = c(3, 8)) +
  theme_bw() +
  labs(x = "-NES", y = "") +  
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", linewidth = 1)
  )