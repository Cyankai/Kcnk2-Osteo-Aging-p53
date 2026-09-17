# ==============================================================================
# Script Name: 03_Kcnk2_aging_scRNAseq.R
# Description: The core single-cell RNA-seq analysis pipeline for the bone niche. 
#              Functions include: data integration (3M vs 16M), cell clustering, 
#              identification of the Intermediate Cluster (IMC), differentiation 
#              trajectory inference (Monocle3), apoptosis scoring, and downstream 
#              pathway/GSEA analysis of the Kcnk2+ subpopulation.
# Input Data:  Raw count matrices (1M, 3M, 16M WT samples).
# Output:      Seurat objects, UMAP embeddings, Monocle3 trajectories, and DEGs.
# ==============================================================================
#### 16M clusters annotation ####
library(Seurat)
library(ggplot2)
library(dplyr)
options(future.globals.maxSize = 16 * 1024^3)
# 1. 创建样本列表（假设已加载3个样本）
pbmc1=as.data.frame(read.table('/Users/kenny/Documents/seq/sc_seq/data/Qin lin/1M/GSM4318799_1M_matrix.txt', row.names = 1, header = T, sep = " "))
pbmc2=as.data.frame(read.table('/Users/kenny/Documents/seq/sc_seq/data/Qin lin/3M/GSM4318801_3M_matrix.txt', row.names = 1, header = T, sep = " "))
pbmc3=as.data.frame(read.table('/Users/kenny/Documents/seq/sc_seq/data/Qin lin/16M/GSM4318802_16M_matrix.txt', row.names = 1, header = T, sep = " "))
pbmc4=Read10X(data.d16M = '/Users/kenny/Desktop/Plpp1/P 骨/scRNA-seq/整体results/matrix/WT/')
pbmc1 =  CreateSeuratObject(counts = pbmc1, project = "1M", min.cells = 2, min.features = 100)
pbmc2 =  CreateSeuratObject(counts = pbmc2,  project = "3M",min.cells = 2, min.features = 100)
pbmc3 =  CreateSeuratObject(counts = pbmc3,  project = "16M",min.cells = 2, min.features = 100)
pbmc4 =  CreateSeuratObject(counts = pbmc4,  project = "WT",min.cells = 2, min.features = 100)

(obj.list <- list(  pbmc1 ,  pbmc2 ,pbmc3 ,pbmc4 ))
rm( pbmc1 , pbmc2 ,pbmc3,pbmc4 )

# 先线粒体质控
(obj.list <- lapply(X = obj.list, FUN = function(x) {
  x[['percent.mt']] <- PercentageFeatureSet(x, pattern = "^mt-")  
  x <- subset(x , subset = nFeature_RNA > 50 & nFeature_RNA < 6000 )
}))

# 再进行 SCTransform
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2")  # 每个样本独立SCT标准化
})  

# 3. 选择整合特征
features <- SelectIntegrationFeatures(
  object.list = obj.list, 
  nfeatures = 2000  # 推荐值
)

# 4. 准备SCTransform整合
obj.list <- PrepSCTIntegration(
  object.list = obj.list,
  anchor.features = features,
  verbose = FALSE
)

# 5. 寻找锚点
anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  normalization.method = "SCT",  # 必须指定
  anchor.features = features,    # 使用上一步选择的特征
  verbose = FALSE
)
rm(obj.list)
# 6. 整合数据
integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  verbose = FALSE,
)
rm(anchors)
# 7. 后续分析
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30) 
integrated <- FindNeighbors(integrated, dims = 1:30)
res.used <- c(0.01,0.03,0.05,0.07,.09)
integrated <- FindClusters(integrated, resolution = res.used)

# Make plot 
library(clustree)
clustree(integrated@meta.data, prefix =
           "integrated_snn_res."
) 

# 确定resolution
final_resolution = 0.09
integrated <- FindClusters(integrated, resolution = final_resolution)

# 确定成骨亚群
MSCs= c('Prrx1','Cxcl12','Lepr','Ebf3','Pdgfrb','Pdgfra') 
Chondrocyte = c('Sox9','Col2a1','Col10a1','Pth1r','Acan','Ihh') 
obot=c("Col1a1","Bglap","Sp7","Alpl","Dmp1","Ptprz1")
VlnPlot(integrated, features = MSCs)
VlnPlot(integrated, features = obot)  #  2\6\14 号亚群
VlnPlot(integrated, features = Chondrocyte)  #  8 号亚群

# 6. 可视化批次效应
DimPlot(integrated, label = T, pt.size = 0.2,label.size = 4,
        #  group.by = 'cell_type',
        split.by= 'orig.ident' )

# Ob亚群
pbmc <- subset(integrated, subset = seurat_clusters %in% c(2,6,14) )
pbmc <- subset(pbmc, subset =  orig.ident %in% c("16M", "3M"))
pbmc=SCTransform(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10,reduction='pca')
res.used <- c(0.01,0.03,0.05,0.07,.09)
pbmc <- FindClusters(pbmc, resolution = res.used)
library(clustree)
clustree(pbmc@meta.data, prefix =
           # "RNA_snn_res.",
           "SCT_snn_res."
         #'integrated_snn_res.'
) 
final_resolution = 0.07
pbmc <- FindClusters(pbmc, resolution = final_resolution)
table(pbmc@meta.data$seurat_clusters)
pbmc <- RunUMAP(pbmc , dims = 1:30,reduction='pca')
DimPlot(pbmc,  label = T, pt.size = 0.7,label.size = 4,
        split.by= 'orig.ident')
as.data.frame.matrix(table(Idents(pbmc), pbmc$orig.ident))
pbmc = PrepSCTFindMarkers(pbmc)
pbmc.allmarkers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
(pbmc.allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10_pbmc)

# 组间比较
cell_type = c("IMC",'LMP',"Osteoblast","Osteocyte")
pbmc[['cell_type']] = unname(cell_type[pbmc@meta.data$seurat_clusters])
pbmc$cell_type <- factor(pbmc$cell_type, levels = c("LMP", "Osteoblast", "IMC", "Osteocyte"))
DotPlot(pbmc, features = unique(c( 'Aspn','Edil3','Tnn','Bglap2',"Pdzd2", 'Bglap','Limch1',"Enpp1",'Kcnk2', 'Ptprz1', "Dmp1" , "Sost")),
        group.by = "cell_type",
        cols = c("lightgrey", "red"),  # 颜色梯度
        scale.min = 0,  # 最小表达阈值
        scale.max = 50,   # 最大表达阈值
        dot.scale = 6            # 点大小
) +
  RotatedAxis() +
  labs(x = "", y = "Sample") +
  theme(
    legend.text = element_text(size = 10, color = "black", family = "Arial"),  # 图例文字
    legend.title = element_text(size = 12, face = "bold")                     # 图例标题
  )

DimPlot(pbmc,  label = T, pt.size = 0.7,label.size = 4,group.by = "cell_type",
        split.by= 'orig.ident')

save(pbmc,file = '/Users/kenny/Desktop/pbmc_yours.RData')

#### Total Expression of Wnt Antagonists (3M vs 16M) #######
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
# 1. 提取所有细胞的表达数据和 Metadata
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc <- JoinLayers(pbmc)
expr_data_all <- GetAssayData(pbmc, assay = "RNA", layer = "data")
meta_data_all <- pbmc@meta.data
# 1. 明确你要提取的基因和 Metadata 列名
# 假设你的 age 存放在 pbmc$age, 细胞类型存放在 pbmc$cell_type
genes_to_get <- c("Dkk1", "Sost")
metadata_to_get <- c("orig.ident", "cell_type")

# 2. 使用 FetchData 提取
expr_df_all <- FetchData(pbmc, vars = c(genes_to_get, metadata_to_get))

# 2. 计算每个年龄组、每个亚群的总表达量
summary_long <- expr_df_all %>%
  group_by(orig.ident, cell_type) %>%
  summarise(Dkk1 = sum(Dkk1), Sost = sum(Sost), .groups = "drop") %>%
  pivot_longer(cols = c(Dkk1, Sost), names_to = "Gene", values_to = "Expression")
# 排序
summary_long$orig.ident <- factor(summary_long$orig.ident, levels = c("3M", "16M"))
summary_long$cell_type <- factor(summary_long$cell_type, 
                                 levels = c("LMP", "Osteoblast", "Kcnk2+Cell", "Osteocyte"))
# 3. 核心计算：计算 16M 的增量 (Delta) 用于图表上的文字悬浮标注
delta_labels <- summary_long %>%
  pivot_wider(names_from = orig.ident, values_from = Expression) %>%
  mutate(
    Delta = `16M` - `3M`,
    # 格式化增量标签：正数加"+"号，负数保留"-"号
    Label = ifelse(Delta > 0, sprintf("Δ = +%.1f", Delta), sprintf("Δ = %.1f", Delta)),
    # 将标签的位置稍微放在两根柱子中最高的那根上面一点点 (防遮挡)
    Y_pos = pmax(`16M`, `3M`) + (max(pmax(`16M`, `3M`)) * 0.05),
    # 关键一步：把标签的 age 设为 "16M"，这样它在画图时就会自动对齐到 16M 的柱子正上方！
    orig.ident = "16M"
  )
# 4. 绘制终极版高分可视化图
ggplot(summary_long, aes(x = cell_type, y = Expression, fill = orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
  facet_wrap(~ Gene, scales = "free_y") + # 按照 Dkk1 和 Sost 分成两张独立的子图
  # 添加 16M 增量悬浮标签
  geom_text(data = delta_labels,
            aes(x = cell_type, y = Y_pos, label = Label),
            position = position_dodge(width = 0.8),
            vjust = 0, fontface = "bold", color = "#C0392B", size = 4) +
  # 使用顶级期刊极其经典的对比色系 (清澈蓝 vs 警示红)
  scale_fill_manual(values = c("3M" = "#3498DB", "16M" = "#E74C3C")) +
  theme_classic() +
  labs(y = "Total Transcript Burden (Sum)", fill = "Age Group") +
  # 极致的美化排版
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, face = "italic", hjust = 0.5, color = "gray30"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  ) 

#### Monocle3 ####
library(Seurat) 
library(monocle3)
library(v16Midis)
library(patchwork)
pbmc = PrepSCTFindMarkers(pbmc)
seurat_obj=NormalizeData(pbmc)

# 提取表达矩阵（建议使用标准化数据）
expr_matrix <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")  # log-normalized

# 提取细胞元数据（必须包含UMAP/tSNE坐标和聚类信息）
cell_metadata <- seurat_obj@meta.data

# 提取基因注释（需为data.frame，包含gene_short_name列）
gene_annotation <- data.frame(
  gene_short_name = rownames(expr_matrix),
  row.names = rownames(expr_matrix))

cds <- new_cell_data_set(
  expression_data = expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)
# 从Seurat中提取UMAP坐标
umap_coords <- Embeddings(seurat_obj, reduction = "umap")

# 确保细胞顺序一致
all.equal(rownames(umap_coords), colnames(cds))  # 应为TRUE

# 将UMAP坐标添加到Monocle对象
reducedDims(cds)$UMAP <- umap_coords

# 可选：导入PCA坐标
pca_coords <- Embeddings(seurat_obj, reduction = "pca")[, 1:50]  # 保留前50个PC
reducedDims(cds)$PCA <- pca_coords

# 使用Seurat的聚类结果（假设元数据中有"seurat_clusters"列）
cds@clusters$UMAP <- list(
  clusters = seurat_obj$cell_type,
  partitions = rep(1, ncol(cds))  # 默认单个分区
)
# 学习轨迹图
cds <- cluster_cells(cds, resolution = 1e-5)  # 低分辨率避免过度分群
cds <- learn_graph(cds)
## 步骤3：定义轨迹起点

# 方法2：自动指定根节点（如选择特定簇）
root_cells <- colnames(cds)[cds$seurat_clusters == "1"]  # 假设簇1为起点
cds <- order_cells(cds, root_cells = root_cells)
pseudotime_values <- pseudotime(cds)
pbmc$pseudotime <- pseudotime_values  # 存回Seurat对象
#### Differentiation Trajectory  ####
library(ggplot2)
library(ggridges)

plot_data <- pbmc@meta.data 
# 1. (关键一步) 按照拟时序中位数，对细胞亚群进行排序，方便看出分化先后
plot_data$cell_type <- reorder(plot_data$cell_type, plot_data$pseudotime, median)

# 2. 画图命令
ggplot(plot_data, aes(x = pseudotime, y = cell_type, fill = cell_type)) +
  # 画脊线图
  geom_density_ridges(scale = 1.5, alpha = 0.8, color = "black") + 
  theme_classic() +
  labs(
    title = "Differentiation Trajectory of Subpopulations",
    x = "Pseudotime (Degree of Differentiation)",
    y = "Subpopulation"
  ) +
  # 隐藏图例（因为Y轴已经标了名字）
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"))


#### 寻找分支点的核心：graph_test
# 寻找在轨迹上显著变化的基因
gene_fits <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

# 筛选出显著基因
significant_genes <- row.names(subset(gene_fits, q_value < 0.05))
# 如果你的轨迹有分支，你可以指定 color_cells_by="branch" 来观察
cds_subset <- cds[c('Aspn' , 'Bglap', 'Dmp1' ,'Phex' ,"Kcnk2","Tnc",'Limch1','Enpp1',"Yap1", "Taz"), ]
p <- plot_genes_in_pseudotime(
  cds_subset,
  color_cells_by = "cell_type", # 按照细胞群上色
  min_expr = 0.5,               # Y轴最小截断值
  ncol = 4,                     # 设置为 2 列 (2x2 排列)
  cell_size = 1.5,              # 调整散点的大小
  trend_formula = "~ splines::ns(pseudotime, df=3)" # 生成黑色平滑拟合线
)

# 叠加 ggplot2 语法进行完美复刻
p_final <- p + 
  # 自定义颜色 (确保这里的名字与你的细胞分群名字完全一致)
  scale_color_manual(values = c(
    "LMP" = "#F8766D",        # 红色
    "Osteoblast" = "#7CAE00", # 绿色
    "Kcnk2+Cell" = "#00BFC4",        # 蓝绿色
    "Osteocyte" = "#C77CFF"   # 紫色
  )) +
  # 调整图例和排版细节
  theme(
    legend.position = "right",           # 图例放在右侧
    legend.title = element_blank(),      # 去除图例的标题
    strip.background = element_blank(),  # 去除每个子图上方基因名字的灰色背景框
    strip.text = element_text(size = 14, face = "bold"), # 放大基因名字
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# 打印最终图像
print(p_final)

plot_cells(cds, 
           color_cells_by = "cluster", 
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,       # 标记叶子节点（终点）
           label_branch_points=TRUE) # 标记分支点

# 将基因聚类成模块
gene_modules <- find_gene_modules(cds[significant_genes,], resolution=1e-2)

# 绘制模块在不同细胞簇中的表达
plot_cells(cds, genes=gene_modules, color_cells_by="cluster")


#### 细胞比例变化 ####
ggplot(pbmc@meta.data, aes(x = orig.ident, fill = cell_type)) +
  # position = "fill" 是自动计算百分比比例的核心
  geom_bar(position = "fill", width = 0.5, color = "black", size = 0.5) +
  # 转换为百分比显示
  scale_y_continuous(labels = scales::percent_format()) +
  # 使用我们一直统一的颜色系
  scale_fill_manual(values = c("LMP" = "#2CA02C",
                               "IMC" = "#D62728",
                               "Osteoblast" = "#1F77B4",
                               "Osteocyte" = "#9467BD")) +
  theme_classic() +
  labs(title = "Osteolineage Composition Shift",
       x = "Age",
       y = "Relative Proportion",
       fill = "Cell Type") +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )

#### 抗凋亡评分 ####
library(ggplot2)
library(dplyr)

# 抗凋亡基因集
anti_apoptotic_genes <- list(c("Bcl2", "Bcl2l1", "Bcl2l2", "Mcl1", "B16Mc5"))
# 基因集打分
pbmc <- AddModuleScore(pbmc, features = anti_apoptotic_genes, name = "Anti_Apoptotic_Score")

# 提取数据画 Split Violin Plot
library(ggplot2)
VlnPlot(pbmc, features = "Anti_Apoptotic_Score1", 
        group.by = "cell_type", 
        split.by = "orig.ident",
        split.plot = T,
        pt.size = 0) +
  scale_fill_manual(values = c("3M" = "#3498DB", "16M" = "#E74C3C")) +
  labs(title = "Anti-Apoptotic Signature in Aging", x = "Osteogenic Trajectory", y = "Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#### export for scvelo analysis ####
library(Seurat)
# 1. 提取 3M 样本的独立 Seurat 对象
# (假设你之前整合的大对象叫 pbmc，如果已经拆分请忽略这步)
pbmc_3M <- subset(pbmc, subset = orig.ident == "3M")
umap_coords <- Embeddings(pbmc_3M, reduction = "umap")
# 将坐标转为数据框并保存
write.csv(as.data.frame(umap_coords), 
          file = "3M_umap_coords.csv", 
          quote = FALSE)
meta_data <- pbmc_3M@meta.data
meta_data$celltype <- meta_data$cell_type
# 导出完整的 metadata
write.csv(meta_data, 
          file = "3M_metadata.csv", 
          quote = FALSE)
#### 16m Kcnk2 亚群差异基因  ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

# 1. 数据集准备与分组 (Data Preparation) 
# 确保使用 RNA assay 并合并 Layer (Seurat V5 标准操作)
DefaultAssay(pbmc) <- "RNA"
pbmc <- JoinLayers(pbmc)

# 提取 16M 的细胞 (假设你的年龄信息存在 orig.ident 或 age 列中)
# 如果 imc_subset 已经是纯 16M 的数据，这行可以跳过
seurat_16m <- subset(pbmc, subset = orig.ident == "16M")

# 创建一个新的 Metadata 列，将细胞严格划分为 Kcnk2_pos 和 Kcnk2_neg
# 假设亚群名称存在 cell_type 列中
seurat_16m$volcano_group <- ifelse(seurat_16m$cell_type == "IMC", "Kcnk2_pos", "Kcnk2_neg")

# 将分组设为默认身份
# 1. 确保在 RNA 模式下，并合并子集可能断开的图层
DefaultAssay(seurat_16m) <- "RNA"
seurat_16m <- JoinLayers(seurat_16m)

# 2. 【核心修复步】：执行基础 Log 标准化，强制生成完整的 'data' 图层！
seurat_16m <- NormalizeData(seurat_16m)

# 3. 再次运行找差异基因的代码（这次一定能顺利跑通）
deg_results <- FindMarkers(seurat_16m, 
                           ident.1 = "Kcnk2_pos", 
                           ident.2 = "Kcnk2_neg", 
                           logfc.threshold = 0.25, 
                           min.pct = 0.1)

# 查看一下算出来的结果（前几行）
head(deg_results)

# 将行名(基因名)提取为单独一列
deg_results$gene <- rownames(deg_results)

# 3. 阈值设定与颜色映射分配 
# 设定显著性阈值 (可根据实际 p_val_adj 范围微调)
p_val_adj_cutoff <- 0.05
logfc_cutoff <- 0.5 

deg_results <- deg_results %>%
  mutate(
    Significance = case_when(
      p_val_adj < p_val_adj_cutoff & avg_log2FC > logfc_cutoff ~ "Up in Kcnk2+",
      p_val_adj < p_val_adj_cutoff & avg_log2FC < -logfc_cutoff ~ "Down in Kcnk2+",
      TRUE ~ "Not Sig"
    )
  )

# 4. 指定需要“强制高亮”的灵魂靶点 
# 这里放上您故事链条里所有的核心枢纽！
genes_to_label <- c("Kcnk2", "Enpp1", "Dkk1", "Sost", 
                    "Trp53", "Cdkn1a", "Cdkn2a", "Serpine1", "Il6")

# 提取这些靶点的数据用于加标签
label_data <- deg_results %>% filter(gene %in% genes_to_label)

# 5. 绘制顶刊级别火山图 
ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  # 绘制背景点
  geom_point(alpha = 0.7, size = 1.5) +
  
  # 使用经典的“警示红”代表病理高表达，“清澈蓝”代表下调
  scale_color_manual(values = c("Up in Kcnk2+" = "#E74C3C", 
                                "Down in Kcnk2+" = "#3498DB", 
                                "Not Sig" = "gray80")) +
  
  # 添加坐标轴的十字辅助虚线
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(p_val_adj_cutoff), linetype = "dashed", color = "black", linewidth = 0.4) +
  
  # 利用 ggrepel 添加不重叠的基因标签 (极其重要的一步)
  geom_text_repel(data = label_data,
                  aes(label = gene),
                  size = 4.5,
                  fontface = "bold.italic", # 基因名标准斜体加粗
                  box.padding = 0.8,
                  point.padding = 0.3,
                  segment.color = "black",
                  segment.size = 0.5,
                  color = "black",
                  max.overlaps = Inf,
                  nudge_x = 0.5, # 稍微向右推一点标签，更具呼吸感
                  nudge_y = 0.5) +
  
  # 极致排版主题
  theme_classic() +
  labs(title = "Pathological Signature of Kcnk2+ vs Kcnk2- Cells",
       subtitle = "16M Aged Niche: Enrichment of p53-driven SASP & Wnt Antagonists",
       x = expression(bold("Average " * log[2] * "(Fold Change)")),
       y = expression(bold("-log"[10] * "(Adjusted P-value)"))) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, face = "italic", hjust = 0.5, color = "gray30"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    # 加粗坐标轴线
    axis.line = element_line(linewidth = 0.8)
  )


##### 富集“蓝色”基因 ####
# 0. 加载富集分析专用包 
library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠基因注释数据库
library(enrichplot)
library(ggplot2)
library(dplyr)

# 1. 提取显著下调的“蓝色”基因 
# 提取 P<0.05 且 LogFC < -0.5 的基因
down_genes <- deg_results %>%
  filter(p_val_adj < 0.05 & avg_log2FC < -0.5) %>%
  pull(gene)

# 2. 基因 ID 转换 (Symbol 转换为 Entrez ID)
# clusterProfiler 在做 KEGG 时严格需要 Entrez ID 
gene_ids <- bitr(down_genes, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Mm.eg.db)

# 3. GO 富集分析 (BP: Biological Process 生物学过程) 
go_down_results <- enrichGO(gene          = gene_ids$ENTREZID,
                            OrgDb         = org.Mm.eg.db,
                            ont           = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.2,
                            readable      = TRUE) # 将结果中的 ID 转回基因名方便阅读

# 4. KEGG 通路富集分析 
kegg_down_results <- enrichKEGG(gene          = gene_ids$ENTREZID,
                                organism      = 'mmu', # mmu 代表小鼠 (Mus musculus)
                                pvalueCutoff  = 0.05)

# 5. 顶刊级别可视化：GO 富集气泡图 (Dotplot) 
dotplot(go_down_results, showCategory = 10, 
        title = "Suppressed Biological Processes in Kcnk2+ Cells") +
  theme_classic() +
  scale_color_gradientn(colors = c("#E74C3C", "#F1C40F", "#3498DB")) + # 调整 p 值渐变色
  theme(
    axis.text.y = element_text(size = 12, face = "bold", color = "black"), # 突出显示通路名称
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

dotplot(kegg_down_results, showCategory = 10, 
        title = "Suppressed Biological Processes in Kcnk2+ Cells") +
  theme_classic() +
  scale_color_gradientn(colors = c("#E74C3C", "#F1C40F", "#3498DB")) + # 调整 p 值渐变色
  theme(
    axis.text.y = element_text(size = 12, face = "bold", color = "black"), # 突出显示通路名称
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "right"
  )
# 确保 rownames 被转化为一列专门的 Gene 名称，防止导出时丢失
deg_export <- deg_results
deg_export$Gene <- rownames(deg_export)

# 调整列的顺序，把 Gene 放在第一列，看起来更专业
deg_export <- deg_export[, c("Gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "Significance")]

# 导出为 CSV 文件，可以直接用 Excel 打开并作为 Supplementary Table 提交
write.csv(deg_export, file = "/Users/kenny/Desktop/2025/运动/sub/re-sub/Table/Supplementary_Table_Kcnk2_DEGs.csv", row.names = FALSE, quote = FALSE)
