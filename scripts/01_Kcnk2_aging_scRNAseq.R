#### 16M clusters annotation ####
library(Seurat)
library(ggplot2)
library(dplyr)
options(future.globals.maxSize = 16 * 1024^3)
# 1. 创建样本列表（假设已加载3个样本）
pbmc1=as.data.frame(read.table('/Users/kenny/Documents/seq/sc_seq/data/Qin lin/1M/GSM4318799_1M_matrix.txt', row.names = 1, header = T, sep = " "))
pbmc2=as.data.frame(read.table('/Users/kenny/Documents/seq/sc_seq/data/Qin lin/3M/GSM4318801_3M_matrix.txt', row.names = 1, header = T, sep = " "))
pbmc3=as.data.frame(read.table('/Users/kenny/Documents/seq/sc_seq/data/Qin lin/16M/GSM4318802_16M_matrix.txt', row.names = 1, header = T, sep = " "))
pbmc4=Read10X(data.dir = '/Users/kenny/Desktop/Plpp1/P 骨/scRNA-seq/整体results/matrix/WT/')
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
res.used <- c(0.01,0.04,0.07)
integrated <- FindClusters(integrated, resolution = res.used)

# Make plot 
library(clustree)
clustree(integrated@meta.data, prefix =
           "integrated_snn_res."
) 

# 确定resolution
final_resolution = 0.07
integrated <- FindClusters(integrated, resolution = final_resolution)

# 确定成骨亚群
Chondrocyte = c('Sox9','Col2a1','Col10a1','Pth1r','Acan','Ihh') 
obot=c("Col1a1","Bglap","Sp7","Alpl","Pdpn","Dmp1")
VlnPlot(integrated, features = obot)  #  2\5\6\10 号亚群
VlnPlot(integrated, features = Chondrocyte)  #  6号亚群

# 6. 可视化批次效应
DimPlot(integrated, label = T, pt.size = 0.2,label.size = 4,
        #  group.by = 'cell_type',
        split.by= 'orig.ident' )

# Ob亚群
pbmc <- subset(integrated, subset = seurat_clusters %in% c(2,5,10) )
pbmc <- subset(pbmc, subset =  orig.ident %in% c("16M", "3M"))
pbmc=SCTransform(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10,reduction='pca')
res.used <- c(0.01,0.04,0.07,.09)
pbmc <- FindClusters(pbmc, resolution = res.used)
library(clustree)
clustree(pbmc@meta.data, prefix =
           # "RNA_snn_res.",
           "SCT_snn_res."
         #'integrated_snn_res.'
) 
final_resolution = 0.04
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
cell_type = c("Kcnk2+Cell",'LMP',"Osteoblast","Osteocyte")
pbmc[['cell_type']] = unname(cell_type[pbmc@meta.data$seurat_clusters])
pbmc$cell_type <- factor(pbmc$cell_type, levels = c("LMP", "Osteoblast", "Kcnk2+Cell", "Osteocyte"))
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

save(pbmc,data,tpms,file = '/Users/kenny/Desktop/pbmc_yours.RData')

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
  labs(title = "Total Expression of Wnt Antagonists (3M vs 16M)",
       subtitle = "Floating annotations indicate the absolute age-induced increment (Δ)",
       x = "Osteogenic Lineage", y = "Total Transcript Expression", fill = "Age Group") +
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

#### Global correlation ####
# 1. 载入必要的 R 包
library(Seurat)
library(pheatmap) 

# 创建“亚群-年龄”组合标签
pbmc$celltype_age <- paste(pbmc$cell_type, pbmc$orig.ident, sep = "-")
Idents(pbmc) <- "celltype_age"
# 2. 计算各亚群的全局平均表达量 (Global transcriptomic profile)
avg_expression_list <- AverageExpression(pbmc, 
                                         assays = "RNA", 
                                         layer = "data", 
                                         return.seurat = FALSE)
# 提取 RNA assay 的表达矩阵 (行是基因，列是你的四个 cell_type 亚群)
expr_matrix <- as.matrix(avg_expression_list$RNA)
# 3. 计算全局转录组的 Spearman 相关性
spearman_cor_matrix <- cor(expr_matrix, method = "spearman")
# 4. (可选) 将相关性矩阵可视化为热图
pheatmap(spearman_cor_matrix, 
         display_numbers = TRUE,          # 在热图上显示具体的数字
         fontsize_number = 12,            # 数字大小
         color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100), # 经典红蓝配色
         main = "Spearman Correlation of 3M Subpopulations")


#### Monocle3 ####
library(Seurat) 
library(monocle3)
library(viridis)
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
root_cells <- colnames(cds)[cds$seurat_clusters == "1"]  # 假设簇0为起点
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
cds_subset <- cds[c("Aspn", "Kcnk2", "Bglap", "Dmp1"), ]
plot_genes_in_pseudotime(cds_subset, 
                         color_cells_by="partition", 
                         min_expr=0.5)
# 在 Monocle 3 中，判断是否偏离最直观的方法是看 Principal Graph（主图结构）。
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
                               "Kcnk2+Cell" = "#D62728",
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
anti_apoptotic_genes <- list(c("Bcl2", "Bcl2l1", "Bcl2l2", "Mcl1", "Birc5"))
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
#### Pathway_Scores_All_Cells ####
pathways_list <- list(
  CaN_Pathway  = c("Ppp3ca", "Ppp3cb", "Ppp3r1", "Nfatc1", "Nfatc2"),
  CaMK_Pathway = c("Camk1", "Camk2a", "Camk2b", "Camk2d", "Camk2g", "Camkk1", "Camkk2"),
  PKC_Pathway  = c("Prkca", "Prkcb", "Prkcg", "Prkcd", "Prkce", "Prkcz"),
  AMPK_Pathway = c("Prkaa1", "Prkaa2", "Prkab1", "Prkab2", "Prkag1")
)

# 过滤掉矩阵中不存在的基因以防报错
valid_pathways <- lapply(pathways_list, function(genes) intersect(genes, rownames(pbmc)))
pbmc <- AddModuleScore(pbmc, features = valid_pathways, name = "Target_Score_")

# Seurat 默认会在名字后加 1, 2, 3, 4，我们将其重命名为直观的名称
score_names <- c("CaN_Score", "CaMK_Score", "PKC_Score", "AMPK_Score")
colnames(pbmc@meta.data)[(ncol(pbmc@meta.data) - 3) : ncol(pbmc@meta.data)] <- score_names
# 建议在全数据集 pbmc 上直接创建新列，方便全局比较
pbmc$Analysis_Group <- "Others" # 初始化

# 定义逻辑：结合年龄 (orig.ident) 和 细胞类型 (cell_type)
pbmc$Analysis_Group[pbmc$orig.ident == "3M" & pbmc$cell_type == "Kcnk2+Cell"] <- "3M_Kcnk2_Pos"
pbmc$Analysis_Group[pbmc$orig.ident == "3M" & pbmc$cell_type != "Kcnk2+Cell"] <- "3M_Others"
pbmc$Analysis_Group[pbmc$orig.ident == "16M" & pbmc$cell_type == "Kcnk2+Cell"] <- "16M_Kcnk2_Pos"
pbmc$Analysis_Group[pbmc$orig.ident == "16M" & pbmc$cell_type != "Kcnk2+Cell"] <- "16M_Others"

# 设置因子顺序
pbmc$Analysis_Group <- factor(pbmc$Analysis_Group, 
                              levels = c("3M_Others", "3M_Kcnk2_Pos", "16M_Others", "16M_Kcnk2_Pos"))
# 1. 明确我们需要提取的列名
# 包括：基础分组信息、Kcnk2表达量、以及四个通路的得分
export_vars <- c("orig.ident",       # 年龄分组 (3M, 16M)
                 "cell_type",        # 细胞类型亚群 (IMC, Kcnk2+Cell 等)
                 "Analysis_Group",   # 我们刚才创建的 4 分组 (3M_Others, 3M_Kcnk2_Pos, 16M_Others, 16M_Kcnk2_Pos)
                 "Kcnk2",            # Kcnk2 的原始表达量 (方便你在 GraphPad 里随时核对)
                 "CaN_Score", 
                 "CaMK_Score", 
                 "PKC_Score", 
                 "AMPK_Score")

# 2. 使用 FetchData 提取这些列的所有细胞数据
# FetchData 是 Seurat 提取数据的最安全方法，能确保基因表达和 metadata 完美对齐
export_data <- FetchData(pbmc, vars = export_vars)

# 3. (可选) 将行名（细胞的 Barcode）提取出来变成单独的一列
export_data$Cell_Barcode <- rownames(export_data)
# 调整列的顺序，把 Cell_Barcode 放到第一列
export_data <- export_data[, c("Cell_Barcode", export_vars)]

# 4. 导出为 CSV 文件
# row.names = FALSE 因为我们已经把行名提取成了 Cell_Barcode 列
write.csv(export_data, file = "/Users/kenny/Desktop/2025/运动/scRNA_seq/Qinlin_sc_data/
          Pathway_Scores_All_Cells.csv", row.names = FALSE)
