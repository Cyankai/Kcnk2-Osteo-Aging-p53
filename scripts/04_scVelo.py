import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np

# 1. 全局设置：强制输出 SVG 矢量图（完美避开 PDF 的报错，且支持 AI 无限放大修改）
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', dpi_save=300, transparent=True, format='svg')

# 2. 读取数据并修复 Barcode
adata = sc.read_loom("3M.loom")
adata.obs_names = [name.split(':')[-1].replace('x', '_2') for name in adata.obs_names]

# 3. 导入 R 语言做好的精美 UMAP 坐标和细胞注释
umap_df = pd.read_csv("3M_umap_coords.csv", index_col=0)
meta_df = pd.read_csv("3M_metadata.csv", index_col=0)

# 4. 提取交集细胞并赋值
common_cells = list(set(adata.obs_names) & set(meta_df.index))
adata = adata[common_cells].copy()
adata.obsm['X_umap'] = umap_df.loc[adata.obs_names].values
adata.obs['cell_clusters'] = meta_df.loc[adata.obs_names, 'cell_type'].astype('category')

# 5. 最新版无错预处理流程
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata, enforce=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)

# 6. 核心动力学计算
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='deterministic')
scv.tl.velocity_graph(adata, show_progress_bar=False)

# 7. 终极排雷补丁：强制清除可能导致画图报错的 NaN 坏点
scv.tl.velocity_embedding(adata, basis='umap')
adata.obsm['velocity_umap'] = np.nan_to_num(adata.obsm['velocity_umap'])

# 8. 可视化
scv.pl.velocity_embedding_stream(
    adata, 
    basis='umap', 
    color='cell_clusters', 
    legend_loc='right margin', 
    title='3M Bone Lineage Velocity',
    density=0.5,      # 【核心】降低流线密度，抹去杂音，凸显主干道
    linewidth=2.0,    # 【核心】线条加粗，增加视觉冲击力
    arrow_size=2.0,   # 【核心】箭头放大，让方向感更清晰
    smooth=True,      # 轨迹平滑化处理
    figsize=(10, 8),  # 稍微放大画布，让细胞群不那么拥挤
    save='3M_velocity_stream.svg' # 保存为 SVG
)
