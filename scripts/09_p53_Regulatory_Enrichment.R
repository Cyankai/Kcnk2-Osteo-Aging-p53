library(dplyr)
library(ggplot2)
library(tibble)

# Part 1. 紧凑型谱系棒棒糖图  #####
# 1. 提取 p53 CUT&Tag 在核心启动子区的结合信号 (忽略大小写以防漏滤)
p53_promoter_anno <- as.data.frame(peakAnno_Nut_anno) %>%
  filter(grepl("Promoter", annotation, ignore.case = TRUE) | abs(distanceToTSS) <= 3000) %>%
  group_by(SYMBOL) %>%
  summarise(p53_binding_score = max(V7, na.rm = TRUE), .groups = "drop")

# 2. 定义四大谱系代表性核心 Identity Markers 
core_identity_markers <- tribble(
  ~SYMBOL,    ~Lineage,
  # LMP 前体特征
  "Aspn",     "LMP (Progenitor)",
  "Edil3",    "LMP (Progenitor)",
  "Tnn",      "LMP (Progenitor)",
  # OB 成熟成骨特征
  "Bglap",    "OB (Osteoblast)",
  "Sp7",      "OB (Osteoblast)",
  "Ifitm5",   "OB (Osteoblast)",
  # IMC 停滞/力学钝化特征 (核心驱动轴)
  "Limch1",   "IMC (Intermediate)",
  "Kcnk2",    "IMC (Intermediate)",
  "Enpp1",    "IMC (Intermediate)",
  "Tnc",      "IMC (Intermediate)",
  # Ot 终末骨细胞成熟特征
  "Dmp1",     "Ot (Osteocyte)",
  "Phex",     "Ot (Osteocyte)", 
  "Fgf23",    "Ot (Osteocyte)"
)

# 动态提取各 Marker 的富集强度
get_identity_fc <- function(gene_name, lineage) {
  marker_df <- switch(lineage,
                      "LMP (Progenitor)"   = markers_LMP,
                      "OB (Osteoblast)"    = markers_OB,
                      "IMC (Intermediate)" = markers_IMC,
                      "Ot (Osteocyte)"     = markers_Ot
  )
  if (!is.null(marker_df) && gene_name %in% rownames(marker_df)) {
    return(marker_df[gene_name, "avg_log2FC"])
  } else {
    return(1.0)
  }
}

df_identity_p53 <- core_identity_markers %>%
  rowwise() %>%
  mutate(Identity_Log2FC = get_identity_fc(SYMBOL, Lineage)) %>%
  ungroup() %>%
  left_join(p53_promoter_anno, by = "SYMBOL") %>%
  mutate(
    p53_binding_score = coalesce(p53_binding_score, 0),
    is_direct_target = ifelse(p53_binding_score >= 4, "Direct Target (Promoter Bound)", "Non-Target / Passive"),
    Lineage = factor(Lineage, levels = c("LMP (Progenitor)", "OB (Osteoblast)", "IMC (Intermediate)", "Ot (Osteocyte)"))
  )

# 指定水平轴排列顺序
gene_horizontal_order <- c(
  "Aspn", "Edil3", "Tnn",
  "Bglap", "Sp7", "Ifitm5",
  "Limch1", "Kcnk2", "Enpp1", "Tnc",
  "Dmp1", "Phex", "Fgf23"
)

df_plot_xy <- df_identity_p53 %>%
  filter(SYMBOL %in% gene_horizontal_order) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = gene_horizontal_order))

# 绘制棒棒糖图
p_identity_compact <- ggplot(df_plot_xy, aes(x = SYMBOL, y = Identity_Log2FC)) +
  geom_segment(aes(xend = SYMBOL, y = 0, yend = Identity_Log2FC), 
               color = "grey60", linewidth = 0.6) +
  geom_point(aes(size = p53_binding_score, fill = is_direct_target, color = is_direct_target), 
             shape = 21, stroke = 0.8) +
  scale_fill_manual(values = c("Direct Target (Promoter Bound)" = "#D9383A", 
                               "Non-Target / Passive"          = "#4E79A7")) +
  scale_color_manual(values = c("Direct Target (Promoter Bound)" = "#990000", 
                                "Non-Target / Passive"          = "#1C4E80")) +
  scale_size_continuous(range = c(2.5, 6.5), breaks = c(0, 2, 4), 
                        name = "p53 CUT&Tag\nScore") +
  facet_grid(. ~ Lineage, scales = "free_x", space = "free_x") +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "Cluster Specificity (Marker avg_log2FC)",
    title = "Lineage-Specific Identity Markers Under Direct p53 Control"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "#EFEFEF", color = NA),
    axis.text.x = element_text(face = "bold.italic", size = 9, angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(5, 5, 5, 12)
  )

print(p_identity_compact)
 
# Part 2. 超几何富集检验柱状图  #####
# 1. 背景基因池与启动子结合基因
bg_genes <- rownames(pbmc)
N_total  <- length(bg_genes)

anno_df <- as.data.frame(peakAnno_Nut_anno)
p53_promoter_genes <- anno_df %>%
  filter(grepl("Promoter", annotation, ignore.case = TRUE) | abs(distanceToTSS) <= 3000) %>%
  pull(SYMBOL) %>%
  unique()
p53_promoter_genes <- intersect(p53_promoter_genes, bg_genes)
M_bound <- length(p53_promoter_genes)

# 2. 超几何检验函数
calc_hyper_enrichment <- function(marker_genes, cluster_name) {
  markers_in_bg <- intersect(marker_genes, bg_genes)
  k_marker <- length(markers_in_bg)
  overlap_genes <- intersect(markers_in_bg, p53_promoter_genes)
  q_overlap <- length(overlap_genes)
  
  p_val <- phyper(q_overlap - 1, m = M_bound, n = N_total - M_bound, k = k_marker, lower.tail = FALSE)
  expected <- (M_bound / N_total) * k_marker
  fold_enrichment <- q_overlap / expected
  
  tibble(
    Cluster = cluster_name,
    Marker_Total = k_marker,
    Overlap = q_overlap,
    Expected = round(expected, 1),
    Fold_Enrichment = round(fold_enrichment, 2),
    P_value = p_val
  )
}

# 3. 提取 OB / IMC / Ot 的 Marker 进行检验
get_sig_markers <- function(df) {
  if (!"SYMBOL" %in% colnames(df)) df$SYMBOL <- rownames(df)
  df %>% filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>% pull(SYMBOL) %>% unique()
}

markers_list <- list( 
  "OB"  = get_sig_markers(markers_OB),
  "IMC" = get_sig_markers(markers_IMC),
  "Ot"  = get_sig_markers(markers_Ot)
)

enrich_results <- bind_rows( 
  calc_hyper_enrichment(markers_list$OB,  "OB (Osteoblast)"),
  calc_hyper_enrichment(markers_list$IMC, "IMC (Intermediate)"),
  calc_hyper_enrichment(markers_list$Ot,  "Ot (Osteocyte)")
) %>%
  mutate(
    neg_log10_P = -log10(P_value),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# 4. 绘制柱状图
enrich_results$Cluster <- factor(enrich_results$Cluster, 
                                 levels = c("OB (Osteoblast)", "IMC (Intermediate)", "Ot (Osteocyte)"))

p_enrich <- ggplot(enrich_results, aes(x = Cluster, y = -log10(P_value), fill = Cluster)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.4, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text(aes(label = paste0(Fold_Enrichment, "x\n(", Significance, ")")), 
            vjust = -0.3, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = c( 
    "OB (Osteoblast)"      = "#F28E2B",
    "IMC (Intermediate)"   = "#D9383A",
    "Ot (Osteocyte)"        = "#76B7B2"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = expression(-log[10](italic(P)~value)),
    title = "p53 Regulatory Enrichment across Lineage Markers",
    subtitle = "Hypergeometric test vs genome-wide background"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(face = "bold", size = 10, color = "black"),
    axis.text.y = element_text(color = "black")
  )

print(p_enrich)

# Part 3. 导出 Source Data (CSV 格式) #####
# 自动创建目标目录（若不存在则自动新建，防止写文件报错）
target_dir <- "/users/kenny/Desktop/2025/运动/sub/re-sub/Table" 

df_source <- enrich_results %>%
  filter(Cluster %in% c("OB (Osteoblast)", "IMC (Intermediate)", "Ot (Osteocyte)")) %>%
  select(
    Lineage_Cluster = Cluster,
    Total_Cluster_Markers = Marker_Total,
    Observed_Overlap_with_p53 = Overlap,
    Expected_Overlap_Background = Expected,
    Fold_Enrichment = Fold_Enrichment,
    P_Value_Hypergeometric = P_value,
    Neg_Log10_P_Value = neg_log10_P,
    Significance_Flag = Significance
  )

df_source <- df_plot_xy %>%
  select(
    Lineage_Facet = Lineage,
    Gene_Symbol = SYMBOL,
    Cluster_Specificity_Marker_avg_log2FC = Identity_Log2FC,
    p53_CUTnTag_Binding_Score_V7 = p53_binding_score,
    Target_Classification = is_direct_target
  )

write.csv(df_source , file = file.path(target_dir, "Source_Data_p53_Regulatory_Enrichment.csv"), row.names = FALSE)
write.csv(df_source , file = file.path(target_dir, "Source_Data_p53_Regulatory_Lollipop.csv"), row.names = FALSE)
