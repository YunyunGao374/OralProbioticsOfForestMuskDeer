# -------------------------------
# NMDS 分析与可视化（按 Periods 分组）
# -------------------------------

# 1. 加载所需包
library(vegan)
library(ggplot2)
library(ggpubr)  #
# 2. 读取 OTU 表格
otu <- read.table("otu_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(otu) <- otu$NAME
otu <- otu[, -1]  # 删除 OTU 名称列

# 3. 转置 OTU 表格，样本为行
otu_t <- t(otu)

# 4. 读取样本元数据
metadata <- read.table("sample_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample

# 5. 计算 Bray-Curtis 距离矩阵
dist_matrix <- vegdist(otu_t, method = "bray")

# 6. 执行 NMDS（二维）
set.seed(123)
nmds <- metaMDS(dist_matrix, k = 2, trymax = 100)
nmds$stress  # 查看 stress 值

# 7. 提取 NMDS 坐标并合并 metadata
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$sample <- rownames(nmds_scores)
nmds_scores <- merge(nmds_scores, metadata, by = "sample")

# 8. 计算 NMDS 拟合优度 R²
fit_dist <- as.dist(vegdist(scores(nmds, display = "sites"), method = "euclidean"))
cor_val <- cor(dist_matrix, fit_dist, method = "spearman")
R2 <- cor_val^2
R2  # 输出 R²

# 9. PERMANOVA 检验 Periods 分组显著性
set.seed(123)
perm <- adonis2(dist_matrix ~ Periods, data = metadata, permutations = 999)
p_value <- perm$`Pr(>F)`[1]
p_value  # 输出 p 值

# 10. NMDS 可视化（按 Periods 分组上色，并显示置信椭圆）
col_groups <- c( "P1"        = "#7f95d3",
                 "P2"       = "#f4a8a8",
                 "P3" = "#68ad5c",
                 "P4"   = "#994d99")

# NMDS 可视化
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Periods, shape = Periods)) +
  geom_point(size = 3) +
  # 置信椭圆，只修改透明度，不修改点的颜色和形状
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.1, aes(fill = Periods), color = NA) +
  scale_color_manual(values = col_groups) +
  scale_fill_manual(values = col_groups) +
  theme_minimal() +
  labs(title = paste("NMDS plot (Stress =", round(nmds$stress, 3),
                     ", R² =", round(R2, 3),
                     ", PERMANOVA p =", signif(p_value, 3), ")"),
       x = "NMDS1", y = "NMDS2") +
  theme(text = element_text(size = 12))


nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Periods, shape = Periods)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.1, aes(fill = Periods), color = NA) +
  scale_color_manual(values = col_groups) +
  scale_fill_manual(values = col_groups) +
  theme_minimal() +
  labs(title = paste("NMDS plot (Stress =", round(nmds$stress, 3),
                     ", R² =", round(R2, 3),
                     ", PERMANOVA p =", signif(p_value, 3), ")"),
       x = "NMDS1", y = "NMDS2") +
  theme(text = element_text(size = 12))

# 2. NMDS1 箱线图 + 两两比较 p 值
# 设置两两比较
comparisons <- list(c("P1", "P2"), c("P1", "P3"), c("P1", "P4"),
                    c("P2", "P3"), c("P2", "P4"), c("P3", "P4"))

nmds1_box <- ggplot(nmds_scores, aes(x = Periods, y = NMDS1, fill = Periods, color = Periods)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.shape = NA, size = 0.8) + # 透明箱体 + 边框加粗
  geom_jitter(width = 0.2, alpha = 0.5) + # 添加散点
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format") +
  coord_flip() + # 横向显示
  scale_fill_manual(values = scales::alpha(col_groups, 0.5)) +
  scale_color_manual(values = col_groups) +
  theme_minimal() +
  labs(x = "Periods", y = "NMDS1") +
  theme(text = element_text(size = 12), legend.position = "none")

nmds1_box


# 3. 将 NMDS 散点图与箱线图合并（垂直排列）
combined_plot <- ggarrange(nmds_plot, nmds1_box,
                           ncol = 1, nrow = 2, heights = c(3, 1))  # 散点图高度比箱线图大

# 4. 输出合并图
print(combined_plot)
