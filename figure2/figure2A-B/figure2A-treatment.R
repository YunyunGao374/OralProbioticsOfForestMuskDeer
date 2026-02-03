# -------------------------------
# NMDS 分析与可视化（按 Treat 分组）
# -------------------------------

library(vegan)
library(ggplot2)
library(ggpubr)

# 读取 OTU 表格
otu <- read.table("otu_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(otu) <- otu$NAME
otu <- otu[, -1]  # 删除 OTU 名称列

# 转置 OTU 表格
otu_t <- t(otu)

# 读取样本元数据
metadata <- read.table("sample_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample

# 计算 Bray-Curtis 距离矩阵
dist_matrix <- vegdist(otu_t, method = "bray")

# 执行 NMDS（二维）
set.seed(123)
nmds <- metaMDS(dist_matrix, k = 2, trymax = 100)
nmds$stress  # 查看 stress 值

# 提取 NMDS 坐标并合并 metadata
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$sample <- rownames(nmds_scores)
nmds_scores <- merge(nmds_scores, metadata, by = "sample")

# 计算 NMDS 拟合优度 R²
fit_dist <- as.dist(vegdist(scores(nmds, display = "sites"), method = "euclidean"))
cor_val <- cor(dist_matrix, fit_dist, method = "spearman")
R2 <- cor_val^2
R2

# PERMANOVA 检验 Treat 分组显著性
set.seed(123)
perm <- adonis2(dist_matrix ~ Treat, data = metadata, permutations = 999)
p_value <- perm$`Pr(>F)`[1]
p_value

# 设置 Treat 分组颜色和形状
col_groups <- c("Control" = "#426aa0", "Treatment" = "#c13328")
shape_groups <- c("Control" = 16, "Treatment" = 17)

# 1. NMDS 散点图（Treat 分组）
nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Treat, shape = Treat)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.1, aes(fill = Treat), color = NA) +
  scale_color_manual(values = col_groups) +
  scale_fill_manual(values = col_groups) +
  scale_shape_manual(values = shape_groups) +
  theme_minimal() +
  labs(title = paste("NMDS plot (Stress =", round(nmds$stress, 3),
                     ", R² =", round(R2, 3),
                     ", PERMANOVA p =", signif(p_value, 3), ")"),
       x = "NMDS1", y = "NMDS2") +
  theme(text = element_text(size = 12))

# 2. NMDS1 横向箱线图 + 两两比较 p 值
# 只有两组 Treat，比较就是 Control vs Treatment
comparisons <- list(c("Control", "Treatment"))

nmds1_box <- ggplot(nmds_scores, aes(x = Treat, y = NMDS1, fill = Treat, color = Treat)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.shape = NA, size = 0.8) +  # 半透明箱体 + 边框
  geom_jitter(width = 0.2, alpha = 0.5) +  # 散点半透明
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format") +  # 显示真实 p 值
  coord_flip() +  # 横向显示
  scale_fill_manual(values = scales::alpha(col_groups, 0.5)) +
  scale_color_manual(values = col_groups) +
  theme_minimal() +
  labs(x = "Treat", y = "NMDS1") +
  theme(text = element_text(size = 12), legend.position = "none")

# 3. 合并 NMDS 散点图和箱线图
combined_plot <- ggarrange(nmds_plot, nmds1_box,
                           ncol = 1, nrow = 2, heights = c(3, 1))  # 散点图高度比箱线图大

# 4. 输出合并图
print(combined_plot)
