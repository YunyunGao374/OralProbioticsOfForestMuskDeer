# -------------------------------
# Alpha 多样性分析 + 小提琴图（含 p 值）
# -------------------------------

# 1. 加载所需包
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)

# 2. 读取 OTU 表格
otu <- read.table("otu_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(otu) <- otu$NAME
otu <- otu[, -1]  # 删除 OTU 名称列

# 3. 转置 OTU 表格（样本为行）
otu_t <- t(otu)

# 4. 读取元数据（包含分组信息）
metadata <- read.table("sample_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample

# 5. 计算 alpha 多样性指数
alpha_div <- data.frame(
  sample = rownames(otu_t),
  Shannon = diversity(otu_t, index = "shannon"),
  Simpson = diversity(otu_t, index = "simpson"),
  Richness = specnumber(otu_t),
  Chao1 = estimateR(otu_t)["S.chao1", ]  # Chao1
)

# 6. 合并 metadata
alpha_div <- merge(alpha_div, metadata, by = "sample")

# 7. 设置比较组（Periods）
group_col <- "Periods" 
groups <- unique(alpha_div[[group_col]])
comparisons <- combn(groups, 2, simplify = FALSE)  # 两两比较

# 8. 定义画小提琴图的函数
plot_violin <- function(df, index, group_col, comparisons, colors = NULL){
  p <- ggplot(df, aes_string(x = group_col, y = index, fill = group_col, color = group_col)) +
    geom_violin(alpha = 0.5, trim = FALSE) +      # 半透明小提琴图
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) + # 散点
    stat_summary(fun = mean, geom = "crossbar", 
                 width = 0.3, color = col_groups , size = 0.6, fatten = 1) + # 均值线
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format") +
    theme_minimal() +
    labs(x = group_col, y = index) +
    theme(text = element_text(size = 12), legend.position = "none")
  
  if(!is.null(colors)){
    p <- p + scale_fill_manual(values = scales::alpha(colors, 0.5)) +
      scale_color_manual(values = colors)
  }
  return(p)
}
# 9. 自定义颜色
col_groups <- c("P1" = "#7f95d3",
                "P2" = "#f4a8a8",
                "P3" = "#68ad5c",
                "P4" = "#994d99")

# 10. 绘制各 alpha 指数小提琴图
p_shannon <- plot_violin(alpha_div, "Shannon", group_col, comparisons, col_groups)
p_simpson <- plot_violin(alpha_div, "Simpson", group_col, comparisons, col_groups)
p_richness <- plot_violin(alpha_div, "Richness", group_col, comparisons, col_groups)
p_chao1 <- plot_violin(alpha_div, "Chao1", group_col, comparisons, col_groups)

# 11. 合并显示（2x2）
combined_plot <- ggarrange(p_shannon, p_simpson, p_richness, p_chao1,
                           ncol = 2, nrow = 2)

# 12. 输出图
print(combined_plot)
