library(tidyverse)
library(ggbreak)
# =========================
# 0. 颜色与顺序设置
# =========================
species_colors <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072",
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9"
)

times_bg <- c(
  "Hiding_phase"        = "#bfc8e6",
  "Conlict_phase"       = "#fbe5e4",
  "Disintegrated_phase" = "#d2e8cf",
  "Independent_phase"   = "#e2cae1"
)

times_order <- c(
  "Hiding_phase",
  "Conlict_phase",
  "Disintegrated_phase",
  "Independent_phase"
)
species_order <- c(
  "Eimeria_aquae",
  "Eimeria_dolichocystis",
  "Eimeria_fengxianensis",
  "Eimeria_helini",
  "Eimeria_jinfengshanenisis",
  "Eimeria_kaii",
  "Eimeria_oocylindrica",
  "Eimeria_sp.1",
  "Eimeria_sp.2"
)

# =========================
# 1. 读入数据
# =========================
df <- read.table(
  "parasite_all.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# =========================
# 2. 过滤：total > 100 + Treatment 组
# =========================
df_filt <- df %>%
  filter(total > 100, Groups == "Treatment")

# =========================
# 3. 物种列
# =========================
species_cols <- names(df_filt)[5:ncol(df_filt)]

# =========================
# 4. 宽表转长表
# =========================
df_long <- df_filt %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "Species",
    values_to = "Count"
  ) %>%
  filter(Count > 0)

# 固定物种顺序，防止颜色乱跳
df_long$Species <- factor(df_long$Species, levels = unique(df_long$Species))
df_long$Species <- factor(df_long$Species, levels = species_order)
# =========================
# 5. 计算百分比
# =========================
df_pct <- df_long %>%
  group_by(Samples) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# 固定 Times 顺序
df_pct$Times <- factor(df_pct$Times, levels = times_order)
df_pct$Species <- factor(df_pct$Species, levels = species_order)
df_long$Times <- factor(df_long$Times, levels = times_order)

# =========================
# 6. 让 Samples 按 Times 连续排列
# =========================
df_pct <- df_pct %>%
  arrange(Times, Samples)

df_long <- df_long %>%
  arrange(Times, Samples)

df_pct$Samples  <- factor(df_pct$Samples,  levels = unique(df_pct$Samples))
df_long$Samples <- factor(df_long$Samples, levels = unique(df_long$Samples))

# =========================
# 7. 计算背景色数据
# =========================
bg_df <- df_pct %>%
  distinct(Samples, Times) %>%
  group_by(Times) %>%
  summarise(
    xmin = min(as.numeric(Samples)) - 0.5,
    xmax = max(as.numeric(Samples)) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(
    bg   = times_bg[as.character(Times)],
    xmid = (xmin + xmax) / 2
  )

# =========================
# 图 1：百分比堆积柱状图
p_pct <- ggplot(df_pct, aes(x = Samples, y = Percent, fill = Species)) +
  
  # 背景分组色块
  geom_rect(
    data = bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = bg_df$bg,
    alpha = 0.25
  ) +
  
  # 百分比堆积柱
  geom_bar(stat = "identity", width = 0.8) +
  
  # 物种配色
  scale_fill_manual(values = species_colors, name = "Species") +
  
  # 纵坐标刻度：0 / 25 / 50 / 75 / 100
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 110),
    expand = c(0, 0)
  ) +
  
  # 分组文字标识
  geom_text(
    data = bg_df,
    aes(x = xmid, y = 105, label = Times),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  
  ylab("Relative abundance (%)") +
  xlab("Samples") +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 20, 10)
  )

print(p_pct)


# =========================
# 图 2：原始丰度堆积柱状图
# =========================
p_cnt <- ggplot(df_long, aes(x = Samples, y = Count, fill = Species)) +
  
  # 背景分组色块
  geom_rect(
    data = bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = bg_df$bg,
    alpha = 0.25
  ) +
  
  # 原始丰度堆积柱
  geom_bar(stat = "identity", width = 0.8) +
  
  # 物种配色
  scale_fill_manual(values = species_colors, name = "Species") +
  
  # 分组文字标识
  geom_text(
    data = bg_df,
    aes(x = xmid, y = max(df_long$Count) * 1.05, label = Times),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  
  ylab("Abundance") +
  xlab("Samples") +
  
  # 断轴（18000–15000）
  scale_y_break(c(8000, 11000)) +
  
  coord_cartesian(ylim = c(0, max(df_long$Count) * 1.15)) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 20, 10)
  )

print(p_cnt)

# =========================
# （可选）保存图片
# =========================
ggsave("Treatment_percent_stacked.pdf", p_pct, width = 8, height = 5, dpi = 300)
ggsave("Treatment_count_stacked.pdf",   p_cnt, width = 8, height = 5, dpi = 300)
