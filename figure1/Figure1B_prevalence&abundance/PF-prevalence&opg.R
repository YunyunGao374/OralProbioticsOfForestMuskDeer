library(tidyverse)
library(patchwork)
library(scales)

col_main <- "#cb3127"

df <- read.table("parasite_all.txt",
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE)

df_treatment <- df %>%
  filter(Groups == "Treatment")
df_treatment$Times <- factor(
  df_treatment$Times,
  levels = c(
    "Hiding_phase",
    "Conlict_phase",
    "Disintegrated_phase",
    "Independent_phase"
  )
)
df_treatment_log <- df_treatment %>%
  mutate(
    total_log = ifelse(total == 0, 0.5, total)
  )
summary_raw <- df_treatment %>%
  group_by(Times) %>%
  summarise(
    mean_raw = mean(total),
    min_raw  = min(total),
    max_raw  = max(total)
  ) %>%
  mutate(
    label = paste0(
      round(mean_raw, 1),
      " (", min_raw, "–", max_raw, ")"
    )
  )
summary_log <- df_treatment_log %>%
  group_by(Times) %>%
  summarise(
    mean_log = mean(total_log),
    sd_log   = sd(total_log)
  )
p1 <- ggplot() +
  geom_jitter(
    data = df_treatment_log,
    aes(x = Times, y = total_log),
    width = 0.15,
    size = 2,
    alpha = 0.7,
    color = col_main
  ) +
  geom_errorbar(
    data = summary_log,
    aes(
      x = Times,
      ymin = mean_log - sd_log,
      ymax = mean_log + sd_log
    ),
    width = 0.25,
    linewidth = 0.9,
    color = col_main
  ) +
  geom_point(
    data = summary_log,
    aes(x = Times, y = mean_log),
    size = 4,
    shape = 21,
    fill = "yellow",
    color = col_main,
    stroke = 1
  ) +
  geom_text(
    data = summary_raw,
    aes(
      x = Times,
      y = summary_log$mean_log[match(Times, summary_log$Times)],
      label = label
    ),
    vjust = -1.3,
    size = 3.5
  ) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    labels = comma
  ) +
  labs(
    y = "Mean faecal oocyst counts (OPG)\nSum of Eimeria spp. (log10)",
    x = NULL
  ) +
  theme_classic(base_size = 14)

rate_df <- df_treatment %>%
  mutate(detected = ifelse(total > 0, 1, 0)) %>%
  group_by(Times) %>%
  summarise(
    detection_rate = mean(detected) * 100
  )

p2 <- ggplot(rate_df, aes(x = Times, y = detection_rate)) +
  geom_col(width = 0.6, fill = col_main) +
  geom_text(
    aes(label = paste0(round(detection_rate, 1), "%")),
    vjust = -0.5,
    size = 4
  ) +
  ylim(0, 105) +
  labs(
    y = "Total oocyst excretion rate (%)",
    x = ""
  ) +
  theme_classic(base_size = 14)

p1 / p2 + plot_layout(heights = c(3, 1))


