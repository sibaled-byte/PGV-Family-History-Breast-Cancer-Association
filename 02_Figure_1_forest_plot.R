library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(patchwork)

# Forest plot of precomputed adjusted estimates from the main analysis.
# This script draws the figure; it does not fit regression models.
# Breast cancer-related PGV forest plot
# Side-by-side: Premenopausal vs Postmenopausal
# Includes:
#   1. Breast cancer-related genes
#   2. BRCA1/2
#   3. HRR

df <- data.frame(
  Group = rep(
    c(
      "Breast cancer-related genes",
      "BRCA1/2",
      "HRR"
    ),
    each = 2
  ),
  
  Menopause = rep(
    c("Premenopausal", "Postmenopausal"),
    times = 3
  ),
  
  # PGV carriers with FH-
  FH_neg_OR = c(
    6.94, 6.14,  # BC-related genes
    7.82, 5.95,  # BRCA1/2
    4.73, 4.58   # HRR
  ),
  
  FH_neg_low = c(
    4.06, 3.27,
    4.17, 2.85,
    2.90, 2.62
  ),
  
  FH_neg_up = c(
    11.88, 11.51,
    14.69, 12.42,
    7.70, 8.02
  ),
  
  # PGV carriers with FH+
  FH_pos_OR = c(
    39.48, 1.36,  # BC-related genes
    36.40, 2.73,  # BRCA1/2
    36.86, 1.58   # HRR
  ),
  
  FH_pos_low = c(
    10.78, 0.29,
    9.80, 0.50,
    9.92, 0.33
  ),
  
  FH_pos_up = c(
    144.62, 6.43,
    135.18, 14.99,
    136.95, 7.63
  ),
  
  # FH x PGV interaction within each menopausal stratum
  p_int = c(
    0.020, 0.031,  # BC-related genes
    0.046, 0.230,  # BRCA1/2
    0.006, 0.097   # HRR
  ),
  
  # FH x PGV x menopausal status interaction
  p_3way = c(
    0.003, 0.003,  # BC-related genes
    0.048, 0.048,  # BRCA1/2
    0.005, 0.005   # HRR
  )
)

df <- df %>%
  mutate(
    Group = factor(
      Group,
      levels = rev(
        c(
          "Breast cancer-related genes",
          "BRCA1/2",
          "HRR"
        )
      )
    ),
    
    # JAMA Network P-value formatting
    p_label = case_when(
      p_int < 0.001 ~ "<.001",
      p_int < 0.01  ~ sub("^0", "", sprintf("%.3f", p_int)),
      p_int <= 0.99 ~ sub("^0", "", sprintf("%.2f", p_int)),
      TRUE          ~ ">.99"
    ),
    
    p_fontface = ifelse(
      p_int < 0.05,
      "bold",
      "plain"
    )
  )

df_long <- df %>%
  pivot_longer(
    cols = c(
      FH_neg_OR,
      FH_neg_low,
      FH_neg_up,
      FH_pos_OR,
      FH_pos_low,
      FH_pos_up
    ),
    names_to = c("FH", ".value"),
    names_pattern = "FH_(.*)_(.*)"
  ) %>%
  mutate(
    FH_label = ifelse(FH == "neg", "FH−", "FH+"),
    
    or_label = sprintf(
      "%.2f (%.2f–%.2f)",
      OR,
      low,
      up
    ),
    
    y_num = as.numeric(Group),
    
    y_plot = ifelse(
      FH_label == "FH+",
      y_num - 0.12,
      y_num + 0.12
    )
  )

plot_colors <- c(
  "FH+" = "#8B1E3F",
  "FH−" = "#1F4E79"
)

make_forest <- function(
    dat,
    dat_long,
    title,
    xmax,
    x_or_text,
    x_pint,
    x_axis_max,
    breaks
) {
  
  dat_long <- dat_long %>%
    mutate(
      low_plot = pmax(low, 1),
      up_plot = pmin(up, xmax),
      arrow_needed = up > xmax,
      left_arrow_needed = low < 1
    )
  
  ggplot() +
    
    geom_segment(
      data = dat_long,
      aes(
        x = low_plot,
        xend = up_plot,
        y = y_plot,
        yend = y_plot,
        color = FH_label
      ),
      linewidth = 0.48
    ) +
    
    geom_segment(
      data = dat_long %>% filter(arrow_needed),
      aes(
        x = xmax * 0.90,
        xend = xmax * 0.98,
        y = y_plot,
        yend = y_plot,
        color = FH_label
      ),
      arrow = arrow(
        length = unit(0.12, "cm"),
        type = "closed"
      ),
      linewidth = 0.48,
      show.legend = FALSE
    ) +
    
    geom_segment(
      data = dat_long %>% filter(left_arrow_needed),
      aes(
        x = 1.02,
        xend = 1.12,
        y = y_plot,
        yend = y_plot,
        color = FH_label
      ),
      arrow = arrow(
        length = unit(0.12, "cm"),
        type = "closed",
        ends = "first"
      ),
      linewidth = 0.48,
      show.legend = FALSE
    ) +
    
    geom_point(
      data = dat_long,
      aes(
        x = OR,
        y = y_plot,
        color = FH_label,
        shape = FH_label
      ),
      size = 2.3
    ) +
    
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "gray60",
      linewidth = 0.4
    ) +
    
    geom_text(
      data = dat_long,
      aes(
        x = x_or_text,
        y = y_plot,
        label = or_label
      ),
      hjust = 0,
      size = 2.25,
      show.legend = FALSE
    ) +
    
    geom_text(
      data = dat,
      aes(
        x = x_pint,
        y = as.numeric(Group),
        label = p_label,
        fontface = "plain"
      ),
      hjust = 0,
      size = 2.4
    ) +
    
    annotate(
      "text",
      x = x_or_text,
      y = length(levels(df$Group)) + 0.40,
      label = "OR (95% CI)",
      hjust = 0,
      size = 2.4,
      fontface = "bold"
    ) +
    
    annotate(
      "text",
      x = x_pint,
      y = length(levels(df$Group)) + 0.40,
      label = expression(bolditalic(P)[int]),
      hjust = 0,
      size = 2.4
    ) +
    
    scale_x_log10(
      breaks = breaks
    ) +
    
    scale_y_continuous(
      breaks = seq_along(levels(df$Group)),
      labels = levels(df$Group),
      expand = expansion(mult = c(0.15, 0.18))
    ) +
    
    scale_color_manual(
      values = plot_colors
    ) +
    
    scale_shape_manual(
      values = c(
        "FH+" = 16,
        "FH−" = 17
      )
    ) +
    
    coord_cartesian(
      xlim = c(1, x_axis_max),
      clip = "off"
    ) +
    
    labs(
      title = title,
      x = "Odds ratio",
      y = NULL
    ) +
    
    theme_bw(base_size = 10) +
    
    theme(
      plot.title = element_text(
        size = 11,
        face = "plain",
        hjust = 0.5
      ),
      axis.text.y = element_text(
        size = 8.5,
        color = "gray20"
      ),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        color = "gray90",
        linewidth = 0.3
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}


# Premenopausal panel


df_pre <- df %>%
  filter(Menopause == "Premenopausal")

df_long_pre <- df_long %>%
  filter(Menopause == "Premenopausal")

p_pre <- make_forest(
  dat = df_pre,
  dat_long = df_long_pre,
  title = "Premenopausal women",
  xmax = 50,
  x_or_text = 65,
  x_pint = 160,
  x_axis_max = 180,
  breaks = c(1, 2, 5, 10, 50)
)


# Postmenopausal panel


df_post <- df %>%
  filter(Menopause == "Postmenopausal")

df_long_post <- df_long %>%
  filter(Menopause == "Postmenopausal")

p_post <- make_forest(
  dat = df_post,
  dat_long = df_long_post,
  title = "Postmenopausal women",
  xmax = 10,
  x_or_text = 12,
  x_pint = 18.7,
  x_axis_max = 20,
  breaks = c(1, 2, 5, 10)
)

# Remove duplicate y-axis labels from postmenopausal panel

p_post_clean <- p_post +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )


# Combine panels


breast_combined <- p_pre + p_post_clean +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(
      face = "bold",
      size = 12
    )
  )

print(breast_combined)


# Save high-resolution TIFF


ggsave(
  filename = "forest_plot_breast_related_genes_BRCA_HRR.tiff",
  plot = breast_combined,
  width = 14.5,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)