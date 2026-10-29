# QC 报告中各图的统一样式，由 replicate_plots.R、fragment_plot.R、profile_plots.R 读入（source）。
# 颜色取自 dataviz skill 的参考色板，已用其 validate_palette.py 校验：
#   分类色：8 个，按固定顺序分配（相邻色对 CVD ΔE ≥ 9.1，正常视觉 ΔE ≥ 19.6）。散点图中任意两色都会相邻，
#           只有前 3 个两两可分，超过 3 组时同时用点形状区分；绿、黄、粉对白底对比度 < 3:1，靠图例和数据表补足
#   顺序色：蓝色单色阶（浅 -> 深），用于相关系数、碱基数等连续数值
# 无网格线；x、y 轴有黑色轴线和朝外的刻度；刻度文字和轴标题为黑色
# 字体为 Arial（缺失时用度量相同的 Liberation Sans）；刻度 8 pt、轴标签 10 pt；图题和说明写在报告的图注中

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"

CAT_COLS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
CAT_SHAPES <- c(16, 15, 17, 18, 3, 4, 8, 1)
SEQ_COLS <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
INK2 <- "#52514e"  # 次要文字（注释）
REF_LINE <- "#898781"  # 参考线（y = x、阈值）

theme_qc <- function() {
  theme_classic(base_size = 10, base_family = font) +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(size = 10, colour = "black"),
      axis.text = element_text(size = 8, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(2.5, "pt"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      strip.background = element_blank(),
      strip.text = element_text(size = 9, hjust = 0),
      panel.grid = element_blank()
    )
}

# 分类按顺序取色；超过 8 个时颜色循环、形状错开，保证颜色与形状的组合不重复
cat_colours <- function(levels) setNames(CAT_COLS[(seq_along(levels) - 1) %% 8 + 1], levels)
cat_shapes <- function(levels) {
  idx <- seq_along(levels) - 1
  setNames(CAT_SHAPES[(idx + idx %/% 8) %% 8 + 1], levels)
}

# 色标用矩形而非位图：位图在 cairo SVG 中靠 mask / feImage 合成，Firefox 不显示
scale_fill_seq <- function(..., guide = guide_colourbar(display = "rectangles", nbin = 50)) {
  scale_fill_gradientn(colours = SEQ_COLS, ..., guide = guide)
}

# 顺序色深端（数值靠近上限）用白字
on_seq_dark <- function(v, lo, hi) (v - lo) / (hi - lo) > 0.55

# 样本名去掉组名前缀作为重复的标签（如 Naive_rep1 -> rep1）
short_label <- function(s, g) {
  g <- as.character(g)
  ifelse(startsWith(s, g) & nchar(s) > nchar(g), sub("^[._-]", "", substring(s, nchar(g) + 1)), s)
}

# 组内重复的颜色：按样本在组内的顺序取分类色。有重复的组的标签一致（如都是 rep1、rep2）时图例用该标签，否则用序号；
# 只有一个样本的组（如 IgG）不参与判断，按其序号取同一标签
replicate_key <- function(samples, groups) {
  idx <- ave(seq_along(samples), groups, FUN = seq_along)
  lab <- short_label(samples, groups)
  multi <- groups %in% groups[duplicated(groups)]
  ref <- if (any(multi)) tapply(lab[multi], idx[multi], unique, simplify = FALSE) else list()
  same <- length(ref) > 0 && all(lengths(ref) == 1) && all(as.character(idx) %in% names(ref))
  key <- if (same) unlist(ref)[as.character(idx)] else as.character(idx)
  levels <- unique(key[order(idx)])
  list(key = factor(key, levels = levels), colours = setNames(CAT_COLS[(seq_along(levels) - 1) %% 8 + 1], levels),
       title = if (same) "Replicate" else "Replicate (order in sample sheet)")
}

save_fig <- function(p, file, w, h) {
  # file 不含扩展名，输出 PNG 300 dpi、PDF、SVG
  ggsave(paste0(file, ".png"), p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(paste0(file, ".pdf"), p, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(paste0(file, ".svg"), p, width = w, height = h, device = svg, bg = "white")
}
