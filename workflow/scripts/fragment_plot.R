# PE 样本去重后的片段长度分布（QC 报告「QC 汇总表」一节）。
# 风格同 replicate_plots.R：Arial / Liberation Sans，刻度 8 pt、轴标签 10 pt；Tol 色盲友好色板；宽 6.9 in；PNG 300 dpi + PDF + SVG
# 每组一个分面（共用坐标轴），组内每个样本一条线，在峰顶标注样本名（去掉组名前缀）；对照用虚线
# 120 bp 处的竖线为亚核小体片段的上限（summary.tsv 的 subnucleosomal_rate）
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

d <- read.delim(snakemake@input[[1]])
groups <- unlist(snakemake@params[["groups"]])
controls <- unlist(snakemake@params[["controls"]])
max_len <- snakemake@params[["max_length"]]
SUBNUC <- 120

d$group <- groups[d$sample]
d$role <- ifelse(d$sample %in% controls, "control", "IP")
# 组的顺序：IP 样本的组在前（按样本表顺序），对照的组在后
all_groups <- unique(c(groups[!names(groups) %in% controls], groups[names(groups) %in% controls]))
all_groups <- all_groups[all_groups %in% d$group]
d$group <- factor(d$group, levels = all_groups)
d <- d[d$length <= max_len, ]

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
theme_apa <- theme_minimal(base_size = 10, base_family = font) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 9, hjust = 0),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(colour = "black", linewidth = 0.3)
  )
cb_palette <- c("#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311", "#EE3377", "#BBBBBB", "#000000")
group_col <- setNames(cb_palette[(seq_along(all_groups) - 1) %% 8 + 1], all_groups)

# 每个样本在峰顶处标注（去掉组名前缀，如 K27AC_A_rep1 -> rep1）
peak <- do.call(rbind, lapply(split(d, d$sample), function(x) x[which.max(x$pct), ]))
peak$label <- mapply(function(s, g) {
  short <- sub(paste0("^", g, "_?"), "", s)
  if (nchar(short)) short else s
}, peak$sample, as.character(peak$group))

ncol <- min(4, length(all_groups))
nrow <- ceiling(length(all_groups) / ncol)
p <- ggplot(d, aes(length, pct, group = sample, colour = group, linetype = role)) +
  annotate("rect", xmin = -Inf, xmax = SUBNUC, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = SUBNUC, colour = "grey55", linewidth = 0.3, linetype = "22") +
  geom_line(linewidth = 0.4) +
  geom_text_repel(data = peak, aes(label = label), size = 7 / .pt, family = font, colour = "grey25",
                  show.legend = FALSE, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey60",
                  box.padding = 0.2, max.overlaps = Inf, seed = 1) +
  facet_wrap(~group, ncol = ncol) +
  scale_colour_manual(values = group_col, guide = "none") +
  scale_linetype_manual(values = c(IP = "solid", control = "22"), name = NULL,
                        guide = if (length(controls) && any(d$role == "control")) "legend" else "none") +
  scale_x_continuous(limits = c(0, max_len), breaks = seq(0, max_len, 100), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Fragment length (bp)", y = "Fragments (%)") +
  theme_apa +
  theme(legend.position = "top", legend.justification = "left", panel.spacing.x = unit(0.8, "lines"),
        plot.margin = margin(5.5, 12, 5.5, 5.5))

out <- unlist(snakemake@output)
w <- 6.9
h <- 0.9 + 1.9 * nrow
for (f in out) {
  dev <- switch(tools::file_ext(f), png = ragg::agg_png, pdf = cairo_pdf, svg = svg)
  ggsave(f, p, width = w, height = h, dpi = 300, device = dev, bg = "white")
}
