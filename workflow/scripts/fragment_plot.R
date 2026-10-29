# PE 样本去重后的片段长度分布（QC 报告「QC 汇总表」一节）。样式见 plot_style.R；宽 6.9 in；PNG 300 dpi + PDF + SVG
# 每组一个分面（共用坐标轴），组内每个重复一条线，颜色按重复（图例说明）；对照用虚线
# 120 bp 处的竖线为亚核小体片段的上限（summary.tsv 的 subnucleosomal_rate）
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages(library(ggplot2))
source(file.path(snakemake@scriptdir, "plot_style.R"))

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

# 重复的颜色：样本按样本表顺序
ss <- names(groups)[names(groups) %in% unique(d$sample)]
rk <- replicate_key(ss, groups[ss])
d$rep <- factor(as.character(rk$key)[match(d$sample, ss)], levels = levels(rk$key))

ncol <- min(4, length(all_groups))
nrow <- ceiling(length(all_groups) / ncol)
p <- ggplot(d, aes(length, pct, group = sample, colour = rep, linetype = role)) +
  annotate("rect", xmin = -Inf, xmax = SUBNUC, ymin = -Inf, ymax = Inf, fill = "#f0efec") +
  geom_vline(xintercept = SUBNUC, colour = REF_LINE, linewidth = 0.3, linetype = "22") +
  geom_line(linewidth = 0.4) +
  facet_wrap(~group, ncol = ncol) +
  scale_colour_manual(values = rk$colours, name = rk$title) +
  scale_linetype_manual(values = c(IP = "solid", control = "22"), name = NULL,
                        guide = if (any(d$role == "control")) "legend" else "none") +
  scale_x_continuous(limits = c(0, max_len), breaks = seq(0, max_len, 100), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Fragment length (bp)", y = "Fragments (%)") +
  theme_qc() +
  theme(legend.position = "top", legend.justification = "left", panel.spacing.x = unit(0.8, "lines"),
        plot.margin = margin(5.5, 12, 5.5, 5.5))

out <- unlist(snakemake@output)
w <- 6.9
h <- 0.9 + 1.9 * nrow
for (f in out) {
  dev <- switch(tools::file_ext(f), png = ragg::agg_png, pdf = cairo_pdf, svg = svg)
  ggsave(f, p, width = w, height = h, dpi = 300, device = dev, bg = "white")
}
