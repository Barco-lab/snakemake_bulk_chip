# 信号 profile 的图（QC 报告「信号 profile」一节）。样式见 plot_style.R；宽 6.9 in；PNG 300 dpi + PDF + SVG
#   metagene_samples    meta gene profile（plotProfile --outFileNameData 的平均信号），每组一个分面，组内每个重复一条线
#   metagene_groups     组 bigwig 的 meta gene profile，每格最多 3 组
#   genome_track        全基因组 bin 信号沿染色体分布：每个样本 log2(bin 平均信号 / 该样本全部 bin 的中位数)，
#                       multiBigwigSummary 的小 bin 先合并为 track_bin 再作图（信号为 0 的 bin，如组装缺口，不画）；每组一个分面
#   genome_distribution 全基因组 bin 信号的分布：信号 > 0 的 bin 的 log2(信号) 的密度曲线
# 只用主要染色体（1–22 / 1–19、X、Y，命名带不带 chr 均可）
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})
source(file.path(snakemake@scriptdir, "plot_style.R"))

figures <- unlist(snakemake@params[["figures"]])
outdir <- snakemake@params[["outdir"]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
groups <- unlist(snakemake@params[["groups"]])  # 样本 -> 组（样本表顺序）
controls <- unlist(snakemake@params[["controls"]])
norm <- snakemake@params[["norm"]]
fig <- function(name) file.path(outdir, name)

# 组的顺序：IP 样本的组在前，对照的组在后（各自按样本表顺序）
group_order <- function(ss) {
  g <- groups[ss]
  unique(c(g[!ss %in% controls], g[ss %in% controls]))
}

# 每组一个分面、组内各重复一条线（颜色按重复，图例说明；对照为虚线）
add_replicates <- function(d) {
  ss <- names(groups)[names(groups) %in% unique(d$sample)]
  rk <- replicate_key(ss, groups[ss])
  d$group <- factor(groups[d$sample], levels = group_order(ss))
  d$rep <- factor(as.character(rk$key)[match(d$sample, ss)], levels = levels(rk$key))
  d$role <- ifelse(d$sample %in% controls, "control", "IP")
  list(d = d, scales = list(
    scale_colour_manual(values = rk$colours, name = rk$title),
    scale_linetype_manual(values = c(IP = "solid", control = "22"), name = NULL,
                          guide = if (any(d$role == "control")) "legend" else "none")
  ))
}

top_legend <- theme(legend.position = "top", legend.justification = "left", panel.spacing.x = unit(1.6, "lines"),
                    plot.margin = margin(5.5, 14, 5.5, 5.5))

kb <- function(x) ifelse(x %% 1000 == 0, sprintf("%d kb", x %/% 1000), sprintf("%.1f kb", x / 1000))

# ------------------------------------------------------------
# meta gene profile
# ------------------------------------------------------------
lay <- snakemake@params[["layout"]]  # upstream、unscaled5、body、unscaled3、downstream、binsize（bp）

# plotProfile 的数据表：前两行为 bin 标签和序号，之后每行一个样本（第 1 列名称、第 2 列区域组、其后为各 bin 的平均信号）
read_profile <- function(path) {
  lines <- readLines(path)[-(1:2)]
  rows <- lapply(strsplit(lines, "\t"), function(f) {
    v <- as.numeric(f[-(1:2)])
    v <- v[!is.na(v)]
    data.frame(sample = f[1], bin = seq_along(v), signal = v)
  })
  d <- do.call(rbind, rows)
  body_len <- lay$unscaled5 + lay$body + lay$unscaled3
  nb <- max(d$bin)
  if (nb * lay$binsize != lay$upstream + body_len + lay$downstream) {
    stop(sprintf("bin 数（%d）与 computeMatrix 参数不符", nb))
  }
  d$x <- (d$bin - 0.5) * lay$binsize - lay$upstream  # TSS 为 0，TES 为 body_len
  d
}

metagene_x <- function() {
  body_len <- lay$unscaled5 + lay$body + lay$unscaled3
  scale_x_continuous(
    limits = c(-lay$upstream, body_len + lay$downstream),
    breaks = c(-lay$upstream, 0, body_len, body_len + lay$downstream),
    labels = c(paste0("-", kb(lay$upstream)), "TSS", "TES", paste0("+", kb(lay$downstream))),
    expand = c(0, 0)
  )
}
metagene_ylab <- sprintf("Mean signal (%s)", norm)

if ("metagene_samples" %in% figures) {
  r <- add_replicates(read_profile(snakemake@input[["metagene_samples"]]))
  ng <- nlevels(r$d$group)
  ncol <- min(3, ng)
  p <- ggplot(r$d, aes(x, signal, group = sample, colour = rep, linetype = role)) +
    geom_line(linewidth = 0.45) +
    facet_wrap(~group, ncol = ncol) +
    r$scales + metagene_x() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    labs(x = "Protein-coding genes (gene body scaled)", y = metagene_ylab) +
    theme_qc() + top_legend
  save_fig(p, fig("metagene_samples"), 6.9, 0.9 + 2.1 * ceiling(ng / ncol))
}

if ("metagene_groups" %in% figures) {
  d <- read_profile(snakemake@input[["metagene_groups"]])
  gs <- unique(d$sample)
  d$group <- factor(d$sample, levels = gs)
  # 每格最多 3 组：按组的顺序每 3 个一格，颜色跟随组（各格共用一个图例）
  d$panel <- factor(ceiling(match(d$sample, gs) / 3))
  np <- nlevels(d$panel)
  ncol <- min(3, np)
  p <- ggplot(d, aes(x, signal, colour = group)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~panel, ncol = ncol) +
    scale_colour_manual(values = cat_colours(gs), name = "Group") +
    metagene_x() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    labs(x = "Protein-coding genes (gene body scaled)", y = metagene_ylab) +
    theme_qc() + top_legend + theme(strip.text = element_blank())
  w <- if (np == 1) 4.6 else 6.9
  save_fig(p, fig("metagene_groups"), w, 0.9 + 2.3 * ceiling(np / ncol))
}

# ------------------------------------------------------------
# 全基因组 bin
# ------------------------------------------------------------
if (any(c("genome_track", "genome_distribution") %in% figures)) {
  bins <- fread(snakemake@input[["bins"]])
  setnames(bins, gsub("^#|'", "", names(bins)))
  setnames(bins, 1:3, c("chr", "start", "end"))
  main <- grepl("^(chr)?([0-9]+|X|Y)$", bins$chr)
  bins <- bins[main]
  ss <- setdiff(names(bins), c("chr", "start", "end"))
  chr_key <- sub("^chr", "", unique(bins$chr))
  chr_rank <- suppressWarnings(as.numeric(chr_key))
  chr_rank[chr_key == "X"] <- 1000
  chr_rank[chr_key == "Y"] <- 1001
  chroms <- unique(bins$chr)[order(chr_rank)]
}

if ("genome_track" %in% figures) {
  tb <- snakemake@params[["track_bin"]]
  long <- melt(bins, id.vars = c("chr", "start", "end"), variable.name = "sample", value.name = "signal",
               variable.factor = FALSE)
  long <- long[!is.na(signal)]
  long[, tbin := start %/% tb]
  tr <- long[, .(signal = mean(signal)), by = .(sample, chr, tbin)]
  tr <- tr[signal > 0]  # 无覆盖的 bin（组装缺口等）不画，也不计入中位数
  tr[, med := median(signal), by = sample]
  tr[, ratio := log2(signal / med)]
  chr_len <- bins[, .(len = max(end)), by = chr][match(chroms, chr)]
  chr_len[, offset := cumsum(as.numeric(len)) - len]
  tr[, x := (tbin + 0.5) * tb + chr_len$offset[match(chr, chr_len$chr)]]
  tr <- tr[order(sample, x)]
  # 相邻染色体之间断开（每条染色体单独成线）
  tr[, line := paste(sample, chr)]
  r <- add_replicates(as.data.frame(tr))
  shade <- chr_len[seq(2, .N, by = 2), .(xmin = offset, xmax = offset + len)]
  lim <- 3
  ng <- nlevels(r$d$group)
  p <- ggplot(r$d, aes(x, pmax(pmin(ratio, lim), -lim), group = line, colour = rep, linetype = role)) +
    geom_rect(data = shade, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf, fill = "#f0efec",
              inherit.aes = FALSE) +
    geom_hline(yintercept = 0, colour = REF_LINE, linewidth = 0.3) +
    geom_line(linewidth = 0.25) +
    facet_wrap(~group, ncol = 1) +
    r$scales +
    scale_x_continuous(breaks = chr_len$offset + chr_len$len / 2, labels = sub("^chr", "", chr_len$chr),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(-lim, lim), breaks = seq(-lim, lim, 1)) +
    labs(x = "Chromosome", y = "log2(signal / sample median)") +
    theme_qc() + top_legend + theme(axis.text.x = element_text(size = 6.5))
  save_fig(p, fig("genome_track"), 6.9, 0.9 + 1.5 * ng)
}

if ("genome_distribution" %in% figures) {
  # 每个样本：信号 > 0 的 bin 的 log2(信号) 的密度（density，默认带宽），共用横轴范围（全部样本的 0.1%–99.9% 分位）
  lv <- lapply(ss, function(s) {
    v <- bins[[s]]
    log2(v[!is.na(v) & v > 0])
  })
  xr <- quantile(unlist(lv), c(0.001, 0.999))
  dens <- rbindlist(lapply(seq_along(ss), function(k) {
    d <- density(lv[[k]], from = xr[1], to = xr[2], n = 512)
    data.table(sample = ss[k], x = d$x, y = d$y)
  }))
  r <- add_replicates(as.data.frame(dens))
  ng <- nlevels(r$d$group)
  ncol <- min(3, ng)
  p <- ggplot(r$d, aes(x, y, group = sample, colour = rep, linetype = role)) +
    geom_line(linewidth = 0.45) +
    facet_wrap(~group, ncol = ncol) +
    r$scales +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    labs(x = sprintf("log2(bin signal, %s)", norm), y = "Density") +
    theme_qc() + top_legend
  save_fig(p, fig("genome_distribution"), if (ncol == 1) 3.6 else 6.9, 0.9 + 2.1 * ceiling(ng / ncol))
}
