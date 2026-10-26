# 重复一致性的图（QC 报告「重复一致性」一节）。
# 按 academic-research-skills 的 statistical_visualization_standards（APA 7.0）：
#   sans-serif（Arial，缺失时用度量相同的 Liberation Sans），刻度 8 pt、轴标签 10 pt；去掉上、右边框
#   分组用 Tol 色盲友好色板，颜色与点形状同时区分；宽度 6.9 in（双栏）；PNG 300 dpi + PDF + SVG
#   图题和说明写在报告的图注中，图内不加标题
# 相关性热图偏离规范的发散色板：样本间 r 均为正，改用 viridis 顺序色板（范围为最小 r 到 1）
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(jsonlite)
})

outdir <- snakemake@params[["outdir"]]
figures <- unlist(snakemake@params[["figures"]])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

info <- fromJSON(snakemake@input[["info"]])
per <- read.delim(snakemake@input[["samples"]], check.names = FALSE)
groups <- names(info$groups)
all_groups <- union(groups, unique(per$group))  # ChIP：control 的组只出现在性别图中
per$group <- factor(per$group, levels = all_groups)

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
theme_apa <- theme_minimal(base_size = 10, base_family = font) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(colour = "black", linewidth = 0.3)
  )

# 分组：颜色与形状组合，超过 8 组时错开形状，保证组合不重复
cb_palette <- c("#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311", "#EE3377", "#BBBBBB", "#000000")
shapes <- c(16, 15, 17, 18, 3, 4, 8, 1)
idx <- seq_along(all_groups) - 1
group_col <- setNames(cb_palette[idx %% 8 + 1], all_groups)
group_shape <- setNames(shapes[(idx + idx %/% 8) %% 8 + 1], all_groups)
scale_group <- list(
  scale_colour_manual(values = group_col, name = "Group", drop = FALSE),
  scale_shape_manual(values = group_shape, name = "Group", drop = FALSE)
)

save_fig <- function(p, name, w, h) {
  f <- file.path(outdir, name)
  ggsave(paste0(f, ".png"), p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(paste0(f, ".pdf"), p, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(paste0(f, ".svg"), p, width = w, height = h, device = svg, bg = "white")
}

# 样本名去掉组名前缀作为点的标签（如 Naive_rep1 -> rep1）
short_label <- function(s, g) {
  g <- as.character(g)
  ifelse(startsWith(s, g) & nchar(s) > nchar(g), sub("^[._-]", "", substring(s, nchar(g) + 1)), s)
}
per$label <- short_label(per$sample, per$group)

# ------------------------------------------------------------
# 相关性热图：按 1 - r 层次聚类（average），上方为聚类树与分组色条
# ------------------------------------------------------------
if ("correlation" %in% figures) {
  cor <- as.matrix(read.delim(snakemake@input[["correlation"]], row.names = 1, check.names = FALSE))
  n <- nrow(cor)
  hc <- hclust(as.dist(1 - cor), method = "average")
  ord <- rownames(cor)[hc$order]

  # 聚类树的线段：叶子的横坐标为其在 ord 中的位置
  pos <- setNames(seq_len(n), ord)
  node_x <- numeric(n - 1)
  node_y <- hc$height
  seg <- list()
  for (i in seq_len(n - 1)) {
    xy <- lapply(hc$merge[i, ], function(k) {
      if (k < 0) c(pos[rownames(cor)[-k]], 0) else c(node_x[k], node_y[k])
    })
    node_x[i] <- mean(c(xy[[1]][1], xy[[2]][1]))
    seg[[i]] <- data.frame(
      x = c(xy[[1]][1], xy[[2]][1], xy[[1]][1]),
      xend = c(xy[[1]][1], xy[[2]][1], xy[[2]][1]),
      y = c(xy[[1]][2], xy[[2]][2], node_y[i]),
      yend = c(node_y[i], node_y[i], node_y[i])
    )
  }
  seg <- do.call(rbind, seg)
  xs <- scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0))
  void <- theme_void(base_size = 10, base_family = font) +
    theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8))

  p_tree <- ggplot(seg) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.3) +
    xs + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + void

  strip <- data.frame(x = seq_len(n), group = factor(per$group[match(ord, per$sample)], levels = groups))
  p_strip <- ggplot(strip, aes(x = x, y = 1, fill = group)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = group_col, name = "Group", drop = FALSE) +
    xs + void

  long <- data.frame(
    x = rep(seq_len(n), each = n), y = rep(seq_len(n), times = n),
    r = as.vector(cor[ord, ord])
  )
  lo <- floor(min(cor) * 100) / 100
  method <- if (info$method == "spearman") "Spearman" else "Pearson"
  p_heat <- ggplot(long, aes(x = x, y = y, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    # 色标用矩形而非位图：位图在 cairo SVG 中靠 mask / feImage 合成，Firefox 不显示
    scale_fill_viridis_c(limits = c(lo, 1), name = paste(method, "r"), option = "viridis",
                         guide = guide_colourbar(display = "rectangles", nbin = 50)) +
    scale_x_continuous(breaks = seq_len(n), labels = ord, expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_len(n), labels = ord, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5)) +
    labs(x = NULL, y = NULL) +
    theme_apa +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()
    )
  if (n <= 15) {  # 格子太小时 8 pt 的数字放不下
    mid <- (lo + 1) / 2
    p_heat <- p_heat + geom_text(aes(label = sprintf("%.2f", r), colour = r < mid), size = 8 / .pt,
                                 family = font, show.legend = FALSE) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"))
  }
  p <- p_tree / p_strip / p_heat + plot_layout(heights = c(0.6, 0.12, 4), guides = "collect")
  save_fig(p, "correlation", 6.9, 6.6)
}

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------
if ("pca" %in% figures) {
  pca <- read.delim(snakemake@input[["pca"]], check.names = FALSE)
  pca$group <- factor(pca$group, levels = groups)
  pca$label <- short_label(pca$sample, pca$group)
  v <- info$pca_variance * 100
  p <- ggplot(pca, aes(PC1, PC2, colour = group, shape = group)) +
    geom_point(size = 2.5) +
    geom_text_repel(aes(label = label), size = 8 / .pt, family = font, show.legend = FALSE,
                    seed = 1, max.overlaps = Inf, min.segment.length = 0.2) +
    scale_group +
    labs(x = sprintf("PC1 (%.1f%% of variance)", v[1]), y = sprintf("PC2 (%.1f%% of variance)", v[2])) +
    theme_apa
  save_fig(p, "pca", 6.9, 4.5)
}

# ------------------------------------------------------------
# peak 重现性：组内每对重复的 Jaccard，横线为中位数
# ------------------------------------------------------------
if ("peak_overlap" %in% figures) {
  ov <- read.delim(snakemake@input[["overlap"]], check.names = FALSE)
  ov$group <- factor(ov$group, levels = groups)
  ov <- ov[!is.na(ov$jaccard), ]
  med <- aggregate(jaccard ~ group, ov, median)
  set.seed(1)
  p <- ggplot(ov, aes(group, jaccard, colour = group, shape = group)) +
    geom_crossbar(data = med, aes(ymin = jaccard, ymax = jaccard), width = 0.5, linewidth = 0.3,
                  colour = "grey30", show.legend = FALSE) +
    geom_jitter(width = 0.12, height = 0, size = 2.2) +
    scale_group +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Group", y = "Jaccard index of replicate peaks") +
    theme_apa +
    theme(legend.position = "none", panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p, "peak_overlap", 6.9, 4)
}

# ------------------------------------------------------------
# peak 重现性 Venn：组内各重复组合独有的碱基数（bedtools multiinter），面积不按比例
#   2、3 个重复用圆，4 个重复用椭圆（经典四椭圆布局）；每个区域用 polyclip 求出多边形，
#   按碱基数着色（viridis），数字放在区域内离边界最远的点（四椭圆用固定位置）
# ------------------------------------------------------------
if ("peak_venn" %in% figures) {
  ellipse <- function(x, y, w, h, angle = 0, k = 360) {
    t <- seq(0, 2 * pi, length.out = k + 1)[-1]
    a <- angle * pi / 180
    ex <- w / 2 * cos(t); ey <- h / 2 * sin(t)
    list(x = x + ex * cos(a) - ey * sin(a), y = y + ex * sin(a) + ey * cos(a))
  }
  layouts <- list(
    `2` = list(shapes = list(ellipse(0.36, 0.5, 0.5, 0.5), ellipse(0.64, 0.5, 0.5, 0.5)),
               names = data.frame(x = c(0.22, 0.78), y = c(0.82, 0.82))),
    `3` = list(shapes = list(ellipse(0.38, 0.6, 0.44, 0.44), ellipse(0.62, 0.6, 0.44, 0.44),
                             ellipse(0.5, 0.39, 0.44, 0.44)),
               names = data.frame(x = c(0.15, 0.85, 0.5), y = c(0.86, 0.86, 0.1))),
    `4` = list(shapes = list(ellipse(0.35, 0.4, 0.72, 0.45, 140), ellipse(0.45, 0.5, 0.72, 0.45, 140),
                             ellipse(0.544, 0.5, 0.72, 0.45, 40), ellipse(0.644, 0.4, 0.72, 0.45, 40)),
               names = data.frame(x = c(0.12, 0.17, 0.83, 0.88), y = c(0.16, 0.85, 0.85, 0.16)),
               # 四椭圆的窄区域自动定位会挤在一起，用手工位置（取自 Intervene / pyvenn 的 venn4），键为各重复是否在内
               pos = list(`0001` = c(0.85, 0.42), `0010` = c(0.68, 0.72), `0011` = c(0.77, 0.59), `0100` = c(0.32, 0.72),
                          `0101` = c(0.71, 0.30), `0110` = c(0.50, 0.66), `0111` = c(0.65, 0.50), `1000` = c(0.14, 0.42),
                          `1001` = c(0.50, 0.17), `1010` = c(0.29, 0.30), `1011` = c(0.39, 0.24), `1100` = c(0.23, 0.59),
                          `1101` = c(0.61, 0.24), `1110` = c(0.35, 0.50), `1111` = c(0.50, 0.38)))
  )
  area <- function(r) abs(sum(r$x * c(r$y[-1], r$y[1]) - c(r$x[-1], r$x[1]) * r$y)) / 2
  # 区域内离边界最远的格点（近似 pole of inaccessibility）
  label_point <- function(poly) {
    gx <- seq(min(poly$x), max(poly$x), length.out = 60)
    gy <- seq(min(poly$y), max(poly$y), length.out = 60)
    g <- expand.grid(x = gx, y = gy)
    g <- g[sp::point.in.polygon(g$x, g$y, poly$x, poly$y) == 1, ]
    d <- vapply(seq_len(nrow(g)), function(i) min((poly$x - g$x[i])^2 + (poly$y - g$y[i])^2), 0)
    g[which.max(d), ]
  }

  venn <- read.delim(snakemake@input[["venn"]], check.names = FALSE)
  vg <- unlist(snakemake@params[["venn_groups"]])
  regions <- list(); labels <- list(); outlines <- list(); set_names <- list()
  for (g in vg) {
    ss <- info$groups[[g]]
    lay <- layouts[[as.character(length(ss))]]
    sub <- venn[venn$group == g, ]
    for (j in seq_len(nrow(sub))) {
      inside <- ss %in% strsplit(sub$samples[j], ",")[[1]]
      reg <- lay$shapes[which(inside)[1]]
      for (k in which(inside)[-1]) reg <- polyclip::polyclip(reg, lay$shapes[k], op = "intersection")
      for (k in which(!inside)) reg <- polyclip::polyclip(reg, lay$shapes[k], op = "minus")
      if (!length(reg)) next
      id <- paste(g, sub$samples[j])
      for (ci in seq_along(reg)) {
        regions[[length(regions) + 1]] <- data.frame(group = g, id = id, part = paste(id, ci),
                                                     x = reg[[ci]]$x, y = reg[[ci]]$y, bp = sub$bp[j])
      }
      key <- paste(as.integer(inside), collapse = "")
      xy <- if (!is.null(lay$pos)) {
        data.frame(x = lay$pos[[key]][1], y = lay$pos[[key]][2])
      } else {
        label_point(reg[[which.max(vapply(reg, area, 0))]])
      }
      labels[[length(labels) + 1]] <- cbind(xy, group = g, bp = sub$bp[j])
    }
    for (k in seq_along(ss)) {
      outlines[[length(outlines) + 1]] <- data.frame(group = g, id = k, x = lay$shapes[[k]]$x, y = lay$shapes[[k]]$y)
    }
    set_names[[g]] <- cbind(lay$names, group = g, label = short_label(ss, g))
  }
  regions <- do.call(rbind, regions); labels <- do.call(rbind, labels)
  outlines <- do.call(rbind, outlines); set_names <- do.call(rbind, set_names)
  for (d in c("regions", "labels", "outlines", "set_names")) {
    x <- get(d); x$group <- factor(x$group, levels = vg); assign(d, x)
  }
  regions$mb <- regions$bp / 1e6
  labels$mb <- labels$bp / 1e6
  hi_mb <- max(regions$mb)
  labels$text <- ifelse(labels$mb >= 10, sprintf("%.1f", labels$mb), sprintf("%.2f", labels$mb))
  labels$dark <- labels$mb < hi_mb / 2  # viridis 低端颜色深，用白字

  ncol <- min(3, length(vg))
  nrow <- ceiling(length(vg) / ncol)
  p <- ggplot() +
    geom_polygon(data = regions, aes(x, y, group = part, fill = mb)) +
    geom_polygon(data = outlines, aes(x, y, group = interaction(group, id)), fill = NA, colour = "black",
                 linewidth = 0.3) +
    geom_text(data = labels, aes(x, y, label = text, colour = dark), size = 8 / .pt, family = font, show.legend = FALSE) +
    geom_text(data = set_names, aes(x, y, label = label), size = 8 / .pt, family = font) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_viridis_c(name = "Mb", limits = c(0, hi_mb),
                         guide = guide_colourbar(display = "rectangles", nbin = 50, barwidth = unit(2.5, "in"),
                                                 barheight = unit(0.12, "in"), title.vjust = 0.9)) +
    facet_wrap(~group, ncol = ncol) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_size = 10, base_family = font) +
    theme(strip.text = element_text(size = 9, margin = margin(b = 2)),
          legend.title = element_text(size = 9), legend.text = element_text(size = 8),
          panel.spacing = unit(0.4, "lines"), legend.position = "bottom")
  # 色标放在下方，横向留给各组的图；高度按每格接近正方形
  save_fig(p, "peak_venn", 6.9, 0.7 + 6.9 / ncol * 1.02 * nrow)
}

# ------------------------------------------------------------
# 性别：chrX 与 chrY 计数区 reads 占全部比对 reads 的百分比（chrY 为对数轴），虚线为判定阈值
# ------------------------------------------------------------
if ("sex" %in% figures) {
  th <- info$y_frac
  floor_y <- 1e-5  # 对数轴上 0 显示为 0.00001%
  per$y_plot <- pmax(per$y_frac, floor_y)
  p <- ggplot(per, aes(x_frac, y_plot, colour = group, shape = group)) +
    geom_hline(yintercept = th, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    annotate("text", x = Inf, y = th, label = c("female below", "male above"), hjust = 1.05,
             vjust = c(1.5, -0.6), size = 8 / .pt, family = font, colour = "grey30") +
    geom_point(size = 2.5) +
    geom_text_repel(data = per[per$sex != "F", ], aes(label = sample), size = 8 / .pt, family = font,
                    show.legend = FALSE, seed = 1, min.segment.length = 0, nudge_x = 0.3) +
    scale_group +
    scale_y_log10(labels = function(x) paste0(format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE), "%")) +
    scale_x_continuous(limits = c(0, max(per$x_frac, na.rm = TRUE) * 1.1),
                       labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "chrX reads (% of mapped reads)", y = sprintf("chrY reads in %s (%% of mapped reads)", info$y_region)) +
    theme_apa
  save_fig(p, "sex", 6.9, 4.5)
}
