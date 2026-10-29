# 重复一致性的图（QC 报告「重复一致性」一节）。样式见 plot_style.R（dataviz 色板、无网格线、黑色轴线与朝外刻度）
#   宽度 6.9 in（双栏）；PNG 300 dpi + PDF + SVG；图题和说明写在报告的图注中，图内不加标题
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(jsonlite)
})
source(file.path(snakemake@scriptdir, "plot_style.R"))

outdir <- snakemake@params[["outdir"]]
figures <- unlist(snakemake@params[["figures"]])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

info <- fromJSON(snakemake@input[["info"]])
per <- read.delim(snakemake@input[["samples"]], check.names = FALSE)
groups <- names(info$groups)
all_groups <- union(groups, unique(per$group))  # ChIP：control 的组只出现在性别图中
per$group <- factor(per$group, levels = all_groups)

# 分组：颜色与形状同时区分（散点图中超过 3 组时颜色不足以两两区分）
group_col <- cat_colours(all_groups)
group_shape <- cat_shapes(all_groups)
scale_group <- list(
  scale_colour_manual(values = group_col, name = "Group", drop = FALSE),
  scale_shape_manual(values = group_shape, name = "Group", drop = FALSE)
)
fig <- function(name) file.path(outdir, name)
per$label <- short_label(per$sample, per$group)

# ------------------------------------------------------------
# 相关性散点矩阵 + 聚类树
#   左侧为按 1 - r 层次聚类（average）的聚类树与分组色条，样本按聚类顺序排列；
#   下三角为两样本在 consensus peak 上的 log2(CPM + 1)（peak 多时随机抽 info$scatter_peaks 个作图），
#   上三角为相关系数（基于全部 consensus peak），对角为样本名。点用 scattermore 栅格化，文件大小与点数无关
# ------------------------------------------------------------
if ("correlation_scatter" %in% figures) {
  logc <- read.delim(snakemake@input[["log2cpm"]], row.names = 1, check.names = FALSE)
  cor <- as.matrix(read.delim(snakemake@input[["correlation"]], row.names = 1, check.names = FALSE))
  g_of <- setNames(as.character(per$group), per$sample)
  hc <- hclust(as.dist(1 - cor), method = "average")
  ss <- rownames(cor)[hc$order]
  n <- length(ss)

  # 聚类树：叶子 ss[k] 在第 k 行（自上而下），纵坐标 y = n - k + 1；横轴为合并高度（1 - r），根在左
  ypos <- setNames(n - seq_len(n) + 1, ss)
  node_y <- numeric(n - 1)
  seg <- list()
  for (i in seq_len(n - 1)) {
    yh <- lapply(hc$merge[i, ], function(k) if (k < 0) c(ypos[rownames(cor)[-k]], 0) else c(node_y[k], hc$height[k]))
    node_y[i] <- mean(c(yh[[1]][1], yh[[2]][1]))
    h <- hc$height[i]
    seg[[i]] <- data.frame(
      y = c(yh[[1]][1], yh[[2]][1], yh[[1]][1]), yend = c(yh[[1]][1], yh[[2]][1], yh[[2]][1]),
      x = c(yh[[1]][2], yh[[2]][2], h), xend = c(h, h, h)
    )
  }
  seg <- do.call(rbind, seg)
  ys <- scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0))
  void <- theme_void(base_size = 10, base_family = font) +
    theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8))
  p_tree <- ggplot(seg) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.3, colour = "black") +
    scale_x_reverse(expand = expansion(mult = c(0.02, 0))) + ys + void
  strip <- data.frame(y = ypos, group = factor(g_of[ss], levels = groups))
  p_strip <- ggplot(strip, aes(x = 1, y = y, fill = group)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = group_col, name = "Group", drop = FALSE) +
    ys + void

  set.seed(1)
  keep <- sample(nrow(logc), info$scatter_peaks)
  sub <- logc[keep, ss, drop = FALSE]
  lim <- range(as.matrix(sub))
  ij <- expand.grid(i = seq_len(n), j = seq_len(n))
  low <- ij[ij$i > ij$j, ]  # 第 i 行（纵轴）、第 j 列（横轴）
  pts <- do.call(rbind, lapply(seq_len(nrow(low)), function(k) {
    data.frame(row = ss[low$i[k]], col = ss[low$j[k]], x = sub[[ss[low$j[k]]]], y = sub[[ss[low$i[k]]]])
  }))
  up <- ij[ij$i < ij$j, ]
  rt <- data.frame(row = ss[up$i], col = ss[up$j], r = cor[cbind(ss[up$i], ss[up$j])])
  lo <- floor(min(cor) * 100) / 100
  rt$dark <- on_seq_dark(rt$r, lo, 1)
  dg <- data.frame(row = ss, col = ss, label = sub("_(rep[0-9]+)$", "\n\\1", ss, ignore.case = TRUE))  # 末尾的 _repN 换行
  lv <- function(d) transform(d, row = factor(row, levels = ss), col = factor(col, levels = ss))
  pts <- lv(pts); rt <- lv(rt); dg <- lv(dg)
  mid <- mean(lim)
  method <- if (info$method == "spearman") "Spearman" else "Pearson"
  p_mat <- ggplot() +
    geom_rect(data = rt, aes(fill = r), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_text(data = rt, aes(label = sprintf("%.2f", r), colour = dark), x = mid, y = mid, size = 8 / .pt,
              family = font, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    geom_abline(data = unique(pts[c("row", "col")]), aes(slope = 1, intercept = 0), colour = REF_LINE,
                linewidth = 0.2, linetype = "dashed") +  # 只画在下三角
    scattermore::geom_scattermore(data = pts, aes(x, y), pointsize = 1.2, alpha = 0.15, colour = "black",
                                  pixels = c(256, 256)) +
    geom_text(data = dg, aes(label = label), x = mid, y = mid, size = 7 / .pt, family = font, lineheight = 0.9,
              colour = "black") +
    scale_fill_seq(limits = c(lo, 1), name = paste(method, "r")) +
    facet_grid(row ~ col) +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(x = "log2(CPM + 1)", y = "log2(CPM + 1)") +
    theme_qc() +
    theme(strip.text = element_blank(), panel.spacing = unit(0.08, "lines"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
          axis.line = element_blank(), axis.text = element_text(size = 6))
  w <- max(6.9, 0.55 * n + 2.2)  # 每格约 0.55 in；样本多时放大，报告中按宽度缩放
  p <- p_tree + p_strip + p_mat + plot_layout(widths = c(0.9, 0.12, 0.55 * n), guides = "collect")
  save_fig(p, fig("correlation_scatter"), w, w - 1.6)
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
    geom_text_repel(aes(label = label), size = 8 / .pt, family = font, colour = "black", show.legend = FALSE,
                    seed = 1, max.overlaps = Inf, min.segment.length = 0.2, segment.colour = REF_LINE) +
    scale_group +
    labs(x = sprintf("PC1 (%.1f%% of variance)", v[1]), y = sprintf("PC2 (%.1f%% of variance)", v[2])) +
    theme_qc()
  save_fig(p, fig("pca"), 6.9, 4.5)
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
                  colour = "black", show.legend = FALSE) +
    geom_jitter(width = 0.12, height = 0, size = 2.2) +
    scale_group +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Group", y = "Jaccard index of replicate peaks") +
    theme_qc() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p, fig("peak_overlap"), 6.9, 4)
}

# ------------------------------------------------------------
# Venn（组内各重复 / 组间各组）：各组合独有的碱基数（bedtools multiinter），面积不按比例
#   2、3 个集合用圆，4 个集合用椭圆（经典四椭圆布局）；每个区域用 polyclip 求出多边形，
#   按碱基数着色（蓝色顺序色阶），数字放在区域内离边界最远的点（四椭圆用固定位置）
# ------------------------------------------------------------
if (any(c("peak_venn", "group_venn") %in% figures)) {
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

  # d：facet（每张小图）、members（逗号分隔的集合组合）、bp；sets：每张小图的集合（按顺序）；set_label：集合名的显示
  draw_venn <- function(d, sets, set_label, name) {
    facets <- names(sets)
    regions <- list(); labels <- list(); outlines <- list(); set_names <- list()
    for (g in facets) {
      ss <- sets[[g]]
      lay <- layouts[[as.character(length(ss))]]
      sub <- d[d$facet == g, ]
      for (j in seq_len(nrow(sub))) {
        inside <- ss %in% strsplit(sub$members[j], ",")[[1]]
        reg <- lay$shapes[which(inside)[1]]
        for (k in which(inside)[-1]) reg <- polyclip::polyclip(reg, lay$shapes[k], op = "intersection")
        for (k in which(!inside)) reg <- polyclip::polyclip(reg, lay$shapes[k], op = "minus")
        if (!length(reg)) next
        id <- paste(g, sub$members[j])
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
      set_names[[g]] <- cbind(lay$names, group = g, label = set_label(ss, g))
    }
    regions <- do.call(rbind, regions); labels <- do.call(rbind, labels)
    outlines <- do.call(rbind, outlines); set_names <- do.call(rbind, set_names)
    for (v in c("regions", "labels", "outlines", "set_names")) {
      x <- get(v); x$group <- factor(x$group, levels = facets); assign(v, x)
    }
    regions$mb <- regions$bp / 1e6
    labels$mb <- labels$bp / 1e6
    hi_mb <- max(regions$mb)
    labels$text <- ifelse(labels$mb >= 10, sprintf("%.1f", labels$mb), sprintf("%.2f", labels$mb))
    labels$dark <- on_seq_dark(labels$mb, 0, hi_mb)

    ncol <- min(3, length(facets))
    nrow <- ceiling(length(facets) / ncol)
    p <- ggplot() +
      geom_polygon(data = regions, aes(x, y, group = part, fill = mb)) +
      geom_polygon(data = outlines, aes(x, y, group = interaction(group, id)), fill = NA, colour = "black",
                   linewidth = 0.3) +
      geom_text(data = labels, aes(x, y, label = text, colour = dark), size = 8 / .pt, family = font, show.legend = FALSE) +
      geom_text(data = set_names, aes(x, y, label = label), size = 8 / .pt, family = font) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
      scale_fill_seq(name = "Mb", limits = c(0, hi_mb),
                     guide = guide_colourbar(display = "rectangles", nbin = 50, barwidth = unit(2.5, "in"),
                                             barheight = unit(0.12, "in"), title.vjust = 0.9)) +
      facet_wrap(~group, ncol = ncol) +
      coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void(base_size = 10, base_family = font) +
      theme(strip.text = element_text(size = 9, margin = margin(b = 2)),
            legend.title = element_text(size = 9), legend.text = element_text(size = 8),
            panel.spacing = unit(0.4, "lines"), legend.position = "bottom")
    # 色标放在下方，横向留给各组的图；高度按每格接近正方形
    save_fig(p, fig(name), 6.9, 0.7 + 6.9 / ncol * 1.02 * nrow)
  }
}

# peak 重现性：组内各重复
if ("peak_venn" %in% figures) {
  venn <- read.delim(snakemake@input[["venn"]], check.names = FALSE)
  vg <- unlist(snakemake@params[["venn_groups"]])
  draw_venn(data.frame(facet = venn$group, members = venn$samples, bp = venn$bp),
            info$groups[vg], short_label, "peak_venn")
}

# 组间：各组的 reproducible peak（只有 1 个样本的组用该样本的 peak）
if ("group_venn" %in% figures) {
  gv <- read.delim(snakemake@input[["group_venn"]], check.names = FALSE)
  sets <- lapply(snakemake@params[["group_venn"]], unlist)
  draw_venn(data.frame(facet = gv$comparison, members = gv$groups, bp = gv$bp),
            sets, function(ss, g) ss, "group_venn")
}

# ------------------------------------------------------------
# 性别：chrX 与 chrY 计数区 reads 占全部比对 reads 的百分比（chrY 为对数轴），虚线为判定阈值
# ------------------------------------------------------------
if ("sex" %in% figures) {
  th <- info$y_frac
  floor_y <- 1e-5  # 对数轴上 0 显示为 0.00001%
  per$y_plot <- pmax(per$y_frac, floor_y)
  p <- ggplot(per, aes(x_frac, y_plot, colour = group, shape = group)) +
    geom_hline(yintercept = th, linetype = "dashed", linewidth = 0.3, colour = REF_LINE) +
    annotate("text", x = Inf, y = th, label = c("female below", "male above"), hjust = 1.05,
             vjust = c(1.5, -0.6), size = 8 / .pt, family = font, colour = INK2) +
    geom_point(size = 2.5) +
    geom_text_repel(data = per[per$sex != "F", ], aes(label = sample), size = 8 / .pt, family = font,
                    colour = "black", segment.colour = REF_LINE, show.legend = FALSE, seed = 1,
                    min.segment.length = 0, nudge_x = 0.3) +
    scale_group +
    scale_y_log10(labels = function(x) paste0(format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE), "%")) +
    scale_x_continuous(limits = c(0, max(per$x_frac, na.rm = TRUE) * 1.1),
                       labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "chrX reads (% of mapped reads)", y = sprintf("chrY reads in %s (%% of mapped reads)", info$y_region)) +
    theme_qc()
  save_fig(p, fig("sex"), 6.9, 4.5)
}
