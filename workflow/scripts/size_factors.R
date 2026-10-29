# DESeq2 size factor（median-of-ratios）：每个区域取全部样本计数的几何平均作参照，
# 样本的 size factor 为各区域 计数 / 参照 的中位数（任一样本计数为 0 的区域不参与）。
# 假设大多数区域在样本间不变；ChIP 只用 IP 样本。bigwig 用 --scaleFactor 1 / size factor。
# 输入：peak_counts.tsv（DiffBind，前 4 列为 peak、Chr、Start、End）或 multiBamSummary --outRawCounts（列名带引号）
# 输出：sample、size_factor、regions（参与计算的区域数）
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages(library(DESeq2))

samples <- unlist(snakemake@params[["samples"]])
d <- read.delim(snakemake@input[["counts"]], check.names = FALSE, comment.char = "")
names(d) <- gsub("^#|'", "", names(d))
missing <- setdiff(samples, names(d))
if (length(missing)) stop("计数表中没有这些样本：", paste(missing, collapse = ", "))
m <- round(as.matrix(d[samples]))
used <- rowSums(m == 0) == 0
sf <- estimateSizeFactorsForMatrix(m)
out <- data.frame(sample = samples, size_factor = round(sf[samples], 6), regions = sum(used))
write.table(out, snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("%d 个区域，%d 个在全部样本中计数 > 0（参与计算）", nrow(m), sum(used)))
print(out)
