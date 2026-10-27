# consensus peak 与计数（DiffBind）：
#   1. 每组：dba.peakset(peaks = 组内样本, minOverlap = min(min_reps, 组内样本数)) 得到组内 consensus
#      （各重复的 peak ≥ 1 bp 重叠即合并，保留至少出现在 minOverlap 个重复中的区域）-> reproducible/{group}.bed
#      只有 1 个样本的组直接用该样本的 peak
#   2. 各组 consensus 取并集（GenomicRanges::reduce）-> consensus peak
#   3. dba.count 在 consensus peak 上计数：summits > 0 时按全部样本的 reads 堆积找峰顶并取 ± summits bp（固定宽度），
#      summits 为 FALSE 时保留原宽度（broad peak）；不扣除对照、不过滤低计数区域；SE 的 reads 按 fragment_size 延长
# 输出：reproducible/{group}.bed、consensus.bed、peak_counts.tsv（reads 数矩阵）、diffbind.rds（DBA 对象，可接着做差异分析）
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressPackageStartupMessages({
  library(DiffBind)
  library(GenomicRanges)
})

samples <- unlist(snakemake@params[["samples"]])
groups <- snakemake@params[["groups"]]
rep_groups <- unlist(snakemake@params[["rep_groups"]])
min_reps <- snakemake@params[["min_reps"]]
summits <- snakemake@params[["summits"]]
group_of <- setNames(rep(names(groups), lengths(groups)), unlist(groups))

ss <- data.frame(
  SampleID = samples,
  Condition = group_of[samples],
  Replicate = ave(seq_along(samples), group_of[samples], FUN = seq_along),
  bamReads = unlist(snakemake@input[["bams"]]),
  Peaks = unlist(snakemake@input[["peaks"]]),
  # narrowPeak / broadPeak 都按 bed 读入，第 7 列 signalValue 作为分数（只影响排序，不影响区域和计数）
  PeakCaller = "bed",
  PeakFormat = "bed",
  ScoreCol = 7,
  LowerBetter = FALSE
)
print(ss)

# blacklist 已在 peak 中去除；greylist 需要对照，这里不用
obj <- dba(sampleSheet = ss, config = data.frame(
  RunParallel = TRUE, fragmentSize = snakemake@params[["fragment_size"]],
  doBlacklist = FALSE, doGreylist = FALSE
))
obj$config$cores <- snakemake@threads

to_bed <- function(gr, path, score = NULL) {
  d <- data.frame(seqnames(gr), start(gr) - 1L, end(gr))
  if (!is.null(score)) d$score <- score
  write.table(d, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
read_peaks <- function(path) {
  d <- read.delim(path, header = FALSE)[, 1:3]
  GRanges(d[[1]], IRanges(d[[2]] + 1L, d[[3]]))
}
peaks_of <- lapply(setNames(ss$Peaks, samples), read_peaks)

# 1. 组内 consensus
group_sets <- list()
rep_paths <- setNames(unlist(snakemake@output[["reproducible"]]), basename(unlist(snakemake@output[["reproducible"]])))
for (g in names(groups)) {
  members <- groups[[g]]
  if (g %in% rep_groups) {
    n <- min(min_reps, length(members))
    gr <- dba.peakset(obj, peaks = dba.mask(obj, DBA_CONDITION, g), minOverlap = n, bRetrieve = TRUE)
    gr <- granges(gr)
    # 第 4 列：区域中出现的重复数
    n_rep <- Reduce(`+`, lapply(members, function(s) overlapsAny(gr, peaks_of[[s]])))
    to_bed(gr, rep_paths[[paste0(g, ".bed")]], n_rep)
    message(sprintf("%s: %d 个重复，至少 %d 个重复中出现的区域 %d 个", g, length(members), n, length(gr)))
  } else {
    gr <- reduce(peaks_of[[members[1]]])
    message(sprintf("%s: 只有 1 个样本，直接用其 %d 个 peak", g, length(gr)))
  }
  group_sets[[g]] <- gr
}

# 2. 并集
consensus <- reduce(do.call(c, unname(group_sets)))
consensus <- sortSeqlevels(consensus)
consensus <- sort(consensus)
message(sprintf("consensus（各组并集）%d 个区域，中位宽度 %d bp", length(consensus), as.integer(median(width(consensus)))))

# 3. 计数
obj <- dba.count(obj, peaks = consensus, summits = summits, filter = 0, minOverlap = 1,
                 score = DBA_SCORE_READS, bSubControl = FALSE, bRemoveDuplicates = FALSE,
                 bParallel = snakemake@threads > 1)
counts <- dba.peakset(obj, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
# DiffBind 取出的列名会把样本名中的 "-" 改成 "."（GRanges 也一样）；按样本顺序换回原名
stopifnot(ncol(counts) == 3 + length(obj$class[DBA_ID, ]))
names(counts) <- c("Chr", "Start", "End", obj$class[DBA_ID, ])
counts <- counts[order(match(counts$Chr, seqlevels(consensus)), counts$Start), ]
# 输出 consensus.bed 为计数所用的区域（summits 时为峰顶 ± summits bp）
to_bed(GRanges(counts$Chr, IRanges(counts$Start, counts$End)), snakemake@output[["consensus"]])
out <- data.frame(peak = sprintf("%s:%d-%d", counts$Chr, counts$Start, counts$End),
                  counts[c("Chr", "Start", "End")], round(counts[samples]), check.names = FALSE)
write.table(out, snakemake@output[["counts"]], sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(obj, snakemake@output[["rds"]])
message(sprintf("计数完成：%d 个区域，中位宽度 %d bp", nrow(out), as.integer(median(out$End - out$Start + 1))))
