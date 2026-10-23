"""
每个样本一行的汇总表。reads 数按单条 read 计（PE 为 R1 + R2）。
rate 以进入该步骤的 reads 为分母（与各软件报告、ENCODE 标准口径一致），分母写在表头中：
  raw_reads           trim 前（trim_galore 报告 Total reads processed，PE 两个报告相加）
  clean_reads         进入 bowtie2 的 reads
  mapped_reads        比对上的 reads（raw BAM flagstat primary mapped）
  filtered_reads      MAPQ / properly paired 过滤后（filter.remove_mito 时另去掉线粒体）
  dedup_reads         去重后，用于 call peak 和 bigwig（PE 除以 2 即 ENCODE 的 usable fragments）
  clean_rate          clean_reads / raw_reads
  mapping_rate        mapped_reads / clean_reads
  mito_rate           raw BAM 中线粒体 reads / mapped_reads
  dedup_rate          dedup_reads / filtered_reads
  duplicates_rate     (filtered_reads - dedup_reads) / filtered_reads
  peaks / frip        blacklist 过滤后的 peak 数与 FRiP（peak 内 reads / dedup_reads；control 样本为空）
  tss_enrichment      TSS 中心信号 / 两端 flank 信号（deeptools profile）
  contamination_rate  trim 后 reads 比对到污染物种的比例（contamination.enabled 时）
表头括号内为对应的 config 值与分母。输出纯数字，MultiQC 与 qc_report 读取。
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def grab(pattern, text, cast=int):
    m = re.search(pattern, text)
    return cast(m.group(1).replace(",", "")) if m else None


def flagstat_reads(path):
    return grab(r"(\d+) \+ \d+ in total", Path(path).read_text())


def flagstat_mapped(path):
    return grab(r"(\d+) \+ \d+ primary mapped", Path(path).read_text())


def bowtie2_stats(path):
    t = Path(path).read_text()
    return grab(r"(\d+) reads; of these:", t), grab(r"([\d.]+)% overall alignment rate", t, float)


def mito_reads(path, mito):
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "len", "mapped", "unmapped"])
    return int(df.loc[df["chrom"] == mito, "mapped"].sum())


def pct(n, total):
    return 100 * n / total if n is not None and total else None


def matrix_window(cmd):
    """从 computeMatrix 命令解析 -b / -a / --binSize（未写时用 deeptools 默认值 500 / 500 / 10）。"""
    def opt(names, default):
        m = re.search(r"(?:^|\s)(?:" + "|".join(names) + r")[\s=]+(\d+)", cmd)
        return int(m.group(1)) if m else default
    return (
        opt(["-b", "--upstream", "--beforeRegionStartLength"], 500),
        opt(["-a", "--downstream", "--afterRegionStartLength"], 500),
        opt(["-bs", "--binSize"], 10),
    )


def tss_enrichment(path, up, bin_size):
    """plotProfile --outFileNameData：每行一个样本，前两列为 label 与 group，之后为各 bin 的均值。"""
    if not path:
        return {}
    out = {}
    lines = Path(path).read_text().splitlines()
    for line in lines[2:]:
        f = line.split("\t")
        vals = [float(x) for x in f[2:] if x != ""]
        if len(vals) < 5:
            continue
        edge = max(1, len(vals) // 20)  # 两端各 5% 的 bin 作为背景
        bg = (sum(vals[:edge]) + sum(vals[-edge:])) / (2 * edge)
        tss = up // bin_size  # TSS 右侧第一个 bin 的下标
        center = max(vals[tss - 1: tss + 1])
        out[f[0]] = center / bg if bg > 0 else None
    return out


meta = snakemake.params.meta
inp = snakemake.input
UP, DOWN, BIN = matrix_window(snakemake.params.tss_matrix)
tss = tss_enrichment(inp.get("tss"), UP, BIN) if "tss" in inp.keys() else {}
frip = {}
for fp in inp.frip:
    r = pd.read_csv(fp, sep="\t").iloc[0]
    frip[r["sample"]] = r["frip"]
peaks = {s: sum(1 for _ in open(p)) for s, p in zip(snakemake.params.peak_samples, inp.peaks)}
contam = {}
for fp in inp.contam:
    contam[Path(fp).stem] = bowtie2_stats(fp)

# 表头括号内为对应的 config 值
C_FILT = f"filtered_reads (q{snakemake.params.mapq})"
C_MITO = f"mito_rate ({snakemake.params.mito}) (% of mapped)"
C_TSS = f"tss_enrichment (±{UP}bp)" if UP == DOWN else f"tss_enrichment (-{UP}/+{DOWN}bp)"
C_CONTAM = "contamination_rate (% of clean)"

rows = []
for i, s in enumerate(snakemake.params.samples):
    m = meta[s]
    div = 2 if m["layout"] == "PE" else 1
    raw = sum(grab(r"Total reads processed:\s+([\d,]+)", Path(fp).read_text()) for fp in snakemake.params.trim[s])
    clean = bowtie2_stats(inp.bowtie2[i])[0]
    clean = clean * div if clean is not None else None
    mapped = flagstat_mapped(inp.raw[i])
    filt = flagstat_reads(inp.filtered[i])
    dedup = flagstat_reads(inp.final[i])
    row = {
        "sample": s,
        "role": m["role"],
        "layout": m["layout"],
        "control": m["control"],
        "peak_type": m["peak_type"],
        "raw_reads": raw,
        "clean_reads": clean,
        "mapped_reads": mapped,
        C_FILT: filt,
        "dedup_reads": dedup,
        "clean_rate (% of raw)": pct(clean, raw),
        "mapping_rate (% of clean)": pct(mapped, clean),
        C_MITO: pct(mito_reads(inp.idxstats[i], snakemake.params.mito), mapped),
        "dedup_rate (% of filtered)": pct(dedup, filt),
        "duplicates_rate (% of filtered)": pct(filt - dedup, filt) if filt is not None and dedup is not None else None,
        "peaks": peaks.get(s),
        "frip": frip.get(s),
        C_TSS: tss.get(s),
        C_CONTAM: contam[s][1] if s in contam else None,
    }
    rows.append(row)

df = pd.DataFrame(rows)
if not contam:
    df = df.drop(columns=[C_CONTAM])
if not tss:
    df = df.drop(columns=[C_TSS])
count_cols = [c for c in df.columns if "_reads" in c or c == "peaks"]
for c in df.columns:
    if "(%" in c:
        df[c] = df[c].astype(float).round(2)
    elif c == "frip" or c == C_TSS:
        df[c] = df[c].astype(float).round(4)
    elif c in count_cols:
        df[c] = df[c].astype("Int64")

df.to_csv(snakemake.output[0], sep="\t", index=False)
