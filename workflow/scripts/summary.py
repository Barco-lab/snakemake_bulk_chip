"""
每个样本一行的汇总表。reads 数按单条 read 计（PE 为 R1 + R2），所有 rate 都以 raw_reads 为分母（%）。
  raw_reads           trim 前（trim_galore 报告 Total reads processed，PE 两个报告相加）
  clean_reads         进入 bowtie2 的 reads
  mapped_reads        比对上的 reads（raw BAM flagstat primary mapped）
  filtered_reads      MAPQ / properly paired / chrM 过滤后
  dedup_reads         去重后，用于 call peak 和 bigwig
  mapping_rate        mapped_reads / raw_reads
  dedup_rate          dedup_reads / raw_reads
  duplicates_rate     (filtered_reads - dedup_reads) / raw_reads
  mito_rate           raw BAM 中线粒体 reads / raw_reads
  peaks / frip        blacklist 过滤后的 peak 数与 FRiP（peak 内 reads / dedup BAM reads；control 样本为空）
  tss_enrichment      TSS 中心信号 / 两端 flank 信号（deeptools profile）
  contamination_rate  比对到污染物种的 reads / raw_reads（contamination.enabled 时）
表头括号内为对应的 config 值。summary.tsv 的 reads 数带千位分隔符，MultiQC 用的表保持纯数字。
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


def tss_enrichment(path, flank):
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
        center = max(vals[len(vals) // 2 - 1: len(vals) // 2 + 1])
        out[f[0]] = center / bg if bg > 0 else None
    return out


meta = snakemake.params.meta
inp = snakemake.input
tss = tss_enrichment(inp.get("tss"), snakemake.params.tss_flank) if "tss" in inp.keys() else {}
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
C_MITO = f"mito_rate ({snakemake.params.mito}) (%)"
C_TSS = f"tss_enrichment (±{snakemake.params.tss_flank}bp)"
C_CONTAM = "contamination_rate (%)"

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
    if s in contam and contam[s][0] is not None:
        n, rate = contam[s]
        contam_reads = n * div * rate / 100
    else:
        contam_reads = None
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
        "mapping_rate (%)": pct(mapped, raw),
        "dedup_rate (%)": pct(dedup, raw),
        "duplicates_rate (%)": pct(filt - dedup, raw) if filt is not None and dedup is not None else None,
        C_MITO: pct(mito_reads(inp.idxstats[i], snakemake.params.mito), raw),
        "peaks": peaks.get(s),
        "frip": frip.get(s),
        C_TSS: tss.get(s),
        C_CONTAM: pct(contam_reads, raw),
    }
    rows.append(row)

df = pd.DataFrame(rows)
if not contam:
    df = df.drop(columns=[C_CONTAM])
if not tss:
    df = df.drop(columns=[C_TSS])
count_cols = [c for c in df.columns if "_reads" in c or c == "peaks"]
for c in df.columns:
    if c.endswith("(%)"):
        df[c] = df[c].astype(float).round(2)
    elif c == "frip" or c == C_TSS:
        df[c] = df[c].astype(float).round(4)
    elif c in count_cols:
        df[c] = df[c].astype("Int64")

out = df.copy()
for c in count_cols:
    out[c] = out[c].map(lambda v: "" if pd.isna(v) else f"{int(v):,}")
out.to_csv(snakemake.output.tsv, sep="\t", index=False)

with open(snakemake.output.mqc, "w") as f:
    f.write("# id: 'pipeline_summary'\n# section_name: 'Pipeline summary'\n"
            "# description: 'results/tables/summary.tsv（reads 按单条计，rate 以 raw_reads 为分母）'\n"
            "# plot_type: 'table'\n")
    df.to_csv(f, sep="\t", index=False)
