"""
每个样本一行的汇总表。reads 数按"片段"计：PE 为 read pairs，SE 为 reads（flagstat 的 reads 数 PE 除以 2）。
  raw_reads           trim 前（trim_galore 报告 Total reads processed）
  clean_reads         进入 bowtie2 的 reads
  mapping_rate        bowtie2 overall alignment rate
  contamination_rate  比对到污染物种的比例（contamination.enabled 时）
  mito_rate           比对上的 reads 中线粒体的比例
  filtered_reads      MAPQ / properly paired / chrM 过滤后
  dup_rate            PCR 重复比例（相对 filtered_reads）
  final_reads         去重后，用于 call peak 和 bigwig
  peaks / frip        blacklist 过滤后的 peak 数与 FRiP（control 样本为空）
  tss_enrichment      TSS 中心信号 / 两端 flank 信号（deeptools profile）
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


def bowtie2_stats(path):
    t = Path(path).read_text()
    return grab(r"(\d+) reads; of these:", t), grab(r"([\d.]+)% overall alignment rate", t, float)


def mito_rate(path, mito):
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "len", "mapped", "unmapped"])
    total = df["mapped"].sum()
    return df.loc[df["chrom"] == mito, "mapped"].sum() / total if total else None


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
    contam[Path(fp).stem] = bowtie2_stats(fp)[1]

rows = []
for i, s in enumerate(snakemake.params.samples):
    m = meta[s]
    div = 2 if m["layout"] == "PE" else 1
    clean, rate = bowtie2_stats(inp.bowtie2[i])
    filt = flagstat_reads(inp.filtered[i])
    final = flagstat_reads(inp.final[i])
    row = {
        "sample": s,
        "role": m["role"],
        "layout": m["layout"],
        "control": m["control"],
        "peak_type": m["peak_type"],
        "raw_reads": grab(r"Total reads processed:\s+([\d,]+)", Path(inp.trim[i]).read_text()),
        "clean_reads": clean,
        "mapping_rate": rate / 100 if rate is not None else None,
        "contamination_rate": contam[s] / 100 if s in contam else None,
        "mito_rate": mito_rate(inp.idxstats[i], snakemake.params.mito),
        "filtered_reads": filt // div if filt is not None else None,
        "dup_rate": 1 - final / filt if filt else None,
        "final_reads": final // div if final is not None else None,
        "peaks": peaks.get(s),
        "frip": frip.get(s),
        "tss_enrichment": tss.get(s),
    }
    rows.append(row)

df = pd.DataFrame(rows)
if not contam:
    df = df.drop(columns=["contamination_rate"])
if not tss:
    df = df.drop(columns=["tss_enrichment"])
for c in df.columns:
    if c.endswith(("_rate", "frip", "enrichment")):
        df[c] = df[c].astype(float).round(4)
    elif c.endswith("_reads") or c == "peaks":
        df[c] = df[c].astype("Int64")
df.to_csv(snakemake.output.tsv, sep="\t", index=False)

with open(snakemake.output.mqc, "w") as f:
    f.write("# id: 'pipeline_summary'\n# section_name: 'Pipeline summary'\n"
            "# description: 'results/tables/summary.tsv（PE 按 read pairs 计）'\n"
            "# plot_type: 'table'\n")
    df.to_csv(f, sep="\t", index=False)
