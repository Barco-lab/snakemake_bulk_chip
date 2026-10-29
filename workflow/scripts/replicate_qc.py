"""
重复一致性：
  peak 重现性  组内重复两两比较：Jaccard（bedtools jaccard，按碱基）与 A 的 peak 中与 B 重叠的比例；
              组内各重复组合独有的碱基数（bedtools multiinter，供 Venn 图）；组间同理，用各组的 reproducible peak
  相关性 / PCA consensus peak 计数（DiffBind，见 diffbind_consensus.R）-> log2(CPM + 1)（CPM 以 consensus peak 内的总计数为分母）
  性别        chrX、chrY 计数区（y_region，chrY 去掉 PAR）的 reads 占全部比对 reads 的百分比；按 chrY 占比判定
ChIP 流程只比较 IP 样本，被用作 control 的样本只做性别检查（samples.tsv 中 role 为 control，重复一致性各列为空）；ATAC 没有对照，role 均为 IP。
结果写成 tsv / json，供 replicate_plots.R 绘图和 qc_report.py 使用。
"""
import itertools
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

p = snakemake.params
out = snakemake.output
samples, groups = list(p.samples), dict(p.groups)
group_of = {s: g for g, ss in groups.items() for s in ss}
peak_of = dict(zip(samples, snakemake.input.peaks))


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


# ------------------------------------------------------------
# peak 重现性
# ------------------------------------------------------------
tmp = Path(tempfile.mkdtemp(dir=Path(out.info).parent))
sorted_peak, n_peaks = {}, {}
for s in samples:
    f = tmp / f"{s}.bed"
    sh(f"cut -f1-3 {peak_of[s]} | sort -k1,1 -k2,2n > {f}")
    sorted_peak[s] = f
    n_peaks[s] = sum(1 for _ in open(f))


def frac_in(a, b):
    """a 的 peak 中与 b 的任一 peak 重叠（≥ 1 bp）的比例。"""
    if n_peaks[a] == 0:
        return np.nan
    return int(sh(f"bedtools intersect -u -a {sorted_peak[a]} -b {sorted_peak[b]} | wc -l")) / n_peaks[a]


rows = []
for g in p.rep_groups:
    for a, b in itertools.combinations(groups[g], 2):
        if n_peaks[a] and n_peaks[b]:
            jac = float(sh(f"bedtools jaccard -a {sorted_peak[a]} -b {sorted_peak[b]}").splitlines()[1].split("\t")[2])
        else:
            jac = np.nan
        rows.append(dict(group=g, sample_a=a, sample_b=b, jaccard=jac,
                         frac_a_in_b=frac_in(a, b), frac_b_in_a=frac_in(b, a)))
overlap = pd.DataFrame(rows, columns=["group", "sample_a", "sample_b", "jaccard", "frac_a_in_b", "frac_b_in_a"])
overlap.to_csv(out.overlap, sep="\t", index=False, float_format="%.4f")

# Venn：multiinter 把组内各重复的 peak 切成按覆盖它的重复组合标记的区段，按组合累加碱基数；与输入顺序无关。
# 每个重复先自身合并：本流程的 MACS3 参数下 peak 互不重叠，此步为防护（如加了 --call-summits，
# 同一区域会按 summit 重复出现，不合并会重复计数）
merged = {}
for s in samples:
    merged[s] = tmp / f"{s}.merged.bed"
    sh(f"bedtools merge -i {sorted_peak[s]} > {merged[s]}")


def combo_bp(names, beds):
    """multiinter：各集合组合独有的碱基数，[(组合, 集合数, bp)]，含 0。"""
    bp = {}
    files = " ".join(str(b) for b in beds)
    for line in sh(f"bedtools multiinter -i {files} -names {' '.join(names)}").splitlines():
        f = line.split("\t")
        key = frozenset(f[4].split(","))
        bp[key] = bp.get(key, 0) + int(f[2]) - int(f[1])
    return [(combo, k, bp.get(frozenset(combo), 0))
            for k in range(1, len(names) + 1) for combo in itertools.combinations(names, k)]


venn = []
for g in p.rep_groups:
    ss = groups[g]
    for combo, k, bp in combo_bp(ss, [merged[s] for s in ss]):
        venn.append(dict(group=g, samples=",".join(combo), n=k, bp=bp))
pd.DataFrame(venn, columns=["group", "samples", "n", "bp"]).to_csv(out.venn, sep="\t", index=False)

# 组间 Venn：各组的 reproducible peak（只有 1 个样本的组用该样本的 peak），同样按碱基
repro_file = {Path(f).stem: f for f in snakemake.input.reproducible}
group_bed = {}
for g in {g for gs in p.group_venn.values() for g in gs}:
    if g in repro_file:
        group_bed[g] = tmp / f"group.{g}.bed"
        sh(f"cut -f1-3 {repro_file[g]} | sort -k1,1 -k2,2n | bedtools merge > {group_bed[g]}")
    else:
        group_bed[g] = merged[groups[g][0]]
gvenn = []
for name, gs in p.group_venn.items():
    for combo, k, bp in combo_bp(gs, [group_bed[g] for g in gs]):
        gvenn.append(dict(comparison=name, groups=",".join(combo), n=k, bp=bp))
pd.DataFrame(gvenn, columns=["comparison", "groups", "n", "bp"]).to_csv(out.group_venn, sep="\t", index=False)
for f in tmp.iterdir():
    f.unlink()
tmp.rmdir()

reproducible = {Path(f).stem: sum(1 for _ in open(f)) for f in snakemake.input.reproducible}

# ------------------------------------------------------------
# consensus peak 计数 -> log2(CPM + 1) -> 相关性、PCA
# ------------------------------------------------------------
counts = pd.read_csv(snakemake.input.counts, sep="\t", index_col=0)[samples]

lib = counts.sum()
logc = np.log2(counts / lib.replace(0, np.nan) * 1e6 + 1).fillna(0)
cor = logc.corr(method=p.method)
cor.to_csv(out.correlation, sep="\t", float_format="%.4f")
logc.to_csv(out.log2cpm, sep="\t", float_format="%.3f")
scatter_peaks = min(20000, len(logc))  # 散点矩阵作图的 peak 数（随机抽取；r 仍基于全部 peak）

top = min(p.pca_top, len(logc))
pca = pd.DataFrame({"sample": samples, "group": [group_of[s] for s in samples]})
var_exp = []
if len(samples) >= 3 and top >= 2:
    sel = logc.loc[logc.var(axis=1).sort_values(ascending=False).index[:top]]
    x = sel.T.values - sel.T.values.mean(axis=0)
    u, d_, _ = np.linalg.svd(x, full_matrices=False)
    k = min(3, len(d_))
    for i in range(k):
        pca[f"PC{i + 1}"] = u[:, i] * d_[i]
    var_exp = (d_ ** 2 / (d_ ** 2).sum())[:k].round(4).tolist()
pca.to_csv(out.pca, sep="\t", index=False, float_format="%.4f")

# ------------------------------------------------------------
# 每个样本的指标
# ------------------------------------------------------------
per = []
for s in samples:
    g = group_of[s]
    mates = [x for x in groups[g] if x != s]
    others = [x for x in samples if group_of[x] != g]
    r = dict(sample=s, group=g, peaks=n_peaks[s], reads_in_consensus=int(lib[s]))
    r["within_cor"] = cor.loc[s, mates].mean() if mates else np.nan
    if others:
        o = cor.loc[s, others].astype(float)
        r["max_other_cor"], r["max_other_sample"] = o.max(), o.idxmax()
        r["max_other_group"] = group_of[o.idxmax()]
    else:
        r["max_other_cor"], r["max_other_sample"], r["max_other_group"] = np.nan, "", ""
    sub = overlap[(overlap.sample_a == s) | (overlap.sample_b == s)]
    r["mean_jaccard"] = sub["jaccard"].mean() if len(sub) else np.nan
    per.append(r)
per = pd.DataFrame(per)
per.insert(2, "role", "IP")

lo, hi = p.y_frac
if snakemake.input.sex:
    sex = pd.concat(pd.read_csv(f, sep="\t") for f in snakemake.input.sex)
    sex["x_frac"] = sex.x_reads / sex.total_reads * 100
    sex["y_frac"] = sex.y_reads / sex.total_reads * 100
    sex["sex"] = np.where(sex.y_frac < lo, "F", np.where(sex.y_frac > hi, "M", "undetermined"))
    per = per.merge(sex[["sample", "total_reads", "x_reads", "y_reads", "x_frac", "y_frac", "sex"]], on="sample", how="outer")
    ctrl = per["role"].isna()
    per.loc[ctrl, "group"] = per.loc[ctrl, "sample"].map(p.group_of)
    per.loc[ctrl, "role"] = "control"
    # IP 在前、control 在后，各自按样本表顺序
    order = {s: i for i, s in enumerate(p.group_of)}
    per = per.iloc[sorted(range(len(per)), key=lambda i: (per["role"].iat[i] == "control", order[per["sample"].iat[i]]))]
per.to_csv(out.samples, sep="\t", index=False, float_format="%.6g")

# ------------------------------------------------------------
# 每组
# ------------------------------------------------------------
grows = []
for g, ss in groups.items():
    sub = overlap[overlap.group == g]
    grows.append(dict(
        group=g, n_samples=len(ss), samples=",".join(ss),
        median_peaks=float(np.median([n_peaks[s] for s in ss])),
        reproducible_peaks=reproducible.get(g, np.nan),
        median_jaccard=sub["jaccard"].median() if len(sub) else np.nan,
        min_jaccard=sub["jaccard"].min() if len(sub) else np.nan,
        mean_within_cor=per.loc[per.group == g, "within_cor"].mean() if len(ss) > 1 else np.nan,
    ))
pd.DataFrame(grows).to_csv(out.groups, sep="\t", index=False, float_format="%.4f")

json.dump(dict(
    method=p.method, min_reps=p.min_reps, consensus_peaks=int(len(counts)), summits=p.summits, pca_peaks=int(top),
    pca_variance=var_exp, sex_check=bool(snakemake.input.sex), y_frac=[lo, hi], y_region=p.y_region, genome=p.genome, mapq=p.mapq,
    groups=groups, rep_groups=list(p.rep_groups), group_venn=dict(p.group_venn), scatter_peaks=int(scatter_peaks),
), open(out.info, "w"), indent=1)
print(f"{len(samples)} 个样本，{len(groups)} 组，consensus peak {len(counts)}", file=sys.stderr)
