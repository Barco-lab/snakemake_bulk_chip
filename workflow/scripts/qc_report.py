"""
单文件 HTML QC 报告：概览、流程与参数、QC 汇总表、重复一致性（中英文可切换）。
所有样式、脚本、流程图和图片都内嵌，不依赖网络。各列说明、参考值和检查规则来自 workflow/report/qc_metrics.yaml：
  checks  有公开阈值的检查（warn 橙 / fail 红）
  outlier 与本批样本比较，|值 - 中位数| > 3 x MAD（x 1.4826）标为橙色
重复一致性的数据与图来自 results/qc/replicates/（replicate_qc.py 计算，replicate_plots.R 绘图）。
"""
import base64
import html
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

sys.stderr = open(snakemake.log[0], "w")

p = snakemake.params
df = pd.read_csv(snakemake.input.summary, sep="\t")
meta = yaml.safe_load(Path(snakemake.input.metrics).read_text())
METRICS, SOURCES = meta["metrics"], meta.get("sources", {})


def base(col):
    """summary.tsv 列名去掉括号部分，如 'mapping_rate (% of clean)' -> 'mapping_rate'。"""
    return col.split(" (")[0]


def esc(x):
    return html.escape(str(x))


def rich(text):
    """转义后把 `...` 显示为代码。"""
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", esc(text))


def L(zh, en=None):
    """中英文各一份，按当前语言显示其一。"""
    en = zh if en is None else en
    if zh == en:
        return rich(zh)
    return f'<span data-l="zh">{rich(zh)}</span><span data-l="en">{rich(en)}</span>'


def T(zh, en):
    """悬停说明（title），切换语言时由脚本替换。"""
    return f' title="{esc(zh)}" data-tz="{esc(zh)}" data-te="{esc(en)}"'


def zh_en(x):
    return (x["zh"], x["en"]) if isinstance(x, dict) else (x, x)


def is_num(col):
    return pd.api.types.is_numeric_dtype(df[col])


def fmt(col, v):
    """reads 等计数加千位分隔符；其余数值按 qc_metrics.yaml 的 digits（默认 2）位小数。"""
    if pd.isna(v):
        return ""
    if base(col).endswith("_reads") or pd.api.types.is_integer_dtype(df[col]):
        return f"{int(v):,}"
    if is_num(col):
        return f"{v:.{METRICS.get(base(col), {}).get('digits', 2)}f}"
    return str(v)


cols = {base(c): c for c in df.columns}
col = lambda key: df[cols[key]] if key in cols else None
n = len(df)

# ------------------------------------------------------------
# 检查：flags[(行号, 列名)] = [(level, 中文, 英文)]
# ------------------------------------------------------------
flags = {}


def flag(i, c, level, zh, en):
    flags.setdefault((i, c), []).append((level, zh, en))


for i, row in df.iterrows():
    env = {base(c): (None if pd.isna(v) else v) for c, v in row.items()}
    for key, m in METRICS.items():
        if key not in cols:
            continue
        for chk in m.get("checks", []):
            if any(env.get(k) != v for k, v in chk.get("when", {}).items()):
                continue
            try:
                hit = eval(chk["expr"], {}, env)
            except TypeError:  # 缺失值
                continue
            if hit:
                flag(i, cols[key], chk["level"], *zh_en(chk["message"]))

outliers = []  # (样本, 列, 高/低, 中位数)
for key, m in METRICS.items():
    if not m.get("outlier") or key not in cols:
        continue
    sel = df
    for k, v in m.get("outlier_when", {}).items():  # 只在部分样本之间比较
        sel = sel[sel[cols[k]] == v]
    s = sel[cols[key]].dropna().astype(float)
    if len(s) < 4:
        continue
    med = s.median()
    mad = 1.4826 * (s - med).abs().median()
    if mad == 0:
        continue
    for i, v in s.items():
        if abs(v - med) > 3 * mad:
            hi = v > med
            outliers.append((df.loc[i, "sample"], key, hi, fmt(cols[key], med)))
            mt = fmt(cols[key], med)
            flag(i, cols[key], "warn",
                 f"明显{'高' if hi else '低'}于本批中位数 {mt}（|差| > 3 × MAD）",
                 f"{'above' if hi else 'below'} the batch median {mt} (|diff| > 3 × MAD)")

# ------------------------------------------------------------
# 重复一致性的检查（均为 warn：没有公开阈值）
# ------------------------------------------------------------
inp = snakemake.input
REP = "rep_info" in inp.keys()
REP_COLS = meta.get("replicates", {})
if REP:
    rinfo = json.loads(Path(inp.rep_info).read_text())
    rs = pd.read_csv(inp.rep_samples, sep="\t")
    rg = pd.read_csv(inp.rep_groups, sep="\t")
    row_of = {s: i for i, s in df["sample"].items()}
    reps = rs[rs["within_cor"].notna()]

    for _, r in reps.iterrows():
        if r.max_other_cor > r.within_cor:
            flag(row_of[r["sample"]], "within_cor", "warn",
                 f"与 {r.max_other_group} 组的 {r.max_other_sample} 相关性（{r.max_other_cor:.3f}）高于组内平均（{r.within_cor:.3f}）："
                 "两组差异小时属正常，否则请确认样本标签",
                 f"correlates more with {r.max_other_sample} of {r.max_other_group} ({r.max_other_cor:.3f}) than with its own "
                 f"replicates on average ({r.within_cor:.3f}): expected when the groups differ little, otherwise check the sample labels")

    for key, zh, en in (("within_cor", "组内相关性", "Within-group correlation"),
                        ("mean_jaccard", "与同组重复的 peak 重叠（Jaccard）", "Peak overlap (Jaccard) with its replicates")):
        s = reps[key].dropna()
        if len(s) < 4:
            continue
        m = s.median()
        mad = 1.4826 * (s - m).abs().median()
        for j, v in s.items():
            if mad > 0 and m - v > 3 * mad:
                flag(row_of[reps.loc[j, "sample"]], key, "warn",
                     f"{zh}明显低于本批中位数 {m:.3f}（差 > 3 × MAD）", f"{en} well below the batch median {m:.3f} (> 3 × MAD)")

    if "sex" in rs:
        lo, hi = rinfo["y_frac"]
        for _, r in rs[rs["sex"] == "undetermined"].iterrows():
            flag(row_of[r["sample"]], "sex", "warn",
                 f"chrY 占比 {r.y_frac:.4g}%，介于阈值 {lo}%–{hi}% 之间，无法判断",
                 f"chrY fraction {r.y_frac:.4g}%, between the thresholds {lo}%–{hi}%; undetermined")
        for g, sub in rs[rs["sex"].isin(["F", "M"])].groupby("group"):
            vc = sub["sex"].value_counts()
            if len(vc) < 2:
                continue
            counts = "、".join(f"{k} {v}" for k, v in vc.items())
            for _, r in sub[sub["sex"].map(vc) <= len(sub) / 2].iterrows():
                flag(row_of[r["sample"]], "sex", "warn",
                     f"推断为 {r.sex}（chrY 占比 {r.y_frac:.4g}%），与同组其他样本不同（{counts}）",
                     f"inferred {r.sex} (chrY fraction {r.y_frac:.4g}%), unlike the rest of the group ({counts.replace('、', ', ')})")
        # chrX 占比随细胞类型变化，只在组内比较：雄性只有一条 X，占比应低于同组雌性
        for g, sub in rs[rs["x_frac"].notna()].groupby("group"):
            if len(sub) < 3:
                continue
            m = sub["x_frac"].median()
            mad = 1.4826 * (sub["x_frac"] - m).abs().median()
            for _, r in sub.iterrows():
                if mad > 0 and m - r.x_frac > 3 * mad:
                    flag(row_of[r["sample"]], "x_frac", "warn",
                         f"chrX 占比 {r.x_frac:.2f}% 明显低于同组中位数 {m:.2f}%（差 > 3 × MAD），可能是雄性或混入了雄性样本",
                         f"chrX fraction {r.x_frac:.2f}% well below the group median {m:.2f}% (> 3 × MAD); possibly male or mixed with a male sample")

flagged = sorted({i for i, _ in flags})

shown = list(df.columns)


# ------------------------------------------------------------
# 概览
# ------------------------------------------------------------
def med(key, scale=1, digits=1, unit=""):
    s = col(key)
    if s is None or s.dropna().empty:
        return None
    return f"{s.median() / scale:.{digits}f}{unit}"


ip = df[col("role") == "chip"] if "role" in cols else df
ctrl = df[col("role") == "control"] if "role" in cols else df.iloc[0:0]


def med_of(sub, key, scale=1, digits=1, unit=""):
    s = sub[cols[key]].dropna() if key in cols else pd.Series(dtype=float)
    return None if s.empty else f"{s.median() / scale:.{digits}f}{unit}"


tiles = [
    (("样本", "Samples"), f"{len(ip)} + {len(ctrl)}", ("IP + 对照", "IP + control")),
    (("raw reads 中位数", "Median raw reads"), med("raw_reads", 1e6, 1, "M"), ("", "")),
    (("mapping rate 中位数", "Median mapping rate"), med("mapping_rate", 1, 1, "%"), ("以 clean reads 为分母", "of clean reads")),
    (("dedup reads 中位数", "Median dedup reads"), med_of(ip, "dedup_reads", 1e6, 1, "M"), ("IP 样本", "IP samples")),
    (("重复率中位数", "Median duplicate rate"), med_of(ip, "duplicates_rate", 1, 1, "%"), ("IP 样本，以 filtered 为分母", "IP samples, of filtered")),
    (("FRiP 中位数", "Median FRiP"), med_of(ip, "frip", 1, 3), ("IP 样本", "IP samples")),
    (("peak 数中位数", "Median peaks"), med_of(ip, "peaks", 1e3, 1, "k"), ("IP 样本", "IP samples")),
    (("组内相关性中位数", "Median within-group r"),
     f"{rs['within_cor'].median():.3f}" if REP and rs["within_cor"].notna().any() else None, ("IP 样本重复一致性", "IP replicate consistency")),
    (("需关注的样本", "Samples to review"), str(len(flagged)), ("见下方列表", "listed below")),
]
tiles_html = "".join(
    f'<div class="tile{" tile-alert" if label[0] == "需关注的样本" and flagged else ""}">'
    f'<div class="tile-label">{L(*label)}</div><div class="tile-value">{esc(value)}</div>'
    f'<div class="tile-note">{L(*note)}</div></div>'
    for label, value, note in tiles
    if value is not None
)

if flagged:
    items = []
    for i in flagged:
        chips = "".join(
            f'<span class="chip chip-{lvl}">{esc(base(c))}：{L(zh, en)}</span>'
            for (j, c), lst in flags.items()
            if j == i
            for lvl, zh, en in lst
        )
        items.append(f'<tr><td class="sample">{esc(df.loc[i, "sample"])}</td><td>{chips}</td></tr>')
    attention = f'<table class="plain"><tbody>{"".join(items)}</tbody></table>'
else:
    attention = f'<p class="ok-note">{L("所有样本均未触发检查，也没有离群值。", "No sample triggered a check or was an outlier.")}</p>'


# ------------------------------------------------------------
# 要点：根据数据自动生成，放在 QC 汇总表下方
# ------------------------------------------------------------
def rng(s, scale=1, digits=1, unit=""):
    s = s.dropna() / scale
    return f"{s.median():.{digits}f}{unit}", f"{s.min():.{digits}f}–{s.max():.{digits}f}{unit}"


def names(lst):
    return ", ".join(lst)


points = []  # (中文, 英文)
raw_m, raw_r = rng(col("raw_reads"), 1e6, 1, "M")
zh, en = f"测序深度：raw reads 中位数 {raw_m}（{raw_r}）。", f"Depth: median raw reads {raw_m} ({raw_r})."
for ptype, need in (("narrow", 20e6), ("broad", 45e6)):
    sub = ip[ip[cols["peak_type"]] == ptype] if "peak_type" in cols else ip.iloc[0:0]
    if sub.empty:
        continue
    frag = sub[cols["dedup_reads"]] / sub[cols["layout"]].map({"PE": 2, "SE": 1})
    f_m, f_r = rng(frag, 1e6, 1, "M")
    low = sub.loc[frag < need, "sample"].tolist()
    zh += (f" {ptype} IP 样本 usable fragments 中位数 {f_m}（{f_r}），"
           + (f"均达到 ENCODE 建议的 {need / 1e6:.0f}M。" if not low else f"{len(low)}/{len(sub)} 个低于 ENCODE 建议的 {need / 1e6:.0f}M（{names(low)}）；CUT&Tag 所需深度远低于此。"))
    en += (f" {ptype} IP samples: median usable fragments {f_m} ({f_r}), "
           + (f"all meeting the ENCODE {need / 1e6:.0f}M." if not low else f"{len(low)}/{len(sub)} below the ENCODE {need / 1e6:.0f}M ({names(low)}); CUT&Tag needs far less."))
if len(ctrl):
    c_m, c_r = rng(ctrl[cols["dedup_reads"]], 1e6, 1, "M")
    zh += f" 对照样本 dedup reads 中位数 {c_m}（{c_r}），对照深度应不低于对应的 IP 样本。"
    en += f" Controls: median dedup reads {c_m} ({c_r}); controls should be at least as deep as their IP samples."
points.append((zh, en))

m_m, m_r = rng(col("mapping_rate"))
mt_m, mt_r = rng(col("mito_rate"))
dp_m, dp_r = rng(ip[cols["duplicates_rate"]])
points.append((
    f"比对与过滤：mapping rate 中位数 {m_m}%（{m_r}%），线粒体比例中位数 {mt_m}%（{mt_r}%），IP 样本重复率中位数 {dp_m}%（{dp_r}%）。"
    "这几项 ENCODE 没有直接阈值，主要看同批是否有离群样本；CUT&Tag 的重复率通常偏高。",
    f"Alignment and filtering: median mapping rate {m_m}% ({m_r}%), mitochondrial fraction {mt_m}% ({mt_r}%), IP duplicate rate {dp_m}% ({dp_r}%). "
    "ENCODE gives no direct cutoffs for these; look for outliers. CUT&Tag usually shows higher duplication.",
))

if len(ip) and "frip" in cols:
    f_m, f_r = rng(ip[cols["frip"]], 1, 3)
    low = ip.loc[ip[cols["frip"]] < 0.01, "sample"].tolist()
    pk_m, pk_r = rng(ip[cols["peaks"]], 1e3, 1, "k")
    zh = (f"信号：IP 样本 FRiP 中位数 {f_m}（{f_r}），"
          + ("均高于 ENCODE 指南的 1%" if not low else f"{len(low)} 个低于 1%（{names(low)}）")
          + f"；peak 数中位数 {pk_m}（{pk_r}）。FRiP 和 peak 数随靶标差别很大，应与同一靶标的样本比较。")
    en = (f"Signal: median IP FRiP {f_m} ({f_r}), "
          + ("all above the 1% of the ENCODE guidelines" if not low else f"{len(low)} below 1% ({names(low)})")
          + f"; median peaks {pk_m} ({pk_r}). Both depend strongly on the target, so compare samples of the same target.")
    if "tss_enrichment" in cols and len(ctrl):
        ti, tc = ip[cols["tss_enrichment"]].median(), ctrl[cols["tss_enrichment"]].median()
        zh += f" TSS enrichment 中位数 IP {ti:.1f}、对照 {tc:.1f}；启动子标记应明显高于对照。"
        en += f" Median TSS enrichment: IP {ti:.1f}, control {tc:.1f}; promoter marks should be well above the control."
    points.append((zh, en))

if "contamination_rate" in cols:
    c_m, c_r = rng(col("contamination_rate"))
    points.append((f"污染：clean reads 比对到污染物种的比例中位数 {c_m}%（{c_r}%）。",
                   f"Contamination: median fraction of clean reads aligning to the contaminant genome {c_m}% ({c_r}%)."))

if outliers:
    zh = "离群值（与本批中位数相差超过 3 × MAD）：" + "；".join(
        f"{s} 的 {k}{'偏高' if hi else '偏低'}（中位数 {m}）" for s, k, hi, m in outliers) + "。"
    en = "Outliers (> 3 × MAD from the batch median): " + "; ".join(
        f"{s} {k} {'high' if hi else 'low'} (median {m})" for s, k, hi, m in outliers) + "."
    points.append((zh, en))

FAIL_NOTE = ("建议重点检查。", "Review them.")
fails = {}  # 列名 -> 低于公开标准的样本
for (i, c), lst in flags.items():
    if any(l == "fail" for l, *_ in lst):
        fails.setdefault(base(c), []).append(df.loc[i, "sample"])
if fails:
    zh = "总体：低于公开标准的样本——" + "；".join(f"{k}：{names(v)}" for k, v in fails.items()) + "。" + FAIL_NOTE[0]
    en = "Overall: samples below published standards — " + "; ".join(f"{k}: {names(v)}" for k, v in fails.items()) + ". " + FAIL_NOTE[1]
else:
    zh = "总体：没有样本低于公开标准；" + ("标为橙色的样本建议结合实验信息复核。" if flagged else "整体质量良好。")
    en = "Overall: no sample falls below published standards; " + ("review the samples marked in orange with the experimental context." if flagged else "overall quality is good.")
points.append((zh, en))
conclusion = "".join(f"<li>{L(zh, en)}</li>" for zh, en in points)
sources = " · ".join(
    f'<a href="{esc(s["url"])}" target="_blank" rel="noopener" title="{esc(s["name"])}">{esc(k)}</a>'
    for k, s in SOURCES.items()
)

# ------------------------------------------------------------
# 流程与参数
# ------------------------------------------------------------
svg = re.sub(r"^<\?xml[^>]*>\s*", "", Path(snakemake.input.flowchart).read_text())
settings_html = "".join(f"<tr><th>{L(kz, ke)}</th><td>{L(vz, ve)}</td></tr>" for kz, ke, vz, ve in p.settings)


# ------------------------------------------------------------
# QC 汇总表
# ------------------------------------------------------------
def head_tip(c):
    m = METRICS.get(base(c), {})
    dz, de = zh_en(m.get("desc", ""))
    if m.get("ref"):
        rz, re_ = zh_en(m["ref"])
        dz += f"\n参考值：{rz}"
        de += f"\nReference: {re_}"
    if m.get("sources"):
        dz += "\n出处：" + ", ".join(m["sources"])
        de += "\nSource: " + ", ".join(m["sources"])
    return T(dz, de) if dz else ""


head = "".join(
    f'<th{head_tip(c)}{" class=num" if is_num(c) else ""}>{esc(c)}<span class="sort"></span></th>' for c in shown
)
body = []
for i, row in df.iterrows():
    cells = []
    for c in shown:
        v = row[c]
        cls = (["num"] if is_num(c) else []) + (["sample"] if c == "sample" else [])
        tip = ""
        if (i, c) in flags:
            lst = flags[(i, c)]
            cls.append("fail" if any(l == "fail" for l, *_ in lst) else "warn")
            tip = T("；".join(z for _, z, _ in lst), "; ".join(e for *_, e in lst))
        data = "" if pd.isna(v) else str(v)
        cells.append(f'<td class="{" ".join(cls)}" data-v="{esc(data)}"{tip}>{esc(fmt(c, v))}</td>')
    body.append("<tr>" + "".join(cells) + "</tr>")
table_html = f'<table class="qc" id="qc" data-key="{esc(p.pipeline)}"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>'

# ------------------------------------------------------------
# 重复一致性
# ------------------------------------------------------------
def inline_svg(path, prefix):
    """内嵌 R 输出的 SVG：去掉 XML 声明和固定尺寸；各图的 id（glyph、clip）加前缀，避免互相冲突。"""
    t = re.sub(r"^<\?xml[^>]*>\s*", "", Path(path).read_text())
    t = re.sub(r'\bid="([^"]+)"', rf'id="{prefix}-\1"', t)
    t = re.sub(r'((?:xlink:)?href)="#([^"]+)"', rf'\1="#{prefix}-\2"', t)
    t = re.sub(r"url\(#([^)]+)\)", rf"url(#{prefix}-\1)", t)
    return re.sub(r'(<svg[^>]*?) width="[^"]*" height="[^"]*"', r"\1", t, count=1)


def fig_img(fig, path, k):
    """散点矩阵的点为位图，内嵌 PNG（位图在 cairo SVG 中的显示因浏览器而异）；其它图内嵌 SVG。"""
    if fig == "correlation_scatter":
        png = base64.b64encode(Path(path).with_suffix(".png").read_bytes()).decode()
        return f'<img src="data:image/png;base64,{png}" alt="">'
    return inline_svg(path, f"f{k}")


def rep_fmt(key, v):
    if pd.isna(v):
        return ""
    if key in ("peaks", "median_peaks", "reproducible_peaks", "n_samples"):
        return f"{int(v):,}"
    if isinstance(v, float):
        return f"{v:.{REP_COLS.get(key, {}).get('digits', 2)}f}"
    return str(v)


def rep_table(frame, keys, show=None):
    """可排序的表；show[key](row) 可改写显示内容；样本行按 flags 标色。"""
    show = show or {}

    def tip(k):
        m = REP_COLS.get(k, {})
        dz, de = zh_en(m.get("desc", ""))
        if m.get("ref"):
            rz, re_ = zh_en(m["ref"])
            dz, de = dz + f"\n参考值：{rz}", de + f"\nReference: {re_}"
        return T(dz, de) if dz else ""

    num = {k: pd.api.types.is_numeric_dtype(frame[k]) for k in keys}
    head = "".join(f'<th{tip(k)}{" class=num" if num[k] else ""}>{esc(k)}<span class="sort"></span></th>' for k in keys)
    body = []
    for _, r in frame.iterrows():
        i = row_of.get(r.get("sample")) if "sample" in keys else None
        cells = []
        for k in keys:
            cls = (["num"] if num[k] else []) + (["sample"] if k in ("sample", "group") and k == keys[0] else [])
            t = ""
            if (i, k) in flags:
                lst = flags[(i, k)]
                cls.append("fail" if any(l == "fail" for l, *_ in lst) else "warn")
                t = T("；".join(z for _, z, _ in lst), "; ".join(e for *_, e in lst))
            v = r[k]
            text = show[k](r) if k in show else rep_fmt(k, v)
            cells.append(f'<td class="{" ".join(cls)}" data-v="{esc("" if pd.isna(v) else v)}"{t}>{esc(text)}</td>')
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f'<div class="scroll"><table class="qc"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


rep_html = ""
if REP:
    method = "Spearman" if rinfo["method"] == "spearman" else "Pearson"
    n_cons = f"{rinfo['consensus_peaks']:,}"
    sm = rinfo.get("summits")
    cons_zh = (f"consensus peak 由 DiffBind 得到：各组组内 consensus（至少 {rinfo['min_reps']} 个重复）取并集，"
               + (f"dba.count 以全部 IP 样本 reads 的峰顶 ± {sm} bp 计数。" if sm else "dba.count 按合并后的原宽度计数（有 broad 的组时不按峰顶重新定位）。"))
    cons_en = (f"Consensus peaks come from DiffBind: the union of the within-group consensus sets (≥ {rinfo['min_reps']} replicates), "
               + (f"counted by dba.count at the summit of the pooled IP reads ± {sm} bp." if sm
                  else "counted by dba.count over the merged widths (no summit recentring when broad groups are present)."))
    lo, hi = rinfo["y_frac"]
    mm10 = rinfo.get("genome") == "mm10"
    big_groups = [g for g in rinfo["rep_groups"] if len(rinfo["groups"][g]) > 4]
    FIGS = {
        "correlation": (
            ("样本间相关性", "Sample-to-sample correlation of ChIP signal"),
            (f"{method} 相关系数，基于 {n_cons} 个 consensus peak 上 IP 计数的 log2(CPM + 1)（不扣除 input），按 1 − r 层次聚类（average linkage），"
             f"上方色条为分组。{cons_zh}色标为 viridis 顺序色板，范围为最小 r 到 1：样本间相关系数均为正，以 0 为中心的发散色板会浪费一半色阶、难以分辨高相关区间的差异。",
             f"{method} correlation of log2(CPM + 1) of IP counts over {n_cons} consensus peaks (input not subtracted), clustered by 1 − r (average linkage); the bar on top "
             f"shows the groups. {cons_en} A sequential viridis scale runs from the lowest r to 1: correlations between samples are all positive, so a "
             "diverging scale centred on 0 would waste half of its range and hide differences among high correlations."),
        ),
        "correlation_scatter": (
            ("样本间相关性散点矩阵", "Scatter-plot matrix of ChIP signal between samples"),
            ("下三角：每个点为一个 consensus peak，坐标为两个样本 IP 计数的 log2(CPM + 1)（不扣除 input）"
             + (f"；共 {n_cons} 个 peak，随机抽取 {rinfo['scatter_peaks']:,} 个作图" if rinfo["scatter_peaks"] < rinfo["consensus_peaks"] else "")
             + f"，虚线为 y = x。上三角：{method} 相关系数（基于全部 consensus peak，与相关性热图相同），底色与热图同一色标。"
             "对角：样本名，颜色表示分组。样本按分组排列。",
             "Lower triangle: each point is one consensus peak, plotted as log2(CPM + 1) of the IP counts of two samples (input not subtracted)"
             + (f"; {rinfo['scatter_peaks']:,} of the {n_cons} peaks are drawn at random" if rinfo["scatter_peaks"] < rinfo["consensus_peaks"] else "")
             + f"; the dashed line is y = x. Upper triangle: {method} correlation over all consensus peaks (as in the correlation heatmap), "
             "shaded with the heatmap's scale. Diagonal: sample names coloured by group. Samples are ordered by group."),
        ),
        "pca": (
            ("主成分分析", "Principal component analysis of ChIP signal"),
            (f"基于方差最大的 {rinfo['pca_peaks']} 个 consensus peak 上 IP 计数的 log2(CPM + 1)（中心化，未标准化）；坐标轴括号内为解释的方差比例，点旁为重复编号。",
             f"log2(CPM + 1) of IP counts on the {rinfo['pca_peaks']} most variable consensus peaks (centred, not scaled); axis labels give the "
             "proportion of variance explained; points are labelled with the replicate."),
        ),
        "peak_overlap": (
            ("组内重复 peak 的重叠", "Pairwise overlap of replicate peaks within groups"),
            ("每个点为组内一对重复的 Jaccard 指数（交集 / 并集，按碱基，bedtools jaccard），横线为中位数。数值随 peak 数和测序深度变化，主要用于组间比较。",
             "Each point is the Jaccard index (intersection / union in bp, bedtools jaccard) of one pair of replicates; bars show the "
             "median. Values depend on peak numbers and depth, so compare across groups."),
        ),
        "peak_venn": (
            ("组内重复 peak 覆盖碱基的 Venn 图", "Venn diagrams of bases covered by replicate peaks within groups"),
            ("每个区域的数字为只被这几个重复的 peak 覆盖的碱基数（Mb，bedtools multiinter），颜色按同一数值着色；"
             "一个圆（椭圆）内各区域之和即该重复 peak 覆盖的总碱基数，与重复的顺序无关。区域面积不按比例。"
             "按碱基而非 peak 个数统计，避免一个重复中的宽 peak 对应另一重复中多个窄 peak 时计数不一致；与 Jaccard 同一口径"
             "（Jaccard = 两者共有碱基 / 并集碱基）。"
             + (f"重复数超过 4 的组（{names(big_groups)}）不画，见 Jaccard 图。" if big_groups else ""),
             "Each number is the number of bases (Mb, bedtools multiinter) covered by the peaks of exactly those replicates, and the "
             "colour encodes the same value; the regions inside a circle (ellipse) add up to the bases covered by that replicate, "
             "independent of replicate order. Areas are not to scale. Bases rather than peak numbers are counted so that one broad "
             "peak matching several narrow peaks in another replicate is counted consistently; this matches the Jaccard index "
             "(shared bases / union)."
             + (f" Groups with more than 4 replicates ({names(big_groups)}) are not drawn; see the Jaccard figure." if big_groups else "")),
        ),
        "group_venn": (
            ("组间 peak 覆盖碱基的 Venn 图", "Venn diagrams of bases covered by the peaks of different groups"),
            (f"每张小图比较几个组。每个区域的数字为只被这几个组的 reproducible peak 覆盖的碱基数（Mb，bedtools multiinter），颜色按同一数值着色；"
             f"reproducible peak 为组内至少 {rinfo['min_reps']} 个重复中出现的区域（DiffBind 组内 consensus），只有 1 个样本的组用该样本的 peak。"
             "一个圆（椭圆）内各区域之和即该组 reproducible peak 覆盖的总碱基数。区域面积不按比例。"
             "peak 集合的差别也受测序深度和重复数影响，定量比较请用差异分析（diffbind.rds）。"
             "比较哪些组由配置 replicates.group_venn 指定（为空时共 2–4 个组则全部比较）。",
             "Each panel compares a few groups. Each number is the number of bases (Mb, bedtools multiinter) covered by the reproducible "
             f"peaks of exactly those groups, and the colour encodes the same value; reproducible peaks are regions found in at least "
             f"{rinfo['min_reps']} replicates of a group (DiffBind within-group consensus), and a single-sample group uses that sample's peaks. "
             "The regions inside a circle (ellipse) add up to the bases covered by that group. Areas are not to scale. Differences "
             "between peak sets also reflect depth and replicate numbers; use differential analysis (diffbind.rds) for quantitative "
             "comparisons. The groups compared are set by replicates.group_venn in the config (when empty, all groups are compared "
             "if there are 2–4)."),
        ),
        "sex": (
            ("根据 chrX、chrY reads 占比推断性别", "Sex inferred from the fraction of reads on chrX and chrY"),
            (f"两轴均为占全部比对 reads（过滤 MAPQ ≥ {rinfo['mapq']} 并去重后）的百分比。纵轴只统计 chrY 计数区 {rinfo['y_region']}"
             "（对数轴，0 显示为 0.00001%）：即 chrY 去掉 PAR（拟常染色体区，X、Y 上序列相同，雌性 reads 也会比对上）后的部分"
             + ("，并在 PAR 起点前多去掉约 745 kb，因为雌性样本在那里也有大量 reads" if mm10 else "") + "。该区域仍含与 chrX 同源的扩增序列，"
             f"雌性样本也有少量 reads，因此按阈值判定：虚线为阈值（< {lo}% 雌性，> {hi}% 雄性，之间无法判断）。"
             "横轴 chrX 占比在 IP 样本中随抗体、细胞类型变化，不用于判定，只在组内比较；input 为全基因组 DNA，其 chrX 占比反映 X 拷贝数。"
             "包括对照样本（按其 group 着色）。只标注未判为雌性的样本。",
             f"Both axes give the percentage of all mapped reads (MAPQ ≥ {rinfo['mapq']}, deduplicated). The y axis counts only the chrY "
             f"region {rinfo['y_region']} (log scale; 0 shown as 0.00001%): chrY without the pseudoautosomal region (PAR, identical on X "
             "and Y, so female reads also map there)"
             + ("; a further ~745 kb before the PAR is excluded because female samples also have many reads there" if mm10 else "")
             + ". The region still contains X-homologous ampliconic sequence, so female samples carry a few reads "
             f"and the call uses thresholds (dashed lines: < {lo}% female, > {hi}% male, undetermined in between). The chrX fraction on "
             "the x axis varies with antibody and cell type in IP samples, so it is only compared within groups and not used for the "
             "call; in input DNA it reflects the X copy number. Control samples are included (coloured by their group). Only samples "
             "not called female are labelled."),
        ),
    }
    figs_html = []
    for k, (fig, path) in enumerate(zip(p.rep_figures, inp.rep_svg), 1):
        (tz, te), (nz, ne) = FIGS[fig]
        links = " · ".join(f'<a href="replicates/figures/{fig}.{x}" target="_blank">{x.upper()}</a>' for x in ("pdf", "png", "svg"))
        figs_html.append(
            f'<figure class="card apa" id="fig-{fig}"><div class="fig-label">{L(f"图 {k}", f"Figure {k}")}</div>'
            f'<div class="fig-title">{L(tz, te)}</div><div class="fig-img">{fig_img(fig, path, k)}</div>'
            f'<figcaption><i>{L("注：", "Note.")}</i> {L(nz, ne)}</figcaption><div class="muted dl">{links}</div></figure>'
        )

    cols_s = ["sample", "group"] + (["role"] if (rs["role"] == "control").any() else []) + [
        "peaks", "within_cor", "max_other_cor", "mean_jaccard"]
    if "sex" in rs:
        cols_s += ["y_frac", "x_frac", "sex"]
    sample_table = rep_table(rs, cols_s, {
        "max_other_cor": lambda r: f"{r.max_other_cor:.3f} ({r.max_other_sample})" if pd.notna(r.max_other_cor) else "",
    })
    group_table = rep_table(rg, ["group", "n_samples", "median_peaks", "reproducible_peaks", "median_jaccard", "min_jaccard", "mean_within_cor"])

    # 要点
    rp = []
    single = [g for g, n_ in zip(rg["group"], rg["n_samples"]) if n_ < 2]
    rp.append((
        f"分组：{len(rg)} 组，{len(rg) - len(single)} 组有 ≥ 2 个重复" + (f"（只有 1 个样本的组不参与组内比较：{names(single)}）" if single else "")
        + f"；各组 reproducible peak 取并集为 {n_cons} 个 consensus peak（DiffBind" + (f"，峰顶 ± {sm} bp" if sm else "，原宽度") + "），用于相关性和 PCA。",
        f"Groups: {len(rg)}, {len(rg) - len(single)} with ≥ 2 replicates" + (f" (single-sample groups are not compared within: {names(single)})" if single else "")
        + f"; the union of the reproducible peaks of all groups gives {n_cons} consensus peaks (DiffBind" + (f", summit ± {sm} bp" if sm else ", merged widths") + ") used for correlation and PCA.",
    ))
    if len(reps):
        wc_m, wc_r = rng(reps["within_cor"], 1, 3)
        worst = rg.loc[rg["mean_within_cor"].idxmin()]
        closer = [r["sample"] for _, r in reps.iterrows() if r.max_other_cor > r.within_cor]
        rp.append((
            f"相关性：样本与同组重复的平均 {method} r 中位数 {wc_m}（{wc_r}），最低的组为 {worst.group}（{worst.mean_within_cor:.3f}）。"
            + (f"{len(closer)} 个样本与其他组某个样本的相关性高于组内平均（{names(closer)}）；如果只发生在差异较小的组之间（如同一细胞类型的不同处理），属正常现象。" if closer else "所有样本都与本组重复最相似。"),
            f"Correlation: median within-group {method} r {wc_m} ({wc_r}); lowest group {worst.group} ({worst.mean_within_cor:.3f}). "
            + (f"{len(closer)} samples correlate more with a sample of another group than with their own replicates on average ({names(closer)}); this is expected between groups that differ little (e.g. one cell type under two conditions)." if closer else "Every sample is most similar to its own replicates."),
        ))
        jg = rg.dropna(subset=["median_jaccard"])
        lowj = jg.loc[jg["median_jaccard"].idxmin()]
        rp.append((
            f"peak 重现性：组内重复 Jaccard 中位数 {jg['median_jaccard'].min():.3f}–{jg['median_jaccard'].max():.3f}，最低的组为 {lowj.group}"
            f"（最低一对 {lowj.min_jaccard:.3f}）。Jaccard 随 peak 数和测序深度变化，深度低或 FRiP 低的样本通常也较低。",
            f"Peak reproducibility: median within-group Jaccard {jg['median_jaccard'].min():.3f}–{jg['median_jaccard'].max():.3f}; lowest group {lowj.group} "
            f"(lowest pair {lowj.min_jaccard:.3f}). Jaccard depends on peak numbers and depth, so shallow or low-FRiP samples tend to score lower.",
        ))
    if rinfo["pca_variance"]:
        v = [x * 100 for x in rinfo["pca_variance"]]
        rp.append((f"PCA：PC1 解释 {v[0]:.1f}% 的方差，PC2 解释 {v[1]:.1f}%。", f"PCA: PC1 explains {v[0]:.1f}% of the variance, PC2 {v[1]:.1f}%."))
    if "sex" in rs:
        vc = rs["sex"].value_counts()
        cnt_z = "、".join(f"{k} {v}" for k, v in vc.items())
        sex_flag = sorted({rs.loc[j, "sample"] for j in rs.index
                           if any((row_of[rs.loc[j, "sample"]], k) in flags for k in ("sex", "x_frac"))})
        rp.append((
            f"性别（按 chrY 计数区 {rinfo['y_region']} 的 reads 占比判定）：{cnt_z}。"
            + (f"需复核：{names(sex_flag)}（无法判断、与同组其他样本不同或 chrX 占比偏低；也可能是混入了另一性别的样本）。" if sex_flag else "各组内一致。"),
            f"Sex (called from the fraction of reads in the chrY region {rinfo['y_region']}): {cnt_z.replace('、', ', ')}. "
            + (f"Review: {names(sex_flag)} (undetermined, unlike the rest of the group or low chrX fraction; may also indicate a mixed sample)." if sex_flag else "Consistent within every group."),
        ))
    rep_keys = ("within_cor", "mean_jaccard", "sex", "x_frac")
    rep_flagged = sorted({df.loc[i, "sample"] for (i, c) in flags if c in rep_keys})
    rp.append((
        "总体：" + (f"建议复核 {names(rep_flagged)}（原因见上表标色单元格）。" if rep_flagged else "各组重复一致性良好。"),
        "Overall: " + (f"review {names(rep_flagged)} (reasons on the colored cells above)." if rep_flagged else "replicates are consistent within every group."),
    ))
    rep_points = "".join(f"<li>{L(zh, en)}</li>" for zh, en in rp)
    tsvs = " · ".join(f'<a href="replicates/{f}" target="_blank">{f}</a>'
                      for f in ("samples.tsv", "groups.tsv", "peak_overlap.tsv", "peak_venn.tsv", "group_venn.tsv", "correlation.tsv",
                                "log2cpm.tsv", "pca.tsv", "peak_counts.tsv"))
    rep_html = f"""
<section id="replicates">
  <h2>{L("重复一致性", "Replicate consistency")}</h2>
  <div class="card">
    <h3>{L("每个样本", "Per sample")}</h3>
    <p class="muted">{L("分组见 samples.tsv 的 group 列（为空时按样本名去掉 _repN 推断）。只比较 IP 样本；对照（input / IgG）只参与性别检查。点击表头排序；鼠标停在表头上显示含义，停在标色单元格上显示原因。",
                        "Groups come from the group column of samples.tsv (inferred by removing _repN when empty). Only IP samples are compared; controls (input / IgG) enter the sex check only. Click a header to sort; hover a header for its definition or a colored cell for the reason.")}</p>
    {sample_table}
  </div>
  <div class="card">
    <h3>{L("每组", "Per group")}</h3>
    {group_table}
  </div>
  <div class="card digest">
    <ul class="points">{rep_points}</ul>
    <p class="muted">{L("数据", "Data")}：{tsvs}</p>
  </div>
  {"".join(figs_html)}
</section>
"""


# ------------------------------------------------------------
# 信号 profile（deeptools plotProfile 的 PNG）
# ------------------------------------------------------------
def mx_opt(cmd, names, default):
    """computeMatrix 命令中某个参数的值（未写时用 deeptools 默认值）。"""
    m = re.search(r"(?:^|\s)(?:" + "|".join(names) + r")[\s=]+(\d+)", cmd)
    return int(m.group(1)) if m else default


def n_lines(path):
    with open(path) as f:
        return sum(1 for _ in f)


prof = p.profile
prof_figs = []
if "prof_tss" in inp.keys():
    up, down = mx_opt(prof["tss"], ["-b", "--upstream", "--beforeRegionStartLength"], 500), mx_opt(prof["tss"], ["-a", "--downstream", "--afterRegionStartLength"], 500)
    n_tss = n_lines(inp.prof_tss_bed)
    prof_figs.append((
        "tss_profile", inp.prof_tss,
        ("TSS profile", "TSS profile"),
        (f"每条线为一个样本（包括对照）在 {n_tss:,} 个 protein_coding 基因 TSS 上游 {up:,} bp 到下游 {down:,} bp 的平均信号"
         f"（bigwig，{prof['norm']} 归一化）。基因按 GTF 的链方向排列，左为上游。启动子标记和开放染色质在 TSS 处形成尖峰，"
         "对照（input / IgG）应接近平坦；汇总表中的 TSS enrichment 由此图的数据计算。",
         f"Each line is the mean signal of one sample (controls included) from {up:,} bp upstream to {down:,} bp downstream of the "
         f"TSS of {n_tss:,} protein_coding genes (bigwig, {prof['norm']} normalized). Genes are oriented by their GTF strand, "
         "upstream on the left. Promoter marks and open chromatin peak sharply at the TSS, while controls (input / IgG) should "
         "stay flat; the TSS enrichment in the summary table is computed from these data."),
    ))
if "prof_samples" in inp.keys():
    mg = prof["metagene"]
    b, a = mx_opt(mg, ["-b", "--upstream", "--beforeRegionStartLength"], 0), mx_opt(mg, ["-a", "--downstream", "--afterRegionStartLength"], 0)
    body = mx_opt(mg, ["-m", "--regionBodyLength"], 1000)
    u5, u3 = mx_opt(mg, ["--unscaled5prime"], 0), mx_opt(mg, ["--unscaled3prime"], 0)
    n_genes = n_lines(inp.prof_genes)
    unscaled_zh = (f"TSS 下游 {u5:,} bp" if u5 else "") + ("、" if u5 and u3 else "") + (f"TES 上游 {u3:,} bp" if u3 else "")
    unscaled_en = (f"the first {u5:,} bp after the TSS" if u5 else "") + (" and " if u5 and u3 else "") + (f"the last {u3:,} bp before the TES" if u3 else "")
    scale_zh = (f"TSS 上游 {b:,} bp、TES 下游 {a:,} bp" + (f"以及{unscaled_zh}" if unscaled_zh else "") + "按真实距离；"
                f"其余基因体统一缩放到 {body:,} bp，这一段横轴是基因内的相对位置，不是真实距离。")
    scale_en = (f"{b:,} bp upstream of the TSS and {a:,} bp downstream of the TES" + (f", as well as {unscaled_en}," if unscaled_en else "")
                + f" are drawn to true distance; the rest of each gene body is scaled to {body:,} bp, so that part of the x axis is the "
                "relative position within the gene, not distance.")
    genes_zh = f"{n_genes:,} 个 protein_coding 基因（TSS 到 TES，长度 ≥ {prof['min_len']:,} bp）"
    genes_en = f"{n_genes:,} protein_coding genes (TSS to TES, ≥ {prof['min_len']:,} bp)"
    prof_figs.append((
        "metagene_samples", inp.prof_samples,
        ("Meta gene profile：每个样本", "Meta gene profile: per sample"),
        (f"每条线为一个样本（包括对照）在 {genes_zh}上的平均信号（bigwig，{prof['norm']} 归一化）。{scale_zh}"
         "基因按 GTF 的链方向排列，左为 5′ 端。对照应接近平坦；启动子标记集中在 TSS，基因体标记（如 H3K36me3、H3K79me2）"
         "在基因体内升高。",
         f"Each line is the mean signal of one sample (controls included) over {genes_en} (bigwig, {prof['norm']} normalized). "
         f"{scale_en} Genes are oriented by their GTF strand, 5′ end on the left. Controls should stay flat; promoter marks "
         "concentrate at the TSS, gene-body marks (e.g. H3K36me3, H3K79me2) rise across the gene body."),
    ))
    if "prof_groups" in inp.keys():
        gz, ge = (("合并组内重复的 BAM 后生成", "replicate BAMs merged") if prof["group_method"] == "merge"
                  else ("各重复的 bigwig 逐 bin 取平均", "replicate bigwigs averaged per bin"))
        prof_figs.append((
            "metagene_groups", inp.prof_groups,
            ("Meta gene profile：每组", "Meta gene profile: per group"),
            (f"同上图，每条线为一组的组 bigwig（{gz}；只包括 IP 样本、有 ≥ 2 个参与合并的重复的组）。",
             f"As above, each line is the bigwig of one group ({ge}; IP samples only, groups with ≥ 2 merged replicates)."),
        ))

prof_html = ""
if prof_figs:
    k0 = len(p.rep_figures) if REP else 0
    items = []
    for k, (name, path, (tz, te), (nz, ne)) in enumerate(prof_figs, k0 + 1):
        png = base64.b64encode(Path(path).read_bytes()).decode()
        stem = Path(path).stem
        links = " · ".join(f'<a href="deeptools/{stem}.{x}" target="_blank">{x.upper()}</a>' for x in ("png", "tsv"))
        items.append(
            f'<figure class="card apa" id="fig-{name}"><div class="fig-label">{L(f"图 {k}", f"Figure {k}")}</div>'
            f'<div class="fig-title">{L(tz, te)}</div><div class="fig-img prof"><img src="data:image/png;base64,{png}" alt=""></div>'
            f'<figcaption><i>{L("注：", "Note.")}</i> {L(nz, ne)}</figcaption><div class="muted dl">{links}</div></figure>'
        )
    prof_html = f"""
<section id="profiles">
  <h2>{L("信号 profile", "Signal profiles")}</h2>
  {"".join(items)}
</section>
"""


# ------------------------------------------------------------
# 页面
# ------------------------------------------------------------
title_zh, title_en = p.title
# 项目：输出目录（snakemake --directory，即当前目录）的上两级，如 Patri_aging/out_atac
OUTDIR = Path.cwd()
PROJECT = "/".join(OUTDIR.parts[-2:])
CSS = """
:root {
  --bg: #f5f5f7; --card: #ffffff; --text: #1d1d1f; --muted: #6e6e73; --line: #d2d2d7; --soft: #f5f5f7;
  --accent: #0071e3; --green: #248a3d; --orange: #c93400; --red: #d70015;
  --orange-bg: #fff4e5; --red-bg: #ffebee; --shadow: 0 2px 12px rgba(0,0,0,.06); --nav: rgba(255,255,255,.72);
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #000000; --card: #1c1c1e; --text: #f5f5f7; --muted: #98989d; --line: #38383a; --soft: #2c2c2e;
    --accent: #2997ff; --green: #30d158; --orange: #ff9f0a; --red: #ff453a;
    --orange-bg: rgba(255,159,10,.16); --red-bg: rgba(255,69,58,.18); --shadow: none; --nav: rgba(28,28,30,.72);
  }
}
* { box-sizing: border-box; }
html { scroll-behavior: smooth; scroll-padding-top: 24px; }
body {
  margin: 0; background: var(--bg); color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, "SF Pro Text", "PingFang SC", "Helvetica Neue", "Noto Sans CJK SC", "Microsoft YaHei", sans-serif;
  font-size: 14px; line-height: 1.55; -webkit-font-smoothing: antialiased;
}
body.en [data-l="zh"], body:not(.en) [data-l="en"] { display: none; }
a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }
code { font-family: "SF Mono", ui-monospace, Menlo, Consolas, monospace; font-size: 12.5px; background: var(--soft); padding: 1px 6px; border-radius: 6px; word-break: break-all; }
nav {
  position: fixed; inset: 0 auto 0 0; width: 232px; padding: 28px 16px; background: var(--nav);
  backdrop-filter: saturate(180%) blur(20px); -webkit-backdrop-filter: saturate(180%) blur(20px); border-right: 1px solid var(--line);
}
.brand { font-size: 17px; font-weight: 600; letter-spacing: -.01em; padding: 0 12px; }
.brand-sub { color: var(--muted); font-size: 12px; padding: 2px 12px 16px; }
.seg { display: inline-flex; background: var(--soft); border-radius: 9px; padding: 2px; margin: 0 12px 16px; }
.seg button { border: 0; background: none; color: var(--text); font: inherit; font-size: 12.5px; padding: 4px 14px; border-radius: 7px; cursor: pointer; }
.seg button.on { background: var(--card); box-shadow: 0 1px 3px rgba(0,0,0,.12); font-weight: 600; }
nav a.item { display: block; padding: 7px 12px; margin: 2px 0; border-radius: 8px; color: var(--text); }
nav a.item:hover { background: var(--soft); text-decoration: none; }
nav a.item.active { background: var(--accent); color: #fff; }
nav .ext { margin-top: 18px; padding-top: 14px; border-top: 1px solid var(--line); }
nav .ext a { display: block; padding: 6px 12px; font-size: 13px; }
main { margin-left: 232px; padding: 40px 48px 80px; max-width: 1480px; }
header.top h1 { font-size: 32px; font-weight: 700; letter-spacing: -.02em; margin: 0; }
header.top .project { margin-top: 8px; font-size: 15px; font-weight: 600; }
header.top .project .path { color: var(--muted); font-weight: 400; font-size: 13px; margin-left: 8px; overflow-wrap: anywhere; }
header.top .meta { color: var(--muted); margin-top: 6px; font-size: 13px; }
section { margin-top: 44px; }
section > h2 { font-size: 24px; font-weight: 600; letter-spacing: -.015em; margin: 0 0 16px; }
h3 { font-size: 17px; font-weight: 600; margin: 0 0 12px; }
.card { background: var(--card); border-radius: 18px; padding: 24px; box-shadow: var(--shadow); margin-bottom: 20px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr)); gap: 14px; margin-bottom: 20px; }
.tile { background: var(--card); border-radius: 16px; padding: 16px 18px; box-shadow: var(--shadow); }
.tile-label { color: var(--muted); font-size: 12.5px; }
.tile-value { font-size: 28px; font-weight: 600; letter-spacing: -.02em; margin-top: 4px; font-variant-numeric: tabular-nums; }
.tile-note { color: var(--muted); font-size: 11.5px; min-height: 1.2em; }
.tile-alert .tile-value { color: var(--orange); }
.chip { display: inline-block; font-size: 12.5px; padding: 3px 10px; border-radius: 999px; margin: 2px 6px 2px 0; }
.chip-warn { background: var(--orange-bg); color: var(--orange); }
.chip-fail { background: var(--red-bg); color: var(--red); }
.ok-note { color: var(--green); margin: 0; }
table { border-collapse: collapse; width: 100%; }
table.plain td { padding: 8px 0; border-bottom: 1px solid var(--line); vertical-align: top; }
table.plain tr:last-child td { border-bottom: 0; }
table.plain td.sample { width: 220px; font-weight: 500; padding-right: 16px; }
table.kv th { text-align: left; font-weight: 500; color: var(--muted); width: 200px; padding: 10px 16px 10px 0; vertical-align: top; }
table.kv td { padding: 10px 0; }
table.kv tr + tr th, table.kv tr + tr td { border-top: 1px solid var(--line); }
.flow { background: #fff; border-radius: 12px; padding: 8px; max-width: 760px; margin: 0 auto; }
.flow svg { width: 100%; height: auto; display: block; }
.muted { color: var(--muted); font-size: 12.5px; }
.legend span.sw { display: inline-block; width: 10px; height: 10px; border-radius: 3px; margin: 0 4px 0 12px; vertical-align: -1px; }
.bar { display: flex; align-items: center; justify-content: space-between; gap: 12px; margin-bottom: 12px; }
.btn { border: 0; background: var(--accent); color: #fff; font: inherit; font-size: 13px; padding: 6px 14px; border-radius: 999px; cursor: pointer; white-space: nowrap; }
.btn:hover { filter: brightness(1.08); }
.tools { display: flex; gap: 8px; position: relative; }
.btn.ghost { background: var(--soft); color: var(--text); border: 1px solid var(--line); }
.colpick { display: none; position: absolute; right: 0; top: calc(100% + 6px); z-index: 10; background: var(--card);
  border: 1px solid var(--line); border-radius: 12px; box-shadow: var(--shadow); padding: 10px 14px;
  max-height: 60vh; overflow: auto; min-width: 220px; font-size: 13px; }
.colpick.open { display: block; }
.colpick label { display: block; white-space: nowrap; padding: 2px 0; cursor: pointer; }
.colpick .all { color: var(--accent); cursor: pointer; margin-bottom: 6px; display: inline-block; }
.scroll { overflow: auto; max-height: 78vh; border-radius: 12px; border: 1px solid var(--line); }
table.qc { font-size: 13px; font-variant-numeric: tabular-nums; }
table.qc th {
  position: sticky; top: 0; z-index: 2; background: var(--soft); text-align: left; font-weight: 600; white-space: nowrap;
  padding: 10px 12px; border-bottom: 1px solid var(--line); cursor: pointer; user-select: none;
}
table.qc th:hover { color: var(--accent); }
table.qc td { padding: 8px 12px; border-bottom: 1px solid var(--line); white-space: nowrap; background: var(--card); }
table.qc tr:last-child td { border-bottom: 0; }
table.qc tbody tr:hover td { background: var(--soft); }
table.qc .num { text-align: right; }
table.qc th:first-child, table.qc td:first-child { position: sticky; left: 0; z-index: 1; }
table.qc th:first-child { z-index: 3; }
table.qc td.sample { font-weight: 500; }
table.qc td.warn { background: var(--orange-bg); color: var(--orange); font-weight: 600; }
table.qc td.fail { background: var(--red-bg); color: var(--red); font-weight: 600; }
.sort::after { content: ""; margin-left: 4px; color: var(--accent); }
th[data-dir="asc"] .sort::after { content: "▲"; font-size: 9px; }
th[data-dir="desc"] .sort::after { content: "▼"; font-size: 9px; }
ul.points { margin: 0; padding-left: 20px; }
ul.points li { margin: 8px 0; }
ul.points li:last-child { font-weight: 600; }
.digest { border-left: 4px solid var(--accent); }
figure.apa { margin: 0 0 20px; }
.fig-label { font-weight: 700; }
.fig-title { font-style: italic; margin-bottom: 12px; }
.fig-img { background: #fff; border-radius: 12px; padding: 8px; max-width: 880px; }
.fig-img svg, .fig-img img { width: 100%; height: auto; display: block; }
.fig-img.prof img { max-width: 820px; margin: 0 auto; background: #fff; border-radius: 6px; }
figure.apa figcaption { margin-top: 12px; font-size: 13px; max-width: 880px; }
.dl { margin-top: 6px; }
@media (max-width: 860px) {
  nav { position: static; width: auto; border-right: 0; border-bottom: 1px solid var(--line); padding: 16px; }
  nav a.item { display: inline-block; }
  nav .ext { display: none; }
  main { margin-left: 0; padding: 24px 16px 60px; }
  table.kv th { width: 120px; }
}
"""
JS = """
// 中英文切换（记住上次的选择）
function setLang(l) {
  document.body.classList.toggle("en", l === "en");
  document.documentElement.lang = l === "en" ? "en" : "zh-CN";
  document.querySelectorAll("[data-tz]").forEach(e => e.title = l === "en" ? e.dataset.te : e.dataset.tz);
  document.querySelectorAll(".seg button").forEach(b => b.classList.toggle("on", b.dataset.lang === l));
  try { localStorage.setItem("qc-lang", l); } catch (e) {}
}
document.querySelectorAll(".seg button").forEach(b => b.addEventListener("click", () => setLang(b.dataset.lang)));
let saved = null;
try { saved = localStorage.getItem("qc-lang"); } catch (e) {}
setLang(saved || "zh");

// 表头排序（按 data-v，空值始终在最后）
document.querySelectorAll("table.qc th").forEach((th, idx) => th.addEventListener("click", () => {
  const table = th.closest("table"), tbody = table.tBodies[0];
  const dir = th.dataset.dir === "asc" ? "desc" : "asc";
  table.querySelectorAll("th").forEach(h => delete h.dataset.dir);
  th.dataset.dir = dir;
  const val = r => { const v = r.cells[idx].dataset.v; const x = Number(v); return v === "" ? null : (isNaN(x) ? v : x); };
  [...tbody.rows].sort((a, b) => {
    const x = val(a), y = val(b);
    if (x === null) return 1; if (y === null) return -1;
    return (x > y ? 1 : x < y ? -1 : 0) * (dir === "asc" ? 1 : -1);
  }).forEach(r => tbody.appendChild(r));
}));

// 显示 / 隐藏列（sample 列始终显示；选择按流程记在浏览器中）
(() => {
  const t = document.getElementById("qc"), box = document.getElementById("colpick");
  const key = "qc-cols-" + t.dataset.key;
  const names = [...t.tHead.rows[0].cells].map(c => c.firstChild.textContent);
  let hidden = new Set();
  try { hidden = new Set(JSON.parse(localStorage.getItem(key) || "[]")); } catch (e) {}
  const apply = () => {
    [...t.rows].forEach(r => [...r.cells].forEach((c, i) => { c.hidden = i > 0 && hidden.has(names[i]); }));
    box.querySelectorAll("input").forEach(x => { x.checked = !hidden.has(x.value); });
    try { localStorage.setItem(key, JSON.stringify([...hidden])); } catch (e) {}
  };
  names.slice(1).forEach(n => {
    const lab = document.createElement("label"), x = document.createElement("input");
    x.type = "checkbox"; x.value = n;
    x.addEventListener("change", () => { x.checked ? hidden.delete(n) : hidden.add(n); apply(); });
    lab.append(x, " " + n); box.appendChild(lab);
  });
  box.querySelector(".all").addEventListener("click", () => { hidden.clear(); apply(); });
  document.getElementById("cols").addEventListener("click", e => { e.stopPropagation(); box.classList.toggle("open"); });
  box.addEventListener("click", e => e.stopPropagation());
  document.addEventListener("click", () => box.classList.remove("open"));
  apply();
})();

// 复制表格：制表符分隔、纯数字，可直接粘贴到 Excel（只复制显示的列）
document.getElementById("copy").addEventListener("click", async e => {
  const t = document.getElementById("qc");
  const rows = [[...t.tHead.rows[0].cells].filter(c => !c.hidden).map(c => c.textContent)];
  [...t.tBodies[0].rows].forEach(r => rows.push([...r.cells].filter(c => !c.hidden).map(c => c.dataset.v)));
  const text = rows.map(r => r.join("\\t")).join("\\n");
  try { await navigator.clipboard.writeText(text); }
  catch (err) {
    const ta = document.createElement("textarea"); ta.value = text; document.body.appendChild(ta);
    ta.select(); document.execCommand("copy"); ta.remove();
  }
  const b = e.currentTarget; b.classList.add("done");
  setTimeout(() => b.classList.remove("done"), 1500);
});

// 左侧导航随滚动高亮
const links = [...document.querySelectorAll("nav a.item")];
const obs = new IntersectionObserver(es => es.forEach(e => {
  if (e.isIntersecting) links.forEach(a => a.classList.toggle("active", a.getAttribute("href") === "#" + e.target.id));
}), { rootMargin: "-30% 0px -60% 0px" });
document.querySelectorAll("section").forEach(s => obs.observe(s));
links[0].classList.add("active");
"""
COPY_CSS = ".btn .ok { display: none; } .btn.done .ok { display: inline; } .btn.done .go { display: none; }"

page = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title_zh)}</title>
<style>{CSS}{COPY_CSS}</style>
</head>
<body>
<nav>
  <div class="brand">{L(title_zh, title_en)}</div>
  <div class="brand-sub">{esc(p.pipeline)}</div>
  <div class="seg"><button data-lang="zh">中文</button><button data-lang="en">EN</button></div>
  <a class="item" href="#overview">{L("概览", "Overview")}</a>
  <a class="item" href="#workflow">{L("流程与参数", "Workflow & parameters")}</a>
  <a class="item" href="#summary">{L("QC 汇总表", "QC summary")}</a>
  {f'<a class="item" href="#replicates">{L("重复一致性", "Replicate consistency")}</a>' if REP else ""}
  {f'<a class="item" href="#profiles">{L("信号 profile", "Signal profiles")}</a>' if prof_html else ""}
  <div class="ext"><a href="{esc(p.multiqc)}" target="_blank">{L("MultiQC 报告", "MultiQC report")} ↗</a></div>
</nav>
<main>
<header class="top">
  <h1>{L(title_zh, title_en)}</h1>
  <div class="project">{esc(PROJECT)}<span class="path">{esc(str(OUTDIR))}</span></div>
  <div class="meta">{esc(p.pipeline)} {esc(p.version)} · {L("生成于", "generated")} {datetime.now():%Y-%m-%d %H:%M} · {L(f"{n} 个样本", f"{n} samples")}</div>
</header>

<section id="overview">
  <h2>{L("概览", "Overview")}</h2>
  <div class="tiles">{tiles_html}</div>
  <div class="card">
    <h3>{L("需关注的样本", "Samples to review")}</h3>
    {attention}
  </div>
</section>

<section id="workflow">
  <h2>{L("流程与参数", "Workflow & parameters")}</h2>
  <div class="card"><div class="flow">{svg}</div></div>
  <div class="card">
    <h3>{L("本次参数", "Parameters of this run")}</h3>
    <table class="kv"><tbody>{settings_html}</tbody></table>
    <p class="muted">{L("完整配置见 `logs/run_<时间>/config.yaml`。", "Full configuration: `logs/run_<time>/config.yaml`.")}</p>
  </div>
</section>

<section id="summary">
  <h2>{L("QC 汇总表", "QC summary")}</h2>
  <div class="card">
    <div class="bar">
      <div class="muted legend">{L("reads 按单条 read 计（PE 为 R1 + R2），rate 以进入该步骤的 reads 为分母（%，见表头括号）。点击表头排序；鼠标停在表头上显示含义与参考值，停在标色单元格上显示原因。",
                                    "Reads are single reads (R1 + R2 for PE); rates are % of the reads entering each step (see header). Click a header to sort; hover a header for its definition and reference, or a colored cell for the reason.")}
        <span class="sw" style="background:var(--orange-bg);border:1px solid var(--orange)"></span>{L("需注意", "review")}
        <span class="sw" style="background:var(--red-bg);border:1px solid var(--red)"></span>{L("低于公开标准", "below standard")}</div>
      <div class="tools">
        <button class="btn ghost" id="cols">{L("显示列", "Columns")} ▾</button>
        <div class="colpick" id="colpick"><span class="all">{L("全部显示", "Show all")}</span></div>
        <button class="btn" id="copy"><span class="go">{L("复制表格", "Copy table")}</span><span class="ok">{L("已复制", "Copied")} ✓</span></button>
      </div>
    </div>
    <div class="scroll">{table_html}</div>
  </div>
  <div class="card digest">
    <ul class="points">{conclusion}</ul>
    <p class="muted">{L("参考来源", "References")}：{sources}</p>
  </div>
</section>
{rep_html}
{prof_html}
</main>
<script>{JS}</script>
</body>
</html>
"""
Path(snakemake.output[0]).write_text(page)
print(f"{n} 个样本，{len(flagged)} 个需关注", file=sys.stderr)
