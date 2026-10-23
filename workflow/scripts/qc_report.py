"""
单文件 HTML QC 报告：概览、流程与参数、QC 汇总表（中英文可切换）。
所有样式、脚本和流程图都内嵌，不依赖网络。各列说明、参考值和检查规则来自 workflow/report/qc_metrics.yaml：
  checks  有公开阈值的检查（warn 橙 / fail 红）
  outlier 与本批样本比较，|值 - 中位数| > 3 x MAD（x 1.4826）标为橙色
"""
import html
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
table_html = f'<table class="qc" id="qc"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>'

# ------------------------------------------------------------
# 页面
# ------------------------------------------------------------
title_zh, title_en = p.title
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

// 复制表格：制表符分隔、纯数字，可直接粘贴到 Excel
document.getElementById("copy").addEventListener("click", async e => {
  const t = document.getElementById("qc");
  const rows = [[...t.tHead.rows[0].cells].map(c => c.textContent)];
  [...t.tBodies[0].rows].forEach(r => rows.push([...r.cells].map(c => c.dataset.v)));
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
  <div class="ext"><a href="{esc(p.multiqc)}" target="_blank">{L("MultiQC 报告", "MultiQC report")} ↗</a></div>
</nav>
<main>
<header class="top">
  <h1>{L(title_zh, title_en)}</h1>
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
      <button class="btn" id="copy"><span class="go">{L("复制表格", "Copy table")}</span><span class="ok">{L("已复制", "Copied")} ✓</span></button>
    </div>
    <div class="scroll">{table_html}</div>
  </div>
  <div class="card digest">
    <ul class="points">{conclusion}</ul>
    <p class="muted">{L("参考来源", "References")}：{sources}</p>
  </div>
</section>
</main>
<script>{JS}</script>
</body>
</html>
"""
Path(snakemake.output[0]).write_text(page)
print(f"{n} 个样本，{len(flagged)} 个需关注", file=sys.stderr)
