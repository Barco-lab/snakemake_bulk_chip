import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


# ============================================================
# 环境与参考
# ============================================================
ENV_MAIN = config["envs"]["main"]
ENV_MACS = config["envs"]["macs3"]

if config["genome"] not in config["references"]:
    raise ValueError(
        f"genome '{config['genome']}' 不在 references 中，可选: {list(config['references'])}"
    )
REF = config["references"][config["genome"]]

CONTAM = config["contamination"]
if CONTAM["enabled"]:
    _idx = config["references"].get(CONTAM["genome"], {}).get("bowtie2_index", "")
    if not _idx:
        raise ValueError(
            f"contamination.enabled=true，但 references.{CONTAM['genome']}.bowtie2_index 未填写"
        )
    CONTAM_INDEX = _idx

THREADS = config["threads"]


def threads_for(step):
    return THREADS.get(step, THREADS["default"])


REMOVE_MITO = config["filter"]["remove_mito"]
MITO = config["filter"]["mito_name"]

# UCSC (chr1) <-> Ensembl (1) 染色体命名转换
RENAME_CHROMS = str(Path(workflow.basedir) / "scripts/rename_chroms.awk")


# ============================================================
# 样本表
# ============================================================
REP_SUFFIX = re.compile(r"[._-]?rep\d+$", re.I)  # 与 make_samplesheet.py 一致


def infer_group(sample):
    """group 列为空时按样本名推断：去掉末尾的 _repN，如 Naive_rep1 -> Naive。"""
    return REP_SUFFIX.sub("", sample) or sample


def load_samples(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"样本表不存在: {path.resolve()}")
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#", index_col=False).fillna("")
    df.columns = df.columns.str.strip()
    for col in ("sample", "fq1"):
        if col not in df.columns:
            raise ValueError(f"样本表缺少必需列 '{col}': {path}")
    for col in ("fq2", "control", "peak_type", "group", "merge"):
        if col not in df.columns:
            df[col] = ""
    df = df.apply(lambda c: c.str.strip())
    df["peak_type"] = df["peak_type"].replace("", config["peak_type"])
    df["group"] = [g or infer_group(s) for s, g in zip(df["sample"], df["group"])]
    df["merge"] = df["merge"].str.lower()

    errors = []
    dup = df.loc[df["sample"].duplicated(), "sample"].tolist()
    if dup:
        errors.append(f"样本名重复: {dup}")
    for s in df["sample"]:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", s):
            errors.append(f"样本名只能包含字母、数字、. _ -（且不能以符号开头）: '{s}'")
    for g in df["group"].unique():
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", g):
            errors.append(f"group 只能包含字母、数字、. _ -（且不能以符号开头）: '{g}'")
    # 同一组的 IP 样本 peak 类型须一致（reproducible peak 按组合并）
    ip = df[~df["sample"].isin(set(df["control"]) - {""})]
    for g, sub in ip.groupby("group", sort=False):
        if sub["peak_type"].nunique() > 1:
            errors.append(f"group '{g}' 中的 IP 样本 peak_type 不一致: {dict(zip(sub['sample'], sub['peak_type']))}")
    for _, r in df.iterrows():
        if r["merge"] not in ("", "yes", "no"):
            errors.append(f"{r['sample']}: merge='{r['merge']}' 无效，可选 yes / no（留空即 yes）")
        if r["peak_type"] not in ("narrow", "broad"):
            errors.append(f"{r['sample']}: peak_type='{r['peak_type']}' 无效，可选 narrow / broad")
        if r["control"]:
            if r["control"] not in set(df["sample"]):
                errors.append(f"{r['sample']}: control '{r['control']}' 不在样本表中")
            elif r["control"] == r["sample"]:
                errors.append(f"{r['sample']}: control 不能是自己")

    for col in ("fq1", "fq2"):
        df[col] = [
            str((path.parent / p).resolve()) if p and not Path(p).is_absolute() else p
            for p in df[col]
        ]
    for _, r in df.iterrows():
        for col in ("fq1", "fq2"):
            if r[col] and not Path(r[col]).exists():
                errors.append(f"{r['sample']}: {col} 文件不存在: {r[col]}")

    if errors:
        raise ValueError("样本表检查未通过:\n  " + "\n  ".join(errors))
    return df.set_index("sample", drop=False)


SAMPLES_DF = load_samples(config["samples"])
SAMPLES = SAMPLES_DF["sample"].tolist()
PE_SAMPLES = [s for s in SAMPLES if SAMPLES_DF.loc[s, "fq2"]]
SE_SAMPLES = [s for s in SAMPLES if not SAMPLES_DF.loc[s, "fq2"]]
# 被别的样本当作 control 的样本只做比对和 QC，不单独 call peak
CONTROLS = sorted(set(SAMPLES_DF["control"]) - {""})
PEAK_SAMPLES = [s for s in SAMPLES if s not in CONTROLS]
# 分组（按样本表中首次出现的顺序），只含 IP 样本；REP_GROUPS 为有 ≥ 2 个重复的组。control 的组只用于性别检查
GROUPS = {
    g: [s for s in PEAK_SAMPLES if SAMPLES_DF.loc[s, "group"] == g]
    for g in SAMPLES_DF.loc[PEAK_SAMPLES, "group"].unique()
}
REP_GROUPS = [g for g, ss in GROUPS.items() if len(ss) >= 2]
REPLICATES = config["replicates"]

# 组 bigwig：组内 merge 列不为 no 的 IP 样本，≥ 2 个时才生成
BW_GROUP = config["bigwig"]["group"]
if BW_GROUP["method"] not in ("merge", "mean"):
    raise ValueError(f"bigwig.group.method='{BW_GROUP['method']}' 无效，可选 merge / mean")
BW_GROUPS = {g: [s for s in ss if SAMPLES_DF.loc[s, "merge"] != "no"] for g, ss in GROUPS.items()}
BW_GROUPS = {g: ss for g, ss in BW_GROUPS.items() if len(ss) >= 2}


def regex_of(names):
    """把样本名列表转为 wildcard 约束；空列表返回永不匹配的正则。"""
    return "|".join(re.escape(n) for n in names) if names else "(?!x)x"


wildcard_constraints:
    sample=regex_of(SAMPLES),
    group=regex_of(list(GROUPS)),
    read="R1|R2",
    ptype="narrow|broad",


def is_pe(sample):
    return sample in PE_SAMPLES


def reads_of(sample):
    return ["R1", "R2"] if is_pe(sample) else ["R1"]


def control_of(sample):
    return SAMPLES_DF.loc[sample, "control"]


def peak_file(sample):
    return f"results/03_peaks/{sample}.{SAMPLES_DF.loc[sample, 'peak_type']}Peak"


def trimmed_reads(wc):
    return [f"results/01_trim/{wc.sample}_{r}.trimmed.fq.gz" for r in reads_of(wc.sample)]


def maybe_temp(path, key):
    return path if config["keep"][key] else temp(path)


# ============================================================
# 目标文件
# ============================================================
def all_targets():
    t = []
    for s in SAMPLES:
        t += expand("results/qc/fastqc_raw/{s}_{r}_fastqc.zip", s=s, r=reads_of(s))
        t += expand("results/qc/fastqc_trimmed/{s}_{r}.trimmed_fastqc.zip", s=s, r=reads_of(s))
        t.append(f"results/04_bigwig/{s}.bw")
        if CONTAM["enabled"]:
            t.append(f"results/qc/contamination/{s}.log")
    for s in PEAK_SAMPLES:
        t.append(peak_file(s))
        if config["motif"]["enabled"] and SAMPLES_DF.loc[s, "peak_type"] == "narrow":
            t.append(f"results/05_motif/{s}/knownResults.txt")
    if BW_GROUP["enabled"]:
        t += expand("results/04_bigwig/groups/{group}.bw", group=BW_GROUPS)
    if replicates_on():
        t += expand("results/03_peaks/reproducible/{group}.bed", group=REP_GROUPS)
        t += replicate_figure_files()
    t += [
        "results/qc/summary.tsv",
        "results/qc/multiqc/multiqc_report.html",
        "results/qc/report.html",
    ]
    return t


# ============================================================
# 运行记录：logs/run_<时间>/
# ============================================================
def pipeline_version(pipe):
    """流程当前的 commit；直接读 .git，运行时不需要 git。"""
    gd = pipe / ".git"
    try:
        head = (gd / "HEAD").read_text().strip()
        if not head.startswith("ref: "):
            return head
        ref = head[5:]
        if (gd / ref).exists():
            return f"{(gd / ref).read_text().strip()} ({ref})"
        for line in (gd / "packed-refs").read_text().splitlines():
            if line.endswith(" " + ref):
                return f"{line.split()[0]} ({ref})"
    except OSError:
        pass
    return "unknown"


def record_run(version_args):
    """onstart 调用：保存本次运行的配置、样本表、流程版本与软件版本。"""
    d = Path(RUN_DIR)
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "config.yaml", "w") as f:
        yaml.safe_dump(dict(config), f, allow_unicode=True, sort_keys=False)
    shutil.copy(config["samples"], d / "samples.tsv")
    pipe = Path(workflow.basedir).parent
    (d / "pipeline.txt").write_text(
        f"pipeline: {pipe}\nversion:  {pipeline_version(pipe)}\ncommand:  {' '.join(sys.argv)}\n"
    )
    with open(d / "software_versions.tsv", "w") as f:
        subprocess.run(
            ["bash", str(pipe / "workflow/scripts/software_versions.sh"), *version_args],
            stdout=f, stderr=subprocess.DEVNULL,
            env={**os.environ, "SMK_PATH": sys.argv[0]},
        )


def finish_run(smk_log):
    """onsuccess/onerror 调用：合并本次各任务的 benchmark 为一张表，并保存 snakemake 主日志。"""
    d = Path(RUN_DIR)
    if not d.is_dir():  # dry-run
        return
    bdir = d / "benchmarks"
    tables = []
    for f in sorted(bdir.rglob("*.tsv")):
        df = pd.read_csv(f, sep="\t")
        df.insert(0, "job", f.stem)
        df.insert(0, "step", str(f.parent.relative_to(bdir)))
        tables.append(df)
    if tables:
        pd.concat(tables).to_csv(d / "benchmarks.tsv", sep="\t", index=False)
    shutil.rmtree(bdir, ignore_errors=True)
    if smk_log:
        shutil.copy(smk_log if isinstance(smk_log, str) else smk_log[0], d / "snakemake.log")
