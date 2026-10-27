# 汇总表、MultiQC、QC 报告


def qc_files():
    """汇总 QC：samtools stats（含片段长度分布）、fingerprint；TSS profile 由 config 控制。"""
    d = {
        "stats": expand("results/02_align/stats/{sample}.dedup.stats", sample=SAMPLES),
        "fingerprint": "results/qc/deeptools/fingerprint.raw_counts.tsv",
    }
    if config["qc"]["tss_profile"]:
        d["tss"] = "results/qc/deeptools/tss_profile.tsv"
    return d


def profile_files():
    """QC 报告「信号 profile」一节：deeptools plotProfile 的 PNG 及基因数（bed），由 config 控制。"""
    d = {}
    if config["qc"]["tss_profile"]:
        d["prof_tss"] = "results/qc/deeptools/tss_profile.png"
        d["prof_tss_bed"] = "results/reference/tss.bed"
    if config["qc"]["metagene_profile"]:
        d["prof_samples"] = "results/qc/deeptools/metagene_samples.png"
        d["prof_genes"] = "results/reference/genes.bed"
        if BW_GROUP["enabled"] and BW_GROUPS:
            d["prof_groups"] = "results/qc/deeptools/metagene_groups.png"
    return d


def trim_report(s, r):
    return f"results/01_trim/reports/{s}_{r}.fastq.gz_trimming_report.txt"


rule summary:
    input:
        unpack(lambda wc: qc_files()),
        trim=[trim_report(s, r) for s in SAMPLES for r in reads_of(s)],
        bowtie2=expand("results/qc/bowtie2/{sample}.log", sample=SAMPLES),
        raw=expand("results/02_align/stats/{sample}.raw.flagstat", sample=SAMPLES),
        idxstats=expand("results/02_align/stats/{sample}.raw.idxstats", sample=SAMPLES),
        filtered=expand("results/02_align/stats/{sample}.filtered.flagstat", sample=SAMPLES),
        final=expand("results/02_align/stats/{sample}.dedup.flagstat", sample=SAMPLES),
        frip=expand("results/qc/frip/{sample}.tsv", sample=PEAK_SAMPLES),
        size_factors=expand("results/04_bigwig/size_factors.{tag}.tsv", tag=[t for t in BW_TAGS if BW_NORM_OF[t] in BW_SF]),
        peaks=[peak_file(s) for s in PEAK_SAMPLES],
        contam=(
            expand("results/qc/contamination/{sample}.log", sample=SAMPLES)
            if CONTAM["enabled"]
            else []
        ),
    output:
        summary="results/qc/summary.tsv",
        fragments="results/qc/fragment_length.tsv",
    params:
        samples=SAMPLES,
        peak_samples=PEAK_SAMPLES,
        meta={
            s: {
                "layout": "PE" if is_pe(s) else "SE",
                "control": control_of(s),
                "peak_type": "" if s in CONTROLS else SAMPLES_DF.loc[s, "peak_type"],
                "role": "control" if s in CONTROLS else "chip",
            }
            for s in SAMPLES
        },
        trim={s: [trim_report(s, r) for r in reads_of(s)] for s in SAMPLES},
        mapq=config["filter"]["mapq"],
        mito=MITO,
        tss_matrix=config["qc"]["tss_matrix"],
        sf_tags=[t for t in BW_TAGS if BW_NORM_OF[t] in BW_SF],
    log:
        "<logs>/summary.log",
    script:
        "../scripts/summary.py"


rule multiqc:
    input:
        unpack(lambda wc: qc_files()),
        fastqc=[
            f"results/qc/{d}/{s}_{r}{suffix}_fastqc.zip"
            for s in SAMPLES
            for r in reads_of(s)
            for d, suffix in (("fastqc_raw", ""), ("fastqc_trimmed", ".trimmed"))
        ],
        trim=[
            f"results/01_trim/reports/{s}_{r}.fastq.gz_trimming_report.txt"
            for s in SAMPLES
            for r in reads_of(s)
        ],
        bowtie2=expand("results/qc/bowtie2/{sample}.log", sample=SAMPLES),
        flagstat=expand("results/02_align/stats/{sample}.dedup.flagstat", sample=SAMPLES),
        peaks=[peak_file(s) for s in PEAK_SAMPLES],
        summary="results/qc/summary.tsv",
        config=workflow.source_path("../assets/multiqc_config.yaml"),
    output:
        "results/qc/multiqc/multiqc_report.html",
    params:
        # peak 文件本身 MultiQC 不识别，传入的是 macs3 的 xls（不是规则输出，按路径传）
        files=lambda wc, input: [
            f
            for k, v in input.items()
            if k not in ("config", "peaks")
            for f in ([v] if isinstance(v, str) else v)
        ]
        + [f"results/03_peaks/macs3/{s}/{s}_peaks.xls" for s in PEAK_SAMPLES],
        title="bulk ChIP-seq",
    log:
        "<logs>/multiqc.log",
    conda:
        ENV_MAIN
    shell:
        """
        multiqc -f -c {input.config} --title '{params.title}' -o results/qc/multiqc \
            -n multiqc_report.html {params.files} > {log} 2>&1
        """



FRAG_FIG = "results/qc/figures/fragment_length"


def fragment_report_inputs():
    """QC 报告中的片段长度分布图（只有 SE 样本时不画）。"""
    return {"frag_svg": f"{FRAG_FIG}.svg"} if PE_SAMPLES else {}


rule fragment_plot:
    """PE 样本去重后的片段长度分布（R，ggplot2）：PNG 300 dpi、PDF、SVG（嵌入 QC 报告）。"""
    input:
        "results/qc/fragment_length.tsv",
    output:
        multiext(FRAG_FIG, ".png", ".pdf", ".svg"),
    params:
        groups={s: SAMPLES_DF.loc[s, "group"] for s in SAMPLES},
        controls=sorted(CONTROLS),
        max_length=600,
    log:
        "<logs>/fragment_plot.log",
    conda:
        ENV_MAIN
    script:
        "../scripts/fragment_plot.R"


def group_bigwig_desc():
    """QC 报告中组 bigwig 的说明：(中文, 英文)。"""
    if not BW_GROUP["enabled"]:
        return "不生成", "disabled"
    if BW_GROUP["method"] == "merge":
        zh, en = "合并组内去重后的 BAM 再 `bamCoverage`（参数同上）", "dedup BAMs of each group merged, then `bamCoverage` (same options)"
    else:
        zh, en = "各重复的 bigwig 逐 bin 取平均（`bigwigAverage`）", "replicate bigwigs averaged per bin (`bigwigAverage`)"
    zh += f"；{len(BW_GROUPS)} 个组"
    en += f"; {len(BW_GROUPS)} groups"
    excluded = SAMPLES_DF.index[SAMPLES_DF["merge"] == "no"].tolist()
    if excluded:
        zh += "；不参与合并：" + "、".join(excluded)
        en += "; excluded: " + ", ".join(excluded)
    return zh, en


def bigwig_desc():
    """设置表的 bigwig 一行：每种归一化一句，注明 QC 图用哪一种。"""
    b = config["bigwig"]
    zh, en = [], []
    for tag, norm in BW_NORM_OF.items():
        if norm == "deseq2_peaks":
            z, e = "DESeq2 size factor（consensus peak 计数）", "DESeq2 size factor (consensus peak counts)"
        elif norm == "deseq2_background":
            z = f"DESeq2 size factor（全基因组 {b['background_bin']} bp bin 计数）"
            e = f"DESeq2 size factor (genome-wide {b['background_bin']} bp bin counts)"
        else:
            z = e = f"`--normalizeUsing {norm}`"
        zh.append(f"`{tag}.bw`：{z}")
        en.append(f"`{tag}.bw`: {e}")
    sf = any(n in BW_SF for n in BW_NORMS)
    return (
        f"`--binSize {b['binsize']}`，排除 blacklist；" + "；".join(zh)
        + ("（size factor 为 median-of-ratios，只做 IP 样本，`--scaleFactor 1 / size factor`）" if sf else "")
        + f"；QC 图用 `{BW_TAG}`",
        f"`--binSize {b['binsize']}`, blacklist excluded; " + "; ".join(en)
        + (" (median-of-ratios size factors, IP samples only, `--scaleFactor 1 / size factor`)" if sf else "")
        + f"; QC figures use `{BW_TAG}`",
    )


def report_settings():
    """QC 报告「流程与参数」中的参数表：(中文名, 英文名, 中文取值, 英文取值)，取值中 `...` 显示为代码。"""
    c, ref, m = config, REF, config["macs3"]
    b = c["bigwig"]
    n_pe, n_se = len(PE_SAMPLES), len(SE_SAMPLES)
    ip = SAMPLES_DF.loc[PEAK_SAMPLES]
    with_ctrl = int((ip["control"] != "").sum())
    types = "、".join(f"{v}（{n}）" for v, n in ip["peak_type"].value_counts().items())
    types_en = ", ".join(f"{v} ({n})" for v, n in ip["peak_type"].value_counts().items())
    mito_zh = f"去掉 {MITO}" if REMOVE_MITO else f"保留 {MITO}"
    mito_en = f"{MITO} removed" if REMOVE_MITO else f"{MITO} kept"
    rows = [
        ("样本数", "Samples", f"{len(SAMPLES)}（IP {len(PEAK_SAMPLES)}，对照 {len(CONTROLS)}；SE {n_se}，PE {n_pe}）",
         f"{len(SAMPLES)} (IP {len(PEAK_SAMPLES)}, control {len(CONTROLS)}; SE {n_se}, PE {n_pe})"),
        ("peak 类型", "Peak type", types, types_en),
        ("对照", "Control", f"{with_ctrl} 个 IP 样本有 input / IgG 对照，{len(PEAK_SAMPLES) - with_ctrl} 个无对照",
         f"{with_ctrl} IP samples with input / IgG, {len(PEAK_SAMPLES) - with_ctrl} without"),
        ("参考基因组", "Genome", c["genome"], c["genome"]),
        ("bowtie2 index", "bowtie2 index", f"`{ref['bowtie2_index']}`", f"`{ref['bowtie2_index']}`"),
        ("GTF（TSS）", "GTF (TSS)", f"`{ref['gtf']}`；只取 protein_coding 基因", f"`{ref['gtf']}`; protein_coding genes only"),
        ("blacklist", "blacklist", f"`{ref['blacklist']}`" if ref.get("blacklist") else "无",
         f"`{ref['blacklist']}`" if ref.get("blacklist") else "none"),
        ("trim_galore", "trim_galore", f"`{c['trim']['params']}`", f"`{c['trim']['params']}`"),
        ("bowtie2", "bowtie2", f"`{c['bowtie2']['params']}`；PE 另加 `--no-mixed --no-discordant`",
         f"`{c['bowtie2']['params']}`; plus `--no-mixed --no-discordant` for PE"),
        ("过滤", "Filtering", f"MAPQ ≥ {c['filter']['mapq']}，`-F 1804`，PE 只保留 properly paired；{mito_zh}；`samtools markdup -r` 去重",
         f"MAPQ ≥ {c['filter']['mapq']}, `-F 1804`, properly paired only for PE; {mito_en}; duplicates removed with `samtools markdup -r`"),
        ("MACS3", "MACS3",
         f"q < {m['qvalue']}；有对照 `{m['chip_control']}`；无对照 `{m['chip_no_control']}`；broad 另加 `{m['broad']}`",
         f"q < {m['qvalue']}; with control `{m['chip_control']}`; without `{m['chip_no_control']}`; broad adds `{m['broad']}`"),
        ("bigwig", "bigwig", *bigwig_desc()),
        ("组 bigwig", "Group bigwig", *group_bigwig_desc()),
        ("TSS profile", "TSS profile", f"`{c['qc']['tss_matrix']}`" if c["qc"]["tss_profile"] else "不计算",
         f"`{c['qc']['tss_matrix']}`" if c["qc"]["tss_profile"] else "disabled"),
        ("Meta gene profile", "Meta gene profile",
         f"`{c['qc']['metagene_matrix']}`；protein_coding 基因，长度 ≥ {c['qc']['metagene_min_length']} bp" if c["qc"]["metagene_profile"] else "不计算",
         f"`{c['qc']['metagene_matrix']}`; protein_coding genes ≥ {c['qc']['metagene_min_length']} bp" if c["qc"]["metagene_profile"] else "disabled"),
        ("污染检查", "Contamination", f"比对到 {CONTAM['genome']}" if CONTAM["enabled"] else "不检查",
         f"aligned to {CONTAM['genome']}" if CONTAM["enabled"] else "disabled"),
        ("motif", "motif", "HOMER（narrow、非对照样本）" if c["motif"]["enabled"] else "不分析",
         "HOMER (narrow, non-control samples)" if c["motif"]["enabled"] else "disabled"),
    ]
    r = REPLICATES
    if replicates_on():
        lo, hi = r["y_frac"]
        sex_zh = (f"；性别（含对照）按 chrY 计数区 `{REF['y_region']}`（chrY 去掉 PAR）的 reads 占全部比对 reads 的比例判定"
                  f"（< {lo}% 雌性，> {hi}% 雄性）" if r["sex_check"] else "")
        sex_en = (f"; sex (controls included) called from the fraction of mapped reads in the chrY region `{REF['y_region']}` "
                  f"(chrY without the PAR; < {lo}% female, > {hi}% male)" if r["sex_check"] else "")
        rows.append((
            "重复一致性", "Replicates",
            f"IP 样本 {len(GROUPS)} 组（{len(REP_GROUPS)} 组有 ≥ 2 个重复）；reproducible peak 至少出现在 {r['min_reps']} 个重复中；"
            f"consensus peak 由 DiffBind 得到（组内 consensus 取并集，" + (f"峰顶 ± {diffbind_summits()} bp" if diffbind_summits() else "原宽度") + "）；"
            f"IP 计数的 log2(CPM + 1)（不扣除 input），{r['method']} 相关，PCA 取方差最大的 {r['pca_top']} 个 peak" + sex_zh,
            f"IP samples in {len(GROUPS)} groups ({len(REP_GROUPS)} with ≥ 2 replicates); reproducible peaks in ≥ {r['min_reps']} replicates; "
            f"consensus peaks from DiffBind (union of within-group consensus sets, " + (f"summit ± {diffbind_summits()} bp" if diffbind_summits() else "merged widths") + "); "
            f"log2(CPM + 1) of IP counts (input not subtracted), {r['method']} correlation, "
            f"PCA on the {r['pca_top']} most variable peaks" + sex_en,
        ))
    else:
        rows.append(("重复一致性", "Replicates", "不分析", "disabled"))
    return rows


rule qc_report:
    """单文件 HTML 报告：概览、流程与参数、QC 汇总表（指标说明与参考值见 workflow/report/qc_metrics.yaml）。"""
    input:
        unpack(lambda wc: replicate_report_inputs()),
        unpack(lambda wc: profile_files()),
        unpack(lambda wc: fragment_report_inputs()),
        summary="results/qc/summary.tsv",
        multiqc="results/qc/multiqc/multiqc_report.html",
        metrics=workflow.source_path("../report/qc_metrics.yaml"),
        flowchart=workflow.source_path("../../docs/images/bulk_chip_workflow.svg"),
    output:
        "results/qc/report.html",
    params:
        title=("ChIP-seq / CUT&Tag QC 报告", "ChIP-seq / CUT&Tag QC report"),
        pipeline="snakemake_bulk_chip",
        version=pipeline_version(Path(workflow.basedir).parent),
        settings=report_settings(),
        qvalue=config["macs3"]["qvalue"],
        multiqc="multiqc/multiqc_report.html",  # 相对 report.html 的路径
        rep_figures=replicate_figures() if replicates_on() else [],
        profile=dict(
            tss=config["qc"]["tss_matrix"],
            metagene=config["qc"]["metagene_matrix"],
            min_len=config["qc"]["metagene_min_length"],
            norm=BW_TAG,
            controls=BW_TAG not in BW_SF.values(),
            group_method=BW_GROUP["method"],
        ),
    log:
        "<logs>/qc_report.log",
    localrule: True
    script:
        "../scripts/qc_report.py"
