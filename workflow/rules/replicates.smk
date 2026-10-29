# 重复一致性：组内 peak 重现性 -> consensus peak（DiffBind）与计数 -> 样本相关性、PCA；性别检查（chrX / chrY reads 占比）
# 分组见 samples.tsv 的 group 列（为空时按样本名去掉 _repN 推断）
# 只比较 IP 样本（PEAK_SAMPLES）；被用作 control 的样本只做性别检查


def replicates_on():
    return REPLICATES["enabled"] and len(PEAK_SAMPLES) >= 2


def replicate_figures():
    """报告中的图（不含扩展名）；条件不满足的图不画。"""
    f = ["correlation_scatter"]
    if len(PEAK_SAMPLES) >= 3:
        f.append("pca")
    if REP_GROUPS:
        f.append("peak_overlap")
    if venn_groups():
        f.append("peak_venn")
    if group_venn_sets():
        f.append("group_venn")
    if REPLICATES["sex_check"]:
        f.append("sex")
    return f


def venn_groups():
    """画 Venn 图的组：2–4 个重复（更多时 Venn 难以辨认，只看 Jaccard）。"""
    return [g for g in REP_GROUPS if len(GROUPS[g]) <= 4]


def group_venn_sets():
    """组间 Venn：{名称: [组]}；replicates.group_venn 为空时，共 2–4 个组则全部比较。"""
    sets = REPLICATES.get("group_venn") or ([list(GROUPS)] if 2 <= len(GROUPS) <= 4 else [])
    out = {}
    for gs in sets:
        bad = [g for g in gs if g not in GROUPS]
        if bad:
            raise ValueError(f"replicates.group_venn 中的组不存在（或不是 IP 样本的组）: {bad}，可选: {list(GROUPS)}")
        if not 2 <= len(gs) <= 4 or len(set(gs)) != len(gs):
            raise ValueError(f"replicates.group_venn 每项须为 2–4 个不同的组: {gs}")
        out[" vs ".join(gs)] = list(gs)
    return out


def replicate_figure_files():
    return expand(
        "results/qc/replicates/figures/{fig}.{ext}", fig=replicate_figures(), ext=["png", "pdf", "svg"]
    )


def replicate_report_inputs():
    """QC 报告「重复一致性」一节的输入；未启用时为空。"""
    if not replicates_on():
        return {}
    d = "results/qc/replicates"
    return dict(
        rep_info=f"{d}/info.json",
        rep_samples=f"{d}/samples.tsv",
        rep_groups=f"{d}/groups.tsv",
        rep_svg=expand(f"{d}/figures/{{fig}}.svg", fig=replicate_figures()),
    )


def region_len(region):
    """chr:start-end（1-based，含两端）的长度。"""
    start, end = region.split(":")[1].replace(",", "").split("-")
    return int(end) - int(start) + 1


rule sex_counts:
    """全部比对 reads、chrX 与 chrY 计数区（references 中的 y_region）上的 reads 数（过滤去重后的 BAM）。"""
    input:
        bam="results/02_align/{sample}.dedup.bam",
        bai="results/02_align/{sample}.dedup.bam.bai",
    output:
        "results/qc/replicates/sex/{sample}.tsv",
    params:
        x=REF["chrX"],
        y=REF["y_region"],
        y_len=region_len(REF["y_region"]),
    conda:
        ENV_MACS
    shell:
        """
        y=$(samtools view -c {input.bam} {params.y})
        samtools idxstats {input.bam} | awk -v OFS='\\t' -v s={wildcards.sample} -v x={params.x} -v y=$y -v yl={params.y_len} '
            {{ t += $3 }}
            $1 == x {{ xr = $3 }}
            END {{
                print "sample", "total_reads", "x_reads", "y_reads", "y_len"
                print s, t, xr + 0, y, yl
            }}' > {output}
        """


def diffbind_summits():
    """narrow 取峰顶 ± summits bp；有 broad 的组时整体不重新定位（宽峰没有单一峰顶）。"""
    if any(SAMPLES_DF.loc[s, "peak_type"] == "broad" for s in PEAK_SAMPLES):
        return False
    return REPLICATES["summits"]


rule diffbind_consensus:
    """DiffBind：组内 consensus（至少 min_reps 个重复）-> 各组并集 -> dba.count 计数（只用 IP 样本，不扣除对照）。"""
    input:
        bams=expand("results/02_align/{sample}.dedup.bam", sample=PEAK_SAMPLES),
        bais=expand("results/02_align/{sample}.dedup.bam.bai", sample=PEAK_SAMPLES),
        peaks=[peak_file(s) for s in PEAK_SAMPLES],
    output:
        reproducible=expand("results/03_peaks/reproducible/{group}.bed", group=REP_GROUPS),
        consensus="results/03_peaks/consensus.bed",
        counts="results/qc/replicates/peak_counts.tsv",
        rds="results/qc/replicates/diffbind.rds",
    params:
        samples=PEAK_SAMPLES,
        groups=GROUPS,
        rep_groups=REP_GROUPS,
        min_reps=REPLICATES["min_reps"],
        summits=diffbind_summits(),
        fragment_size=config["bigwig"]["se_fragment_length"],
    threads: threads_for("deeptools")
    log:
        "<logs>/diffbind_consensus.log",
    benchmark:
        "<benchmarks>/diffbind_consensus.tsv"
    conda:
        ENV_MAIN
    script:
        "../scripts/diffbind_consensus.R"


rule replicate_qc:
    """peak 重叠（Jaccard）、相关性、PCA、性别；结果为 tsv，绘图见 replicate_plots。"""
    input:
        counts="results/qc/replicates/peak_counts.tsv",
        peaks=[peak_file(s) for s in PEAK_SAMPLES],
        reproducible=expand("results/03_peaks/reproducible/{group}.bed", group=REP_GROUPS),
        consensus="results/03_peaks/consensus.bed",
        sex=(
            expand("results/qc/replicates/sex/{sample}.tsv", sample=SAMPLES)
            if REPLICATES["sex_check"]
            else []
        ),
    output:
        samples="results/qc/replicates/samples.tsv",
        groups="results/qc/replicates/groups.tsv",
        overlap="results/qc/replicates/peak_overlap.tsv",
        venn="results/qc/replicates/peak_venn.tsv",
        group_venn="results/qc/replicates/group_venn.tsv",
        correlation="results/qc/replicates/correlation.tsv",
        log2cpm="results/qc/replicates/log2cpm.tsv",
        pca="results/qc/replicates/pca.tsv",
        info="results/qc/replicates/info.json",
    params:
        samples=PEAK_SAMPLES,
        group_venn=group_venn_sets(),
        groups=GROUPS,
        controls=CONTROLS,
        group_of=dict(zip(SAMPLES, SAMPLES_DF["group"])),
        rep_groups=REP_GROUPS,
        method=REPLICATES["method"],
        pca_top=REPLICATES["pca_top"],
        min_reps=REPLICATES["min_reps"],
        summits=diffbind_summits(),
        y_frac=REPLICATES["y_frac"],
        y_region=REF["y_region"],
        genome=config["genome"],
        mapq=config["filter"]["mapq"],
    log:
        "<logs>/replicate_qc.log",
    conda:
        ENV_MACS
    script:
        "../scripts/replicate_qc.py"


rule replicate_plots:
    """R（ggplot2）绘图：PNG 300 dpi、PDF、SVG（嵌入 QC 报告）。"""
    input:
        samples="results/qc/replicates/samples.tsv",
        groups="results/qc/replicates/groups.tsv",
        overlap="results/qc/replicates/peak_overlap.tsv",
        venn="results/qc/replicates/peak_venn.tsv",
        group_venn="results/qc/replicates/group_venn.tsv",
        correlation="results/qc/replicates/correlation.tsv",
        log2cpm="results/qc/replicates/log2cpm.tsv",
        pca="results/qc/replicates/pca.tsv",
        info="results/qc/replicates/info.json",
    output:
        replicate_figure_files(),
    params:
        figures=replicate_figures(),
        venn_groups=venn_groups(),
        group_venn=group_venn_sets(),
        outdir="results/qc/replicates/figures",
    log:
        "<logs>/replicate_plots.log",
    conda:
        ENV_MAIN
    script:
        "../scripts/replicate_plots.R"
