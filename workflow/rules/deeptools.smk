# bigwig 与 deeptools QC：fingerprint（ChIP）、TSS profile / TSS enrichment、meta gene profile
# bigwig 每种归一化一份（bigwig.normalize，见 common.smk 的 BW_TAGS）；DESeq2 size factor 的两种方式先算 size factor


def extend_args(samples):
    # PE 按真实片段延伸；SE 延伸到估计片段长度（PE、SE 合并时 PE 仍按真实片段）
    b = config["bigwig"]
    return "--extendReads" if all(is_pe(s) for s in samples) else f"--extendReads {b['se_fragment_length']}"


def bigwig_args(samples, tag, sf_table=None):
    """bamCoverage 参数；size factor 方式为 --normalizeUsing None --scaleFactor 1 / Σ size factor（组 bigwig 为组内之和，
    即组内各样本归一化信号的平均），size factor 在运行时从 sf_table 读取。"""
    b = config["bigwig"]
    norm = BW_NORM_OF[tag]
    if norm in BW_SF:
        names = ",".join(samples)
        args = ("--normalizeUsing None --scaleFactor "
                f"$(awk -F'\\t' -v s={names} 'BEGIN {{n = split(s, a, \",\"); for (i = 1; i <= n; i++) w[a[i]] = 1}} "
                f"w[$1] {{t += $2}} END {{printf \"%.6g\", 1 / t}}' {sf_table})")
    else:
        args = f"--normalizeUsing {norm}"
        if norm == "RPGC":
            args += f" --effectiveGenomeSize {REF['effective_genome_size']}"
    args = f"--binSize {b['binsize']} {args}"
    if REF.get("blacklist"):
        args += f" --blackListFileName {REF['blacklist']}"
    return f"{args} {extend_args(samples)} {b['extra']}"


def sf_input(wc):
    """size factor 方式的 bigwig 需要 size factor 表。"""
    return {"sf": f"results/04_bigwig/size_factors.{wc.norm}.tsv"} if BW_NORM_OF[wc.norm] in BW_SF else {}


def group_bigwig_cmd(wc, input):
    """merge：合并后的 BAM 按单样本参数 bamCoverage；mean：各重复的 bigwig 逐 bin 取平均。"""
    if BW_GROUP["method"] == "merge":
        return f"bamCoverage -b {input.bam} {bigwig_args(BW_GROUPS[wc.group], wc.norm, input.get('sf'))}"
    args = f"--binSize {config['bigwig']['binsize']}"
    if REF.get("blacklist"):
        args += f" --blackListFileName {REF['blacklist']}"
    return f"bigwigAverage -b {' '.join(input.bw)} {args}"


def profile_plot_args(n):
    """plotProfile 画在一张图上；图高随曲线数增加（cm），图例放在右上角（TES 一侧信号低）。"""
    return (f"--perGroup --plotHeight {max(7, round(0.8 * n, 1))} --plotWidth 12 "
            "--legendLocation upper-right --regionsLabel 'protein-coding genes'")


def group_bigwig_input(wc):
    if BW_GROUP["method"] == "merge":
        bam = f"results/04_bigwig/groups/{wc.group}.merged.bam"
        return {"bam": bam, "bai": f"{bam}.bai", **sf_input(wc)}
    return {"bw": sample_bw(wc.norm, BW_GROUPS[wc.group])}


rule bigwig:
    input:
        unpack(sf_input),
        bam="results/02_align/{sample}.dedup.bam",
        bai="results/02_align/{sample}.dedup.bam.bai",
    output:
        "results/04_bigwig/{sample}.{norm}.bw",
    params:
        args=lambda wc, input: bigwig_args([wc.sample], wc.norm, input.get("sf")),
    threads: threads_for("deeptools")
    log:
        "<logs>/bigwig/{sample}.{norm}.log",
    benchmark:
        "<benchmarks>/bigwig/{sample}.{norm}.tsv"
    conda:
        ENV_MACS
    shell:
        "bamCoverage -b {input.bam} -o {output} -p {threads} {params.args} > {log} 2>&1"


rule group_bam:
    """组 bigwig（method: merge）：合并组内参与合并的去重 BAM，生成 bigwig 后自动删除。"""
    input:
        bams=lambda wc: expand("results/02_align/{sample}.dedup.bam", sample=BW_GROUPS[wc.group]),
    output:
        bam=temp("results/04_bigwig/groups/{group}.merged.bam"),
        bai=temp("results/04_bigwig/groups/{group}.merged.bam.bai"),
    threads: threads_for("filter")
    log:
        "<logs>/bigwig/groups/{group}.merge.log",
    benchmark:
        "<benchmarks>/bigwig/groups/{group}.merge.tsv"
    conda:
        ENV_MACS
    shell:
        """
        samtools merge -f -@ {threads} -o {output.bam} {input.bams} > {log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """


rule group_bigwig:
    """每组一个 bigwig：样本表 merge 列为 no 的重复不参与；方法见 bigwig.group.method。"""
    input:
        unpack(group_bigwig_input),
    output:
        "results/04_bigwig/groups/{group}.{norm}.bw",
    params:
        cmd=group_bigwig_cmd,
    threads: threads_for("deeptools")
    log:
        "<logs>/bigwig/groups/{group}.{norm}.log",
    benchmark:
        "<benchmarks>/bigwig/groups/{group}.{norm}.tsv"
    conda:
        ENV_MACS
    shell:
        "{params.cmd} -o {output} -p {threads} > {log} 2>&1"


rule bin_counts:
    """deseq2_background：IP 样本在全基因组 bin（bigwig.background_bin）上的 reads 数（排除 blacklist，片段延伸同 bigwig）。"""
    input:
        bam=expand("results/02_align/{sample}.dedup.bam", sample=PEAK_SAMPLES),
        bai=expand("results/02_align/{sample}.dedup.bam.bai", sample=PEAK_SAMPLES),
    output:
        counts=temp("results/04_bigwig/bin_counts.SFbg.tsv"),
        npz=temp("results/04_bigwig/bin_counts.SFbg.npz"),
    params:
        labels=PEAK_SAMPLES,
        args=" ".join([
            f"--binSize {config['bigwig']['background_bin']}",
            f"--blackListFileName {REF['blacklist']}" if REF.get("blacklist") else "",
            extend_args(PEAK_SAMPLES),
        ]),
    threads: threads_for("deeptools")
    log:
        "<logs>/bigwig/bin_counts.SFbg.log",
    benchmark:
        "<benchmarks>/bigwig/bin_counts.SFbg.tsv"
    conda:
        ENV_MACS
    shell:
        """
        multiBamSummary bins -b {input.bam} --labels {params.labels} {params.args} -p {threads} \
            -o {output.npz} --outRawCounts {output.counts} > {log} 2>&1
        """


def sf_counts(wc):
    # SFpeaks：DiffBind 在 consensus peak 上的计数；SFbg：全基因组 bin 的计数
    if wc.sf == "SFpeaks":
        return "results/qc/replicates/peak_counts.tsv"
    return "results/04_bigwig/bin_counts.SFbg.tsv"


rule size_factors:
    """DESeq2 size factor（median-of-ratios，只用 IP 样本）；bigwig 除以它（--scaleFactor 1 / size factor）。"""
    input:
        counts=sf_counts,
    output:
        "results/04_bigwig/size_factors.{sf}.tsv",
    params:
        samples=PEAK_SAMPLES,
    log:
        "<logs>/bigwig/size_factors.{sf}.log",
    conda:
        ENV_MAIN
    script:
        "../scripts/size_factors.R"


rule fingerprint:
    """ChIP 富集程度：IP 曲线越偏离对角线（相对 input）越好。"""
    input:
        bam=expand("results/02_align/{sample}.dedup.bam", sample=SAMPLES),
    output:
        png="results/qc/deeptools/fingerprint.png",
        counts="results/qc/deeptools/fingerprint.raw_counts.tsv",
        metrics="results/qc/deeptools/fingerprint.qc_metrics.tsv",
    params:
        labels=SAMPLES,
    threads: threads_for("deeptools")
    log:
        "<logs>/deeptools/fingerprint.log",
    conda:
        ENV_MACS
    shell:
        """
        plotFingerprint -b {input.bam} --labels {params.labels} --skipZeros -p {threads} \
            --plotFile {output.png} --outRawCounts {output.counts} \
            --outQualityMetrics {output.metrics} > {log} 2>&1
        """


rule tss_bed:
    """
    从 GTF 提取 protein_coding 基因的 TSS（+ 链取 start，- 链取 end），染色体名转换为 bowtie2 index 的命名。
    只用 protein_coding，TSS enrichment 不受 GTF 中非编码基因数量的影响。
    """
    input:
        gtf=REF["gtf"],
        chroms="results/reference/index_chroms.txt",
    output:
        "results/reference/tss.bed",
    params:
        rename=RENAME_CHROMS,
    log:
        "<logs>/tss_bed.log",
    localrule: True
    shell:
        """
        if [[ {input.gtf} == *.gz ]]; then zcat {input.gtf}; else cat {input.gtf}; fi \
        | awk -F'\\t' 'BEGIN {{OFS = "\\t"}} $3 == "gene" && $9 ~ /gene_(bio)?type "protein_coding"/ {{
            match($9, /gene_id "[^"]+"/); id = substr($9, RSTART + 9, RLENGTH - 10);
            tss = ($7 == "+") ? $4 - 1 : $5 - 1;
            print $1, tss, tss + 1, id, ".", $7 }}' \
        | awk -F'\\t' -v OFS='\\t' -v col=1 -f {params.rename} {input.chroms} - 2> {log} \
        | sort -k1,1 -k2,2n -u > {output}
        """


rule tss_profile:
    input:
        bw=sample_bw(),
        tss="results/reference/tss.bed",
    output:
        matrix="results/qc/deeptools/tss_matrix.gz",
        png="results/qc/deeptools/tss_profile.png",
        data="results/qc/deeptools/tss_profile.tsv",
    params:
        labels=bw_samples(BW_TAG),
        plot=profile_plot_args(len(bw_samples(BW_TAG))),
        cmd=config["qc"]["tss_matrix"],
    threads: threads_for("deeptools")
    log:
        "<logs>/deeptools/tss_profile.log",
    conda:
        ENV_MACS
    shell:
        """
        {params.cmd} -S {input.bw} -R {input.tss} --samplesLabel {params.labels} \
            -p {threads} -o {output.matrix} > {log} 2>&1
        plotProfile -m {output.matrix} -o {output.png} --outFileNameData {output.data} \
            {params.plot} >> {log} 2>&1
        """


rule gene_bed:
    """meta gene profile 用的 protein_coding 基因（TSS 到 TES，长度 ≥ qc.metagene_min_length），染色体名转换同 tss_bed。"""
    input:
        gtf=REF["gtf"],
        chroms="results/reference/index_chroms.txt",
    output:
        "results/reference/genes.bed",
    params:
        rename=RENAME_CHROMS,
        min_len=config["qc"]["metagene_min_length"],
    log:
        "<logs>/gene_bed.log",
    localrule: True
    shell:
        """
        if [[ {input.gtf} == *.gz ]]; then zcat {input.gtf}; else cat {input.gtf}; fi \
        | awk -F'\\t' -v min={params.min_len} 'BEGIN {{OFS = "\\t"}}
            $3 == "gene" && $9 ~ /gene_(bio)?type "protein_coding"/ && $5 - $4 + 1 >= min {{
            match($9, /gene_id "[^"]+"/); id = substr($9, RSTART + 9, RLENGTH - 10);
            print $1, $4 - 1, $5, id, ".", $7 }}' \
        | awk -F'\\t' -v OFS='\\t' -v col=1 -f {params.rename} {input.chroms} - 2> {log} \
        | sort -k1,1 -k2,2n -u > {output}
        """


def metagene_bw(wc):
    return group_bw() if wc.set == "groups" else sample_bw()


rule metagene_profile:
    """meta gene profile：samples 为每个样本的 bigwig（含对照；size factor 方式不含），groups 为组 bigwig。"""
    input:
        bw=metagene_bw,
        genes="results/reference/genes.bed",
    output:
        matrix="results/qc/deeptools/metagene_{set}.matrix.gz",
        png="results/qc/deeptools/metagene_{set}.png",
        data="results/qc/deeptools/metagene_{set}.tsv",
    wildcard_constraints:
        set="samples|groups",
    params:
        labels=lambda wc: list(BW_GROUPS) if wc.set == "groups" else bw_samples(BW_TAG),
        plot=lambda wc: profile_plot_args(len(BW_GROUPS) if wc.set == "groups" else len(bw_samples(BW_TAG))),
        cmd=config["qc"]["metagene_matrix"],
    threads: threads_for("deeptools")
    log:
        "<logs>/deeptools/metagene_{set}.log",
    conda:
        ENV_MACS
    shell:
        """
        {params.cmd} -S {input.bw} -R {input.genes} --samplesLabel {params.labels} \
            -p {threads} -o {output.matrix} > {log} 2>&1
        plotProfile -m {output.matrix} -o {output.png} --outFileNameData {output.data} \
            {params.plot} >> {log} 2>&1
        """
