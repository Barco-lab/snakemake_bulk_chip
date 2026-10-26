# bigwig 与 deeptools QC：fingerprint（ChIP）、TSS profile / TSS enrichment、meta gene profile


def bigwig_args(samples):
    b = config["bigwig"]
    args = f"--binSize {b['binsize']} --normalizeUsing {b['normalize']}"
    if b["normalize"] == "RPGC":
        args += f" --effectiveGenomeSize {REF['effective_genome_size']}"
    if REF.get("blacklist"):
        args += f" --blackListFileName {REF['blacklist']}"
    # PE 按真实片段延伸；SE 延伸到估计片段长度（PE、SE 合并时 PE 仍按真实片段）
    args += " --extendReads" if all(is_pe(s) for s in samples) else f" --extendReads {b['se_fragment_length']}"
    return f"{args} {b['extra']}"


def group_bigwig_cmd(wc, input):
    """merge：合并后的 BAM 按单样本参数 bamCoverage；mean：各重复的 bigwig 逐 bin 取平均。"""
    if BW_GROUP["method"] == "merge":
        return f"bamCoverage -b {input.bam} {bigwig_args(BW_GROUPS[wc.group])}"
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
        return {"bam": bam, "bai": f"{bam}.bai"}
    return {"bw": expand("results/04_bigwig/{sample}.bw", sample=BW_GROUPS[wc.group])}


rule bigwig:
    input:
        bam="results/02_align/{sample}.dedup.bam",
        bai="results/02_align/{sample}.dedup.bam.bai",
    output:
        "results/04_bigwig/{sample}.bw",
    params:
        args=lambda wc: bigwig_args([wc.sample]),
    threads: threads_for("deeptools")
    log:
        "<logs>/bigwig/{sample}.log",
    benchmark:
        "<benchmarks>/bigwig/{sample}.tsv"
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
        "results/04_bigwig/groups/{group}.bw",
    params:
        cmd=group_bigwig_cmd,
    threads: threads_for("deeptools")
    log:
        "<logs>/bigwig/groups/{group}.log",
    benchmark:
        "<benchmarks>/bigwig/groups/{group}.tsv"
    conda:
        ENV_MACS
    shell:
        "{params.cmd} -o {output} -p {threads} > {log} 2>&1"


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
        bw=expand("results/04_bigwig/{sample}.bw", sample=SAMPLES),
        tss="results/reference/tss.bed",
    output:
        matrix="results/qc/deeptools/tss_matrix.gz",
        png="results/qc/deeptools/tss_profile.png",
        data="results/qc/deeptools/tss_profile.tsv",
    params:
        labels=SAMPLES,
        plot=profile_plot_args(len(SAMPLES)),
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
    if wc.set == "groups":
        return expand("results/04_bigwig/groups/{group}.bw", group=BW_GROUPS)
    return expand("results/04_bigwig/{sample}.bw", sample=SAMPLES)


rule metagene_profile:
    """meta gene profile：samples 为每个样本的 bigwig（含对照），groups 为组 bigwig。"""
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
        labels=lambda wc: list(BW_GROUPS) if wc.set == "groups" else SAMPLES,
        plot=lambda wc: profile_plot_args(len(BW_GROUPS) if wc.set == "groups" else len(SAMPLES)),
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
