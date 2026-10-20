# bigwig 与 deeptools QC：fingerprint（ChIP）、TSS profile / TSS enrichment


def bigwig_args(wc):
    b = config["bigwig"]
    args = f"--binSize {b['binsize']} --normalizeUsing {b['normalize']}"
    if b["normalize"] == "RPGC":
        args += f" --effectiveGenomeSize {REF['effective_genome_size']}"
    if REF.get("blacklist"):
        args += f" --blackListFileName {REF['blacklist']}"
    # PE 按真实片段延伸；SE 延伸到估计片段长度
    args += " --extendReads" if is_pe(wc.sample) else f" --extendReads {b['se_fragment_length']}"
    return f"{args} {b['extra']}"


rule bigwig:
    input:
        bam="results/02_align/{sample}.bam",
        bai="results/02_align/{sample}.bam.bai",
    output:
        "results/04_bigwig/{sample}.bw",
    params:
        args=bigwig_args,
    threads: threads_for("deeptools")
    log:
        "<logs>/bigwig/{sample}.log",
    benchmark:
        "<benchmarks>/bigwig/{sample}.tsv"
    conda:
        ENV_MACS
    shell:
        "bamCoverage -b {input.bam} -o {output} -p {threads} {params.args} > {log} 2>&1"


rule fingerprint:
    """ChIP 富集程度：IP 曲线越偏离对角线（相对 input）越好。"""
    input:
        bam=expand("results/02_align/{sample}.bam", sample=SAMPLES),
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
    """从 GTF 的 gene 行提取 TSS（+ 链取 start，- 链取 end），染色体名转换为 bowtie2 index 的命名。"""
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
        | awk -F'\\t' 'BEGIN {{OFS = "\\t"}} $3 == "gene" {{
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
        flank=config["qc"]["tss_flank"],
    threads: threads_for("deeptools")
    log:
        "<logs>/deeptools/tss_profile.log",
    conda:
        ENV_MACS
    shell:
        """
        computeMatrix reference-point -S {input.bw} -R {input.tss} --samplesLabel {params.labels} \
            -a {params.flank} -b {params.flank} --binSize 50 --skipZeros -p {threads} \
            -o {output.matrix} > {log} 2>&1
        plotProfile -m {output.matrix} -o {output.png} --outFileNameData {output.data} \
            --perGroup >> {log} 2>&1
        """
