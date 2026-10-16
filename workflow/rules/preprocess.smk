# 原始数据 -> trim_galore ；raw / trimmed 各做一次 FastQC


rule link_raw:
    """把原始 fastq 统一链接为 {sample}_R1/R2.fastq.gz，后续文件名、报告中的样本名都与之一致。"""
    input:
        lambda wc: SAMPLES_DF.loc[wc.sample, "fq1" if wc.read == "R1" else "fq2"],
    output:
        "results/00_raw/{sample}_{read}.fastq.gz",
    localrule: True
    shell:
        "ln -sf {input:q} {output:q}"


TRIM_ARGS = config["trim"]["params"]


rule trim_pe:
    wildcard_constraints:
        sample=regex_of(PE_SAMPLES),
    input:
        r1="results/00_raw/{sample}_R1.fastq.gz",
        r2="results/00_raw/{sample}_R2.fastq.gz",
    output:
        r1=maybe_temp("results/01_trim/{sample}_R1.trimmed.fq.gz", "trimmed_fastq"),
        r2=maybe_temp("results/01_trim/{sample}_R2.trimmed.fq.gz", "trimmed_fastq"),
        report_r1="results/01_trim/reports/{sample}_R1.fastq.gz_trimming_report.txt",
        report_r2="results/01_trim/reports/{sample}_R2.fastq.gz_trimming_report.txt",
    params:
        args=TRIM_ARGS,
        tmp="results/01_trim/{sample}.tmp",
    threads: threads_for("trim")
    log:
        "<logs>/trim/{sample}.log",
    benchmark:
        "<benchmarks>/trim/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {params.tmp} results/01_trim/reports
        trim_galore {params.args} --cores {threads} --paired --basename {wildcards.sample} \
            -o {params.tmp} {input.r1} {input.r2} > {log} 2>&1
        mv {params.tmp}/{wildcards.sample}_val_1.fq.gz {output.r1}
        mv {params.tmp}/{wildcards.sample}_val_2.fq.gz {output.r2}
        mv {params.tmp}/*_trimming_report.txt results/01_trim/reports/
        rm -rf {params.tmp}
        """


rule trim_se:
    wildcard_constraints:
        sample=regex_of(SE_SAMPLES),
    input:
        r1="results/00_raw/{sample}_R1.fastq.gz",
    output:
        r1=maybe_temp("results/01_trim/{sample}_R1.trimmed.fq.gz", "trimmed_fastq"),
        report_r1="results/01_trim/reports/{sample}_R1.fastq.gz_trimming_report.txt",
    params:
        args=TRIM_ARGS,
        tmp="results/01_trim/{sample}.tmp",
    threads: threads_for("trim")
    log:
        "<logs>/trim/{sample}.log",
    benchmark:
        "<benchmarks>/trim/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {params.tmp} results/01_trim/reports
        trim_galore {params.args} --cores {threads} --basename {wildcards.sample} \
            -o {params.tmp} {input.r1} > {log} 2>&1
        mv {params.tmp}/{wildcards.sample}_trimmed.fq.gz {output.r1}
        mv {params.tmp}/*_trimming_report.txt results/01_trim/reports/
        rm -rf {params.tmp}
        """


FASTQC_JAVA = f"--java {config['fastqc_java']}" if config.get("fastqc_java") else ""


rule fastqc_raw:
    input:
        "results/00_raw/{sample}_{read}.fastq.gz",
    output:
        html="results/qc/fastqc_raw/{sample}_{read}_fastqc.html",
        zip="results/qc/fastqc_raw/{sample}_{read}_fastqc.zip",
    params:
        java=FASTQC_JAVA,
    threads: threads_for("fastqc")
    log:
        "<logs>/fastqc_raw/{sample}_{read}.log",
    conda:
        ENV_MAIN
    shell:
        "fastqc {params.java} -t {threads} -o results/qc/fastqc_raw {input} > {log} 2>&1"


rule fastqc_trimmed:
    input:
        "results/01_trim/{sample}_{read}.trimmed.fq.gz",
    output:
        html="results/qc/fastqc_trimmed/{sample}_{read}.trimmed_fastqc.html",
        zip="results/qc/fastqc_trimmed/{sample}_{read}.trimmed_fastqc.zip",
    params:
        java=FASTQC_JAVA,
    threads: threads_for("fastqc")
    log:
        "<logs>/fastqc_trimmed/{sample}_{read}.log",
    conda:
        ENV_MAIN
    shell:
        "fastqc {params.java} -t {threads} -o results/qc/fastqc_trimmed {input} > {log} 2>&1"
