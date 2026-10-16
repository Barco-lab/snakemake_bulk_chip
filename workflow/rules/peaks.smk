# macs3 call peak -> 去除 blacklist 区域 -> FRiP ；可选: HOMER motif


def macs3_input(wc):
    d = {"treat": f"results/02_align/{wc.sample}.bam"}
    if control_of(wc.sample):
        d["control"] = f"results/02_align/{control_of(wc.sample)}.bam"
    return d


def macs3_args(wc):
    m = config["macs3"]
    s = wc.sample
    if control_of(s):
        args = f"{m['chip_control']} -c results/02_align/{control_of(s)}.bam"
    else:
        args = m["chip_no_control"]
    if wc.ptype == "broad":
        args += f" {m['broad']}"
    fmt = "BAMPE" if is_pe(s) else "BAM"
    return f"-f {fmt} -g {REF['effective_genome_size']} -q {m['qvalue']} {args}"


rule macs3:
    input:
        unpack(macs3_input),
    output:
        peak="results/03_peaks/macs3/{sample}/{sample}_peaks.{ptype}Peak",
    params:
        args=macs3_args,
        outdir="results/03_peaks/macs3/{sample}",
    log:
        "<logs>/macs3/{sample}.{ptype}.log",
    benchmark:
        "<benchmarks>/macs3/{sample}.{ptype}.tsv"
    conda:
        ENV_MACS
    shell:
        """
        macs3 callpeak -t {input.treat} {params.args} \
            -n {wildcards.sample} --outdir {params.outdir} > {log} 2>&1
        """


rule filter_blacklist:
    input:
        "results/03_peaks/macs3/{sample}/{sample}_peaks.{ptype}Peak",
    output:
        "results/03_peaks/{sample}.{ptype}Peak",
    params:
        blacklist=REF.get("blacklist", ""),
    conda:
        ENV_MACS
    shell:
        """
        if [ -n "{params.blacklist}" ]; then
            bedtools intersect -v -a {input} -b {params.blacklist} > {output}
        else
            cp {input} {output}
        fi
        """


rule frip:
    """FRiP = 落在 peak 内的 reads / 全部 reads（过滤去重后的 BAM）。"""
    input:
        bam="results/02_align/{sample}.bam",
        peak=lambda wc: peak_file(wc.sample),
    output:
        "results/qc/frip/{sample}.tsv",
    conda:
        ENV_MACS
    shell:
        """
        total=$(samtools view -c {input.bam})
        if [ -s {input.peak} ]; then
            inpeak=$(bedtools intersect -u -a {input.bam} -b {input.peak} -ubam | samtools view -c -)
        else
            inpeak=0
        fi
        awk -v s={wildcards.sample} -v a=$inpeak -v t=$total 'BEGIN {{
            printf "sample\\treads_in_peaks\\ttotal_reads\\tfrip\\n%s\\t%d\\t%d\\t%.4f\\n", s, a, t, (t ? a / t : 0)
        }}' > {output}
        """


rule homer:
    """以 summit ± size/2 做 motif 分析；需要 references.<genome>.fasta。"""
    input:
        "results/03_peaks/{sample}.narrowPeak",
    output:
        "results/05_motif/{sample}/knownResults.txt",
    params:
        fasta=REF.get("fasta", ""),
        half=config["motif"]["size"] // 2,
        extra=config["motif"]["params"],
        outdir="results/05_motif/{sample}",
    threads: threads_for("motif")
    log:
        "<logs>/homer/{sample}.log",
    benchmark:
        "<benchmarks>/homer/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        bed={params.outdir}/summits.bed
        awk 'BEGIN {{OFS = "\\t"}} {{ s = $2 + $10; b = s - {params.half}; if (b < 0) b = 0;
            print $1, b, s + {params.half}, $4, ".", "+" }}' {input} > $bed
        if [ -s $bed ]; then
            findMotifsGenome.pl $bed {params.fasta} {params.outdir} -size given -p {threads} \
                -preparsedDir {params.outdir}/preparsed {params.extra} > {log} 2>&1
            rm -rf {params.outdir}/preparsed
        else
            echo "没有 peak，跳过 motif 分析" > {log}
            touch {output}
        fi
        """
