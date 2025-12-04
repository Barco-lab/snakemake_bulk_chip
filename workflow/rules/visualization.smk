rule bam2bw:
    input:
        dedu_bam = get_dedu_list,
        uniq_bam = get_uniq_list
    output:
        dedu_bw = "Results/06_visualization/{sample}_dedu_{binsize}.bw",
        uniq_bw="Results/06_visualization/{sample}_unique_{binsize}.bw"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        thread = config["threads"]
    log:
        log = "Results/06_visualization/logs/{sample}_{binsize}.log"
    shell:
        """
        # bam2bw
        bamCoverage --bam {input.dedu_bam} --outFileName {output.dedu_bw} --outFileFormat bigwig \
              --binSize {wildcards.binsize} --normalizeUsing RPKM --numberOfProcessors {params.thread} >> {log.log} 2>&1
        bamCoverage --bam {input.uniq_bam} --outFileName {output.uniq_bw} --outFileFormat bigwig \
              --binSize {wildcards.binsize} --normalizeUsing RPKM --numberOfProcessors {params.thread} >> {log.log} 2>&1
        """

rule bam2bed:
    input:
        dedu_bam = get_dedu_list,
        uniq_bam = get_uniq_list
    output:
        dedu_bed = "Results/06_visualization/{sample}_dedu_{binsize}.bed",
        uniq_bed = "Results/06_visualization/{sample}_unique_{binsize}.bed"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        thread = config["threads"]
    log:
        log = "Results/06_visualization/logs/{sample}_{binsize}.log"
    shell:
        """
        # bam2bed
        bamCoverage --bam {input.dedu_bam} --outFileName {output.dedu_bed} --outFileFormat bedgraph \
              --binSize {wildcards.binsize} --normalizeUsing RPKM --numberOfProcessors {params.thread} >> {log.log} 2>&1
        bamCoverage --bam {input.uniq_bam} --outFileName {output.uniq_bed} --outFileFormat bedgraph \
              --binSize {wildcards.binsize} --normalizeUsing RPKM --numberOfProcessors {params.thread} >> {log.log} 2>&1
        """