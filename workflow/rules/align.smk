# 参考检查 -> bowtie2 -> 过滤（MAPQ、properly paired、chrM）-> 去重 ；可选: 污染检查


rule check_reference:
    """bowtie2 index 与 blacklist 的染色体命名必须一致（chr1 vs 1）；GTF / fasta 使用时转换命名，
    这里检查转换后能对上。不一致时尽早报错。"""
    output:
        check="results/qc/reference_check.txt",
        index_chroms="results/reference/index_chroms.txt",
        fasta_chroms="results/reference/fasta_chroms.txt",
    params:
        index=REF["bowtie2_index"],
        blacklist=REF.get("blacklist", ""),
        gtf=REF["gtf"],
        fasta=REF.get("fasta", ""),
        rename=RENAME_CHROMS,
    log:
        "<logs>/check_reference.log",
    conda:
        ENV_MAIN
    shell:
        """
        bowtie2-inspect -n {params.index} | awk '{{print $1}}' | sort -u > {output.index_chroms}
        if [ -n "{params.fasta}" ]; then
            grep '^>' {params.fasta} | awk '{{print substr($1, 2)}}' | sort -u > {output.fasta_chroms}
        else
            : > {output.fasta_chroms}
        fi
        cat_any() {{ if [[ $1 == *.gz ]]; then zcat "$1"; else cat "$1"; fi; }}
        rename() {{ awk -F'\\t' -v OFS='\\t' -v col=1 -f {params.rename} "$1" -; }}
        {{
            echo "bowtie2 index: {params.index} ($(wc -l < {output.index_chroms}) 条序列)"
            for f in {params.blacklist}; do
                n=$(cat_any $f | grep -v '^#' | cut -f1 | sort -u | comm -12 {output.index_chroms} - | wc -l)
                echo "$f: 共有序列 $n"
                if [ "$n" -eq 0 ]; then
                    echo "错误: $f 与 bowtie2 index 没有共同的染色体名（检查 chr 前缀）" >&2
                    exit 1
                fi
            done
            n=$(cat_any {params.gtf} | grep -v '^#' | cut -f1 | sort -u | rename {output.index_chroms} | wc -l)
            echo "{params.gtf}: 转换为 index 命名后共有序列 $n"
            if [ "$n" -eq 0 ]; then
                echo "错误: GTF 与 bowtie2 index 的染色体名无法对应" >&2
                exit 1
            fi
            if [ -n "{params.fasta}" ]; then
                n=$(rename {output.fasta_chroms} < {output.index_chroms} | wc -l)
                echo "{params.fasta}: 转换为 fasta 命名后共有序列 $n"
                if [ "$n" -eq 0 ]; then
                    echo "错误: fasta 与 bowtie2 index 的染色体名无法对应" >&2
                    exit 1
                fi
            fi
        }} > {output.check}.tmp 2> {log}
        mv {output.check}.tmp {output.check}
        """


def bowtie2_reads(wc, input):
    if is_pe(wc.sample):
        return f"-1 {input.reads[0]} -2 {input.reads[1]} --no-mixed --no-discordant"
    return f"-U {input.reads[0]}"


rule bowtie2:
    input:
        reads=trimmed_reads,
        check="results/qc/reference_check.txt",
    output:
        bam=maybe_temp("results/02_align/{sample}.raw.bam", "raw_bam"),
        bai=maybe_temp("results/02_align/{sample}.raw.bam.bai", "raw_bam"),
        flagstat="results/02_align/stats/{sample}.raw.flagstat",
        idxstats="results/02_align/stats/{sample}.raw.idxstats",
        stats="results/qc/bowtie2/{sample}.log",  # bowtie2 比对统计，summary 与 MultiQC 使用
    params:
        index=REF["bowtie2_index"],
        reads=bowtie2_reads,
        extra=config["bowtie2"]["params"],
        sort_mem=config["samtools_sort_mem"],
    threads: threads_for("align")
    log:
        "<logs>/bowtie2/{sample}.log",
    benchmark:
        "<benchmarks>/bowtie2/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        bowtie2 -p {threads} {params.extra} -x {params.index} {params.reads} 2> {output.stats} \
        | samtools sort -@ {threads} -m {params.sort_mem} -T {output.bam}.tmp -o {output.bam} - 2> {log}
        samtools index -@ {threads} {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        samtools idxstats {output.bam} > {output.idxstats}
        """


rule filter_dedup:
    """
    保留 MAPQ >= filter.mapq 的主比对（-F 1804）；PE 只保留两端都合格的 properly paired；
    filter.remove_mito=true 时去掉线粒体；samtools markdup -r 删除 PCR 重复。
    """
    input:
        bam="results/02_align/{sample}.raw.bam",
        bai="results/02_align/{sample}.raw.bam.bai",
    output:
        bam="results/02_align/{sample}.bam",
        bai="results/02_align/{sample}.bam.bai",
        filtered=temp("results/02_align/{sample}.filtered.bam"),
        flagstat_filtered="results/02_align/stats/{sample}.filtered.flagstat",
        flagstat="results/02_align/stats/{sample}.flagstat",
        markdup="results/02_align/stats/{sample}.markdup.txt",
    params:
        mapq=config["filter"]["mapq"],
        pe=lambda wc: int(is_pe(wc.sample)),
        mito=f"-e {MITO}" if REMOVE_MITO else "",
        sort_mem=config["samtools_sort_mem"],
        tmp="results/02_align/{sample}.sorttmp",
    threads: threads_for("filter")
    log:
        "<logs>/filter_dedup/{sample}.log",
    benchmark:
        "<benchmarks>/filter_dedup/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        chroms=$(samtools idxstats {input.bam} | cut -f1 | grep -v -x -F -e '*' {params.mito} | tr '\\n' ' ')
        if [ {params.pe} -eq 1 ]; then
            # fixmate -r 会把 mate 被过滤掉的 read 标为非 properly paired，再过滤一次去掉孤立 read
            samtools view -u -q {params.mapq} -F 1804 -f 2 {input.bam} $chroms \
            | samtools sort -n -u -@ {threads} -T {params.tmp}.n - \
            | samtools fixmate -r -m -u - - \
            | samtools view -u -F 1804 -f 2 - \
            | samtools sort -@ {threads} -m {params.sort_mem} -T {params.tmp} -o {output.filtered} - 2> {log}
        else
            samtools view -b -q {params.mapq} -F 1804 -o {output.filtered} {input.bam} $chroms 2> {log}
        fi
        samtools flagstat {output.filtered} > {output.flagstat_filtered}
        samtools markdup -@ {threads} -r -s -f {output.markdup} {output.filtered} {output.bam} 2>> {log}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """


rule samtools_stats:
    """最终 BAM 的 samtools stats；MultiQC 据此绘制 insert size（片段长度）分布。"""
    input:
        "results/02_align/{sample}.bam",
    output:
        "results/02_align/stats/{sample}.stats",
    threads: threads_for("filter")
    conda:
        ENV_MAIN
    shell:
        "samtools stats -@ {threads} {input} > {output}"


rule contamination:
    """把 trim 后的 reads 比对到其它物种（默认 human），只记录比对率。"""
    input:
        reads=trimmed_reads,
    output:
        "results/qc/contamination/{sample}.log",
    params:
        index=CONTAM_INDEX if CONTAM["enabled"] else "",
        reads=bowtie2_reads,
        extra=config["bowtie2"]["params"],
    threads: threads_for("align")
    benchmark:
        "<benchmarks>/contamination/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        "bowtie2 -p {threads} {params.extra} -x {params.index} {params.reads} -S /dev/null 2> {output}"
