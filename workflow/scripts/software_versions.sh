#!/usr/bin/env bash
# 输出流程用到的软件版本与实际路径（流程启动时由 onstart 调用）
#   software_versions.sh <main 环境> <macs3/deeptools 环境> [fastqc 的 --java 参数]
#   环境变量 SMK_PATH：正在运行的 snakemake
main=$1 macs=$2 java=${3:-}
export PATH=$main/bin:$PATH

v() { printf "%s\t%s\t%s\n" "$1" "$2" "$(command -v "$3" || echo "$3")"; }
printf "tool\tversion\tpath\n"
v snakemake   "$("$SMK_PATH" --version)" "$SMK_PATH"
v trim_galore "$(trim_galore --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)" trim_galore
v cutadapt    "$(cutadapt --version)" cutadapt
v fastqc      "$(fastqc $java --version 2>/dev/null | grep -oE 'v[0-9.]+')" fastqc
v bowtie2     "$(bowtie2 --version | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')" bowtie2
v samtools    "$(samtools --version | head -1 | cut -d' ' -f2)" samtools
v multiqc     "$(multiqc --version | grep -oE '[0-9.]+$')" multiqc
v macs3       "$("$macs/bin/macs3" --version | grep -oE '[0-9.]+$')" "$macs/bin/macs3"
v deeptools   "$("$macs/bin/bamCoverage" --version | grep -oE '[0-9.]+$')" "$macs/bin/bamCoverage"
v bedtools    "$("$macs/bin/bedtools" --version | grep -oE '[0-9.]+$')" "$macs/bin/bedtools"
