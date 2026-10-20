# snakemake_bulk_chip

ChIP-seq（以及 CUT&Tag 等类似数据）本地分析流程（Snakemake ≥ 9.12）。
ATAC-seq 请用 `pipeline/snakemake_bulk_atac`。

![workflow](docs/images/bulk_chip_workflow.svg)

<sub>实线框为主流程，虚线框为可选或按样本条件执行的步骤，绿色为 QC，★ 为交付结果（路径相对 `results/`）。PNG 版本：`docs/images/bulk_chip_workflow.png`。</sub>

<details>
<summary>文字版流程</summary>

```
fastq ─► trim_galore ─► bowtie2 ─► 过滤（MAPQ≥30、properly paired、chrM）─► 去重 ─► macs3 ─► 去除 blacklist ─► FRiP
             └─ FastQC (raw / trimmed)                                       ├─► bigwig ─► TSS profile / enrichment
                                                                             ├─► fingerprint、samtools stats（片段长度）
                                                                             └─► HOMER motif（可选）
                                     可选：污染检查（比对到其它物种）
                                     汇总：summary.tsv + MultiQC 报告
```

</details>

---

## 快速开始

```bash
P=/home/LAB/lab_ab_1/pipeline/snakemake_bulk_chip
cd /path/to/project

# 1. 生成项目配置模板（-c 指定的文件不存在时自动生成），修改其中的 genome 等
bash $P/run.sh -c config.yaml

# 2. 扫描 fastq 目录，生成 <outdir>/samples.tsv，并预览任务
bash $P/run.sh -c config.yaml -i /path/to/rawdata -o chipseq_out -n

# 3. 检查 chipseq_out/samples.tsv 后，后台运行
nohup bash $P/run.sh -c config.yaml -o chipseq_out -j 48 -y > run.log 2>&1 &
```

在 samples.tsv 中填写 control / peak_type。

- 不需要激活任何环境，`run.sh` 直接调用 `~/.conda/envs/omics` 中的 snakemake（可用 `SMK_ENV=<环境路径>` 指定其它环境）；macs3 和 deeptools 也使用 `omics` 环境。
- 流程代码只有一份，所有项目共用；**每个项目用 `-c` 指定自己的配置文件（必填）**。
  流程先读默认配置 `config/config.yaml`，再用 `-c` 的文件覆盖（嵌套的键逐层合并），所以项目配置可以只写与默认不同的参数。
  每次运行实际使用的（合并后的）配置会存档到 `<outdir>/logs/run_<时间>/config.yaml`。
  默认配置是所有项目的基础，一般不要改它；改了会影响所有项目的续跑。
- 结果、日志都写在 `-o` 指定的目录，多个项目可以同时跑。
- 中断或出错后，**重新运行同一条命令即可断点续跑**。

`run.sh` 的参数与 RNA 流程相同：`-c` 项目配置（必填，不存在时生成模板）、`-i` fastq 目录、`-o` 输出目录（默认 `./chipseq_out`）、`-s` 已有样本表、`-j` 核数、`-n` 预览、`-y` 直接运行、`--` 之后的参数原样传给 snakemake。

---

## 输入：samples.tsv

```
sample        fq1                      fq2                      control     peak_type
H3K4me3_rep1  /data/K4_rep1_1.fq.gz    /data/K4_rep1_2.fq.gz    input_rep1
H3K27me3_rep1 /data/K27_rep1_1.fq.gz   /data/K27_rep1_2.fq.gz   input_rep1  broad
input_rep1    /data/input_rep1_1.fq.gz /data/input_rep1_2.fq.gz
```

| 列 | 说明 |
|---|---|
| `fq2` | 为空即单端。识别 `_1/_2`、`_R1/_R2`、`_R1_001/_R2_001` |
| `control` | ChIP 的 input / IgG 对照，填写**样本表中另一个样本的名字**。<br>被引用为 control 的样本只做比对、QC 和 bigwig，不单独 call peak。<br>留空则不用对照（使用旧版的 `--nolambda` 参数）。 |
| `peak_type` | `narrow` 或 `broad`，留空则使用 config 中的 `peak_type`（默认 narrow） |

---

## 主要参数

| 步骤 | 设置 |
|---|---|
| trim | trim_galore（自动识别接头） |
| bowtie2 | `-N 1`；PE 加 `--no-mixed --no-discordant` |
| 线粒体 | 默认保留（`filter.remove_mito: true` 时去除），summary 中报告 `mito_rate` |
| macs3 | 有对照：`-c control --keep-dup all`<br>无对照：`--nomodel --nolambda --extsize 200 --min-length 100`<br>broad：`--broad --broad-cutoff 0.1`；PE 用 `-f BAMPE` |
| QC | FRiP、fingerprint、TSS profile / enrichment、片段长度分布 |

其它设置：
- 过滤条件为 MAPQ ≥ 30，并去掉 secondary、supplementary、unmapped reads；PE 只保留两端都合格的 properly paired；
- 用 `samtools markdup -r` 删除 PCR 重复；
- peak 去除 ENCODE blacklist 区域；
- bigwig 默认 RPGC 归一化、bin 10bp，并排除 blacklist 区域。

---

## 配置 config.yaml 要点

| 配置项 | 说明 |
|---|---|
| `genome` | `mm10` / `hg38`；每个物种有 `bowtie2_index`、`fasta`（HOMER 用）、`gtf`（TSS 用）、`blacklist`、`effective_genome_size` |
| `contamination` | 污染检查：`enabled: true`，`genome` 默认 `hg38`。**本机目前没有人的 bowtie2 index**，打开前需先构建并填入 `references.hg38.bowtie2_index` |
| `filter.mapq` | 默认 30 |
| `filter.remove_mito` | 默认 false |
| `macs3.qvalue` | 默认 0.01 |
| `bigwig.normalize` | `RPGC`（默认）/ `CPM` / `RPKM` / `BPM` |
| `motif.enabled` | HOMER，仅对 narrow peak 的非 control 样本；`motif.params: "-nomotif"` 只做已知 motif，速度快很多 |
| `keep.raw_bam` | 默认 false：只保留过滤去重后的 BAM |

参考文件的染色体命名：bowtie2 mm10 index、blacklist 是 UCSC 命名（`chr1`）；fasta / gtf 与 RNA 流程同源
（`/data/genome` 下的 Ensembl release-100，命名为 `1` / `MT`），生成 TSS 和 HOMER 输入时自动互转命名。
流程启动时会自动检查能否对应，结果写在 `results/qc/reference_check.txt`。
blacklist 的来源和下载记录见 `~/reference/blacklist/download.log`。

---

## 输出

```
project/
├── samples.tsv
├── results/
│   ├── tables/summary.tsv               # ★ 每个样本的 reads 数、比对率、重复率、peak 数、FRiP 等
│   ├── 03_peaks/{sample}.narrowPeak     # ★ 已去除 blacklist 的 peak（broad 样本为 .broadPeak）
│   ├── 03_peaks/macs3/{sample}/         # macs3 原始输出（xls、summits 等）
│   ├── 04_bigwig/{sample}.bw            # ★
│   ├── 05_motif/{sample}/               # HOMER（可选）
│   ├── 02_align/{sample}.bam            # 过滤去重后的 BAM；stats/ 下为 flagstat、idxstats、markdup、samtools stats
│   ├── 01_trim/reports/                 # trim_galore 报告
│   ├── 00_raw/                          # 原始 fastq 的软链接（统一命名为 {sample}_R1/R2）
│   └── qc/
│       ├── multiqc/multiqc_report.html  # ★ 所有 QC 汇总
│       ├── bowtie2/{sample}.log         # 比对统计
│       ├── reference_check.txt          # index、blacklist、GTF 染色体命名检查
│       └── fastqc_*/  frip/  deeptools/  contamination/
└── logs/
    └── run_<时间>/                    # 每执行一次 snakemake 新建一个；续跑时只含这次实际跑了的任务
        ├── <步骤>/<样本>.log          # 出错时先看这里
        ├── benchmarks.tsv             # 本次所有任务的耗时与内存（一张表）
        ├── snakemake.log              # snakemake 完整输出
        ├── config.yaml、samples.tsv   # 本次使用的配置与样本表
        ├── pipeline.txt               # 流程路径、commit 版本、运行命令
        └── software_versions.tsv      # 各软件的版本与实际路径
```

**summary.tsv 的各列**（PE 样本按 read pairs 计，修复了旧版 SE 样本也除以 2 的问题）：

| 列 | 含义 | 参考 |
|---|---|---|
| `raw_reads` / `clean_reads` | trim 前 / 进入 bowtie2 的 reads | |
| `mapping_rate` | bowtie2 overall alignment rate | |
| `contamination_rate` | 比对到污染物种的比例（可选） | |
| `mito_rate` | 比对上的 reads 中 chrM 的比例 | |
| `filtered_reads` | MAPQ、配对、chrM 过滤后 | |
| `dup_rate` | PCR 重复比例 | |
| `final_reads` | 去重后，用于 call peak 和 bigwig | |
| `peaks` / `frip` | peak 数与 FRiP | ENCODE：TF ChIP 的 FRiP > 1% |
| `tss_enrichment` | TSS 中心信号 / 两端背景（按 deeptools profile 估算） | |

---

## 常用操作

```bash
bash $P/run.sh -c config.yaml -o chipseq_out -y                          # 断点续跑
bash $P/run.sh -c config.yaml -o chipseq_out -y -- --until macs3         # 只跑到 call peak
bash $P/run.sh -c config.yaml -o chipseq_out -y -- --forcerun macs3      # 改了 macs3 参数后强制重跑（及其下游）
bash $P/run.sh -c config.yaml -o chipseq_out -y -- --unlock              # 强行中断后解锁
```

在 samples.tsv 中追加样本后直接运行即可，已完成的样本不会重跑。

---

## 与旧版的主要区别

- 不再提交 PBS，改为本地运行；不再通过 `--config reads=<JSON>` 传参，改用 samples.tsv。
- SE/PE 规则合并，文件名不再带 `_pe` / `_se`。
- ATAC-seq 拆分为独立的 `snakemake_bulk_atac` 流程。
- 新增 input/IgG 对照、broad peak、blacklist 过滤、FRiP、fingerprint、TSS enrichment、片段长度分布，以及可选的污染检查。
- 唯一比对从 `grep -v "XS:i:"`（PE 会留下孤立 mate）改为 MAPQ ≥ 30 加 properly paired 过滤。
- 去重从 sambamba 改为 `samtools markdup`（本机没有 sambamba）；call peak 使用 `omics` 环境中的 macs3。
- HOMER 改为传入基因组 FASTA（本机 HOMER 没有安装基因组），以 summit ± 100bp 做 motif 分析。
- human 比对改为可选的污染检查，默认关闭，并且可以选择其它物种。
- 删除了不属于流程本身的独立脚本（`idr_peak.sh`、`mean.py`、`download_sra.sh` 等）；下载数据请用 `pipeline/utilities/download_sra_iseq.sh`。
