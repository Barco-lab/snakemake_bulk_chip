# snakemake_bulk_chip

[中文](README.md) | **English**

Local pipeline for ChIP-seq (and similar data such as CUT&Tag), Snakemake ≥ 9.12.
For ATAC-seq use `pipeline/snakemake_bulk_atac`.

![workflow](docs/images/bulk_chip_workflow.svg)

<sub>Solid boxes are the main workflow, dashed boxes are optional or per-sample conditional steps, green boxes are QC, and ★ marks deliverables (paths relative to `results/`). The diagram labels are in Chinese. PNG version: `docs/images/bulk_chip_workflow.png`.</sub>

<details>
<summary>Text version</summary>

```
fastq ─► trim_galore ─► bowtie2 ─► filtering (MAPQ≥30, properly paired, chrM) ─► dedup ─► macs3 ─► blacklist removal ─► FRiP
             └─ FastQC (raw / trimmed)                                       ├─► bigwig ─► TSS profile / enrichment
                                                                             ├─► fingerprint, samtools stats (fragment length)
                                                                             └─► HOMER motif (optional)
                                     Optional: contamination check (alignment to another species)
                                     Replicates (IP samples): peak reproducibility, consensus-peak correlation and PCA; sex (fraction of reads on chrX / chrY, controls included)
                                     Summary: summary.tsv + QC report + MultiQC
```

</details>

---

## Quick start

```bash
P=/home/LAB/lab_ab_1/pipeline/snakemake_bulk_chip
cd /path/to/project

# 1. Create a project config template (generated when the -c file does not exist); edit genome, etc.
bash $P/run.sh -c config.yaml

# 2. Scan the fastq directory, write <outdir>/samples.tsv and preview the jobs
bash $P/run.sh -c config.yaml -i /path/to/rawdata -o chipseq_out -n

# 3. Check chipseq_out/samples.tsv, then run in the background
nohup bash $P/run.sh -c config.yaml -o chipseq_out -j 48 -y > run.log 2>&1 &
```

Fill in control / peak_type in samples.tsv and check group.

- No environment needs to be activated: `run.sh` calls snakemake from `~/.conda/envs/omics` (use `SMK_ENV=<env path>` for another environment); macs3 and deeptools also use the `omics` environment;
  the replicate figures are drawn with R (ggplot2) from `in202506`.
- There is one copy of the pipeline code shared by all projects; **each project passes its own config with `-c` (required)**.
  The pipeline reads the default `config/config.yaml` first and then overrides it with the `-c` file (nested keys are merged level by level), so a project config only needs the parameters that differ from the defaults.
  The merged configuration actually used by each run is archived in `<outdir>/logs/run_<time>/config.yaml`.
  The default config is the base for every project, so normally leave it unchanged; editing it affects resumed runs of all projects.
- Results and logs are written to the `-o` directory, so several projects can run at the same time.
- After an interruption or error, **rerun the same command to resume**.

`run.sh` options are the same as in the RNA pipeline: `-c` project config (required; a template is generated if missing), `-i` fastq directory, `-o` output directory (default `./chipseq_out`), `-s` existing sample sheet, `-j` cores, `-n` preview, `-y` run without asking; everything after `--` is passed to snakemake.

---

## Input: samples.tsv

```
sample        fq1                      fq2                      group     control     peak_type
H3K4me3_rep1  /data/K4_rep1_1.fq.gz    /data/K4_rep1_2.fq.gz    H3K4me3   input_rep1
H3K4me3_rep2  /data/K4_rep2_1.fq.gz    /data/K4_rep2_2.fq.gz    H3K4me3   input_rep2
H3K27me3_rep1 /data/K27_rep1_1.fq.gz   /data/K27_rep1_2.fq.gz   H3K27me3  input_rep1  broad
input_rep1    /data/input_rep1_1.fq.gz /data/input_rep1_2.fq.gz input
input_rep2    /data/input_rep2_1.fq.gz /data/input_rep2_2.fq.gz input
```

| Column | Description |
|---|---|
| `fq2` | Empty means single-end. Recognises `_1/_2`, `_R1/_R2`, `_R1_001/_R2_001` |
| `control` | input / IgG control of the ChIP, given as **the name of another sample in the sheet**.<br>Samples used as a control are only aligned, QC'd and turned into bigwigs; no peaks are called for them.<br>Leave empty to run without a control (using the old `--nolambda` options). |
| `peak_type` | `narrow` or `broad`; empty uses `peak_type` from the config (default narrow). Must be the same for the IP samples of a group |
| `group` | Optional group a replicate belongs to, used for the replicate analysis; when the column is missing or empty the group is inferred by removing a trailing `_repN` from the sample name (`run.sh -i` fills it in; edit as needed).<br>Only IP samples are compared; the group of a control sample is used for the sex check only |
| `merge` | Optional, whether the replicate goes into the group bigwig: empty or `yes` to include, `no` to leave it out (e.g. a replicate flagged by the replicate QC). Only affects the group bigwig, not the reproducible / consensus peaks or the replicate analysis |

---

## Main settings

| Step | Setting |
|---|---|
| trim | trim_galore (adapter auto-detection) |
| bowtie2 | `-N 1`; plus `--no-mixed --no-discordant` for PE |
| mitochondria | kept by default (removed with `filter.remove_mito: true`); summary reports `mito_rate (chrM) (% of mapped)` |
| macs3 | with control: `-c control --keep-dup all`<br>without control: `--nomodel --nolambda --extsize 200 --min-length 100`<br>broad: `--broad --broad-cutoff 0.1`; `-f BAMPE` for PE |
| QC | FRiP, fingerprint, TSS profile / enrichment, meta gene profile (TSS to TES), fragment length distribution |
| replicates | IP samples only. Pairwise Jaccard of replicate peaks and Venn diagrams of the bases covered by each replicate combination (bedtools multiinter, groups of 2–4 replicates), and between groups using each group's reproducible peaks (`replicates.group_venn` sets the groups compared); within-group consensus sets from DiffBind (reproducible peaks present in ≥ 2 replicates), merged across groups into consensus peaks and counted with `dba.count` on IP reads (input not subtracted; summit ± 200 bp for narrow, merged widths when broad groups are present), log2(CPM + 1) for correlation (heatmap and scatter-plot matrix), clustering and PCA; sex of every sample (controls included) inferred from the fraction of mapped reads on chrX and chrY (called from the chrY region; chrX only compared within groups) |

Other settings:
- filtering keeps MAPQ ≥ 30 and removes secondary, supplementary and unmapped reads; PE keeps properly paired reads with both mates passing;
- PCR duplicates are removed with `samtools markdup -r`;
- ENCODE blacklist regions are removed from the peaks;
- bigwigs use RPGC normalisation and 10 bp bins by default, with blacklist regions excluded.

---

## Key config.yaml settings

| Key | Description |
|---|---|
| `genome` | `mm10` / `hg38`; each has `bowtie2_index`, `fasta` (for HOMER), `gtf` (TSSs and gene bodies of protein_coding genes), `blacklist`, `effective_genome_size` |
| `contamination` | Contamination check: `enabled: true`, `genome` defaults to `hg38`. **There is currently no human bowtie2 index on this machine**; build one and fill in `references.hg38.bowtie2_index` before enabling |
| `filter.mapq` | Default 30 |
| `filter.remove_mito` | Default false |
| `macs3.qvalue` | Default 0.01 |
| `bigwig.normalize` | One or more methods (a list), one `{sample}.{method}.bw` each; the first is used for the QC figures: `RPGC` (default) / `CPM` / `RPKM` / `BPM` / `None` (`bamCoverage --normalizeUsing`); `deseq2_peaks` (`.SFpeaks.bw`, consensus peak counts) / `deseq2_background` (`.SFbg.bw`, genome-wide bins of `bigwig.background_bin`, default 10 kb) use DESeq2 size factors (median-of-ratios, assuming most regions are unchanged between samples; IP samples only), coverage divided by the size factor; values in `04_bigwig/size_factors.{SFpeaks,SFbg}.tsv` and summary.tsv. `bigwig.binsize` defaults to 10 |
| `replicates.group_venn` | Groups compared in the between-group Venn diagrams: each item is 2–4 groups, e.g. `[[CBP_KA, CBP_NE, CBP_Saline], [p300_KA, p300_NE, p300_Saline]]`, one panel per item; when empty, all groups are compared if there are 2–4, otherwise none |
| `bigwig.group` | One bigwig per group (IP groups only, controls are not merged), on by default; `method: merge` (default) merges the group's dedup BAMs and runs `bamCoverage`, `mean` averages the replicate bigwigs with `bigwigAverage` (equal weight per replicate). Only groups with ≥ 2 replicates left after `merge: no` are made.<br>You can first run with `enabled: false`, check the replicate QC, fill in the `merge` column, set `true` and run again: only the merging step is added. Changing the `merge` column later only redoes the affected groups |
| `motif.enabled` | HOMER, only for non-control samples with narrow peaks; `motif.params` defaults to `"-nomotif"` (known motifs only); set it to `""` to add de novo motifs (lengths 8/10/12, much slower) |
| `keep.raw_bam` | Default false: only the filtered, deduplicated BAM is kept |
| `replicates` | Replicate consistency: `enabled`; `min_reps` (replicates a reproducible peak must appear in, default 2); `summits` (DiffBind counts the summit ± this many bp, default 200; not used when broad groups are present); `method` (`pearson` / `spearman`); `pca_top` (default 500); `sex_check` and `y_frac` (default `[0.002, 0.01]`: reads in the chrY region as a percentage of all mapped reads < 0.002% female, > 0.01% male) |
| `references.<genome>.y_region` | chrY region counted for the sex check, i.e. chrY without the pseudoautosomal region (PAR, identical on X and Y): `chrY:1-90000000` for mm10 (the PAR is chrY:90745845-91644698; females also have many reads just upstream of it, so a further ~745 kb is left out), between PAR1 and PAR2 for hg38 (`chrY:2781480-56887902`). The region still contains X-homologous ampliconic sequence, so female samples carry a few reads and the call uses thresholds |

Chromosome naming of the references: the bowtie2 mm10 index and the blacklist use UCSC names (`chr1`); the fasta / gtf come from the same source as the RNA pipeline
(Ensembl release-100 under `/data/genome`, named `1` / `MT`), and names are converted automatically when building the TSS set and the HOMER input.
At start-up the pipeline checks that the names can be matched and writes the result to `results/qc/reference_check.txt`.
The source and download record of the blacklist are in `~/project/reference/blacklist/download.log`.

---

## Output

```
project/
├── samples.tsv
├── results/
│   ├── qc/report.html                   # ★ QC report (single file): overview, workflow & parameters, QC summary, replicate consistency, signal profiles; definitions and references in workflow/report/qc_metrics.yaml
│   ├── qc/summary.tsv                   # reads, mapping rate, duplicate rate, peaks, FRiP, etc. per sample (plain numbers)
│   ├── qc/fragment_length.tsv           # fragments and share per fragment length for PE samples; plot in qc/figures/fragment_length.(pdf / png / svg)
│   ├── 03_peaks/{sample}.narrowPeak     # ★ peaks with blacklist regions removed (.broadPeak for broad samples)
│   ├── 03_peaks/macs3/{sample}/         # raw macs3 output (xls, summits, etc.)
│   ├── 03_peaks/reproducible/{group}.bed  # peaks present in at least min_reps replicates (column 4: number of replicates)
│   ├── 03_peaks/consensus.bed           # union of the reproducible peaks of all groups (DiffBind; summit ± summits bp for narrow)
│   ├── 04_bigwig/{sample}.{norm}.bw     # ★ {norm} is the normalization (RPGC, CPM, SFpeaks, SFbg, ...)
│   ├── 04_bigwig/groups/{group}.{norm}.bw  # ★ merged bigwig per group (for size factors: mean of the members' normalized signal)
│   ├── 04_bigwig/size_factors.{SFpeaks,SFbg}.tsv  # DESeq2 size factors (when used)
│   ├── 05_motif/{sample}/               # HOMER (optional)
│   ├── 02_align/{sample}.dedup.bam      # filtered, deduplicated BAM ({sample}.raw.bam is the raw bowtie2 alignment); stats/ holds flagstat, idxstats, markdup, samtools stats
│   ├── 01_trim/reports/                 # trim_galore reports
│   ├── 00_raw/                          # symlinks to the raw fastq files (renamed to {sample}_R1/R2)
│   └── qc/
│       ├── multiqc/multiqc_report.html  # ★ all QC in one report
│       ├── bowtie2/{sample}.log         # alignment statistics
│       ├── reference_check.txt          # chromosome naming check of index, blacklist and GTF
│       ├── replicates/                  # replicate consistency: samples.tsv, groups.tsv, peak_overlap.tsv, peak_venn.tsv, group_venn.tsv, correlation.tsv, log2cpm.tsv, pca.tsv,
│       │                                #   peak_counts.tsv (consensus peak count matrix), diffbind.rds (DBA object for differential analysis); figures/ holds PDF, PNG (300 dpi) and SVG
│       ├── deeptools/                   # fingerprint, tss_profile, metagene_samples / metagene_groups (meta gene profile; PNG, TSV, matrix)
│       └── fastqc_*/  frip/  contamination/
└── logs/
    └── run_<time>/                    # one per snakemake run; a resumed run only contains the jobs it actually ran
        ├── <step>/<sample>.log        # look here first when something fails
        ├── benchmarks.tsv             # runtime and memory of all jobs in this run (one table)
        ├── snakemake.log              # full snakemake output
        ├── config.yaml, samples.tsv   # config and sample sheet used in this run
        ├── pipeline.txt               # pipeline path, commit and command
        └── software_versions.tsv      # versions and paths of all tools
```

**Columns of summary.tsv** (reads are counted as single reads, R1 + R2 for PE; each rate uses the reads entering that step as the denominator, in %, matching ENCODE and other published standards; the header shows the config value and the denominator; reference values come from the ENCODE ChIP-seq standards, see `results/qc/report.html`):

| Column | Meaning | Reference |
|---|---|---|
| `raw_reads` / `clean_reads` | reads before trimming / entering bowtie2 | |
| `mapped_reads` | mapped reads | |
| `filtered_reads (q30)` | after MAPQ (`filter.mapq`) and pairing filters (chrM also removed with `filter.remove_mito: true`) | |
| `dedup_reads` | after deduplication, used for peaks and bigwig; divided by 2 for PE gives usable fragments | TF / narrow ≥ 20M, broad ≥ 45M fragments |
| `clean_rate (% of raw)` | clean_reads / raw_reads | |
| `mapping_rate (% of clean)` | mapped_reads / clean_reads | compare within the batch |
| `dedup_rate (% of filtered)` | dedup_reads / filtered_reads | |
| `duplicates_rate (% of filtered)` | (filtered_reads − dedup_reads) / filtered_reads | compare within the batch; usually higher for CUT&Tag |
| `mito_rate (chrM) (% of mapped)` | mitochondrial reads / mapped_reads (`filter.mito_name`) | compare within the batch |
| `peaks` / `frip` | number of peaks and FRiP (reads in peaks / dedup_reads) | FRiP ≥ 1% for TF ChIP (ENCODE guidelines) |
| `tss_enrichment (±3000bp)` | TSS centre signal / flank background (estimated from the deeptools profile; window set by `qc.tss_matrix`) |  |
| `contamination_rate (% of clean)` | fraction of clean reads aligning to the contaminant genome (only with `contamination.enabled`) | |
| `subnucleosomal_rate (% of fragments < 120bp)` | share of deduplicated PE fragments shorter than 120 bp (sub-nucleosomal; samtools stats; empty for SE). In CUT&Tag a high share close to the IgG control points to an open-chromatin (Tn5) bias; distribution in `qc/figures/fragment_length.*` | |

---

## Common tasks

```bash
bash $P/run.sh -c config.yaml -o chipseq_out -y                          # resume
bash $P/run.sh -c config.yaml -o chipseq_out -y -- --until macs3         # run up to peak calling
bash $P/run.sh -c config.yaml -o chipseq_out -y -- --forcerun macs3      # rerun macs3 (and downstream) after changing its options
bash $P/run.sh -c config.yaml -o chipseq_out -y -- --unlock              # unlock after a killed run
```

To add samples, append rows to samples.tsv and run; finished samples are not rerun.

---

## Main changes from the old version

- Runs locally instead of submitting PBS jobs; parameters are no longer passed with `--config reads=<JSON>` but through samples.tsv.
- SE/PE rules are merged and file names no longer carry `_pe` / `_se`.
- ATAC-seq has been split into the separate `snakemake_bulk_atac` pipeline.
- New: input/IgG controls, broad peaks, blacklist filtering, FRiP, fingerprint, TSS enrichment, fragment length distribution, a QC report and an optional contamination check.
- Unique alignment filtering changed from `grep -v "XS:i:"` (which left orphan mates in PE data) to MAPQ ≥ 30 plus properly paired.
- Deduplication changed from sambamba to `samtools markdup` (sambamba is not installed here); peaks are called with macs3 from the `omics` environment.
- HOMER now receives the genome FASTA (HOMER on this machine has no genomes installed) and analyses summit ± 100 bp.
- Human alignment became an optional contamination check, off by default, and another species can be chosen.
- Stand-alone scripts that are not part of the pipeline (`idr_peak.sh`, `mean.py`, `download_sra.sh`, etc.) were removed; download data with `pipeline/utilities/download_sra_iseq.sh`.
