# Chienlab-RNAseq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**chienlab-rnaseq** is a Nextflow pipeline for performing bacterial RNA-Seq analysis.

This pipeline has been heavily inspired by the [BactSeq](https://github.com/adamd3/BactSeq) pipeline.

## Pipeline summary

The pipeline will perform the following steps:

1. Trim adaptors from reads, performs QC, and filters reads ([`Fastp`](https://github.com/OpenGene/fastp))
2. Align reads to reference genome ([`BWA-MEM`](https://github.com/lh3/bwa/))
3. Performs read quantificantion ([`Rsubread`](https://bioconductor.org/packages/release/bioc/html/Rsubread.html))
4. Generates BigWig files for visualization in genome browsers ([`deeptools`](https://deeptools.readthedocs.io/en/develop/)) 
5. Size-factor scaling and gene length (RPKM) scaling of counts (TMM from [`edgeR`](http://bioconductor.org/packages/release/bioc/html/edgeR.html))
6. Principal component analysis (PCA) of normalised expression values
7. Differential gene expression ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) (optional)

## Installation

You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).

```
Usage:
nextflow run baldikacti/chienlab-rnaseq --data_dir [dir] --sample_file [file] --ref_genome [file] --ref_ann [file] -profile conda [other_options]

Mandatory arguments:
  --data_dir [file]               Path to directory containing FastQ files.
  --ref_genome [file]             Path to FASTA file containing reference genome sequence (bwa) or multi-FASTA file containing coding gene sequences (kallisto).
  --ref_ann [file]                Path to GFF file containing reference genome annotation.
  --sample_file [file]            Path to file containing sample information.
  -profile [str]                  Configuration profile to use.
                                  Available: conda

Other options:
  --cont_tabl [file]              Path to tsv file containing contrasts to be performed for differential expression.
  --l2fc_thresh [str]             Absolute log2(FoldChange) threshold for identifying differentially expressed genes. Default = 1.
  --outdir [file]                 The output directory where the results will be saved (Default: './results').
  --p_thresh [str]                Adjusted p-value threshold for identifying differentially expressed genes. Default = 0.05.
  --max_memory ['32.GB']          Maximum available memory in the system
  --max_cpus [int]                Maximum available cpu's in the system
  --max_time ['10.h']             Maximum time requested time for the pipeline
  -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
```

Explanation of parameters:

- `ref_genome`: genome sequence for mapping reads.
- `ref_ann`: annotation of genes/features in the reference genome.
- `sample_file`: TSV file containing sample information (see below)
- `data_dir`: path to directory containing FASTQ files.
- `cont_tabl`: (optional) table of contrasts to be performed for differential expression.
- `p_thresh`: adjusted p-value threshold for identifying differentially expressed genes. Default = 0.05.
- `l2fc_thresh`: absolute log2(FoldChange) threshold for identifying differentially expressed genes. Default = 1.
- `outdir`: the output directory where the results will be saved (Default: `./results`).
- `-resume`: will re-start the pipeline if it has been previously run.

## Required input

- **Genome sequence**: FASTA file containing the genome sequence. Can be retrieved from NCBI.
- **Gene annotation file**: GFF file containing the genome annotation. Can be retrieved from NCBI.
- **Sample file**: TSV file containing sample information. Must contain the following columns:

  - `sample`: sample ID
  - `file1`: name of the first FASTQ file.
  - `file2`: name of the second FASTQ file. (For single-end sequences, leave blank)
  - `group`: grouping factor for differential expression and exploratory plots.
  - `rep_no`: repeat number (if more than one sample per group).
  - `paired`: data are paired-end? (0 = single-end, 1 = paired-end).
  - `strandedness` : Is data stranded? Options: `unstranded`, `forward`, `reverse`.

  **Example:**

  If data are single-end, leave the `file2` column blank.
  
  Sample file can contain a mix of single-end and paired-end, and a mix of stranded and unstranded samples.

**reference.tsv**
  | sample | file1                  | file2 | group             | rep_no | paired | strandedness |
  | :----: | :-------------------:  | :---: | :---------------: | :----: | :----: | :----------: |
  | AS_1   | SRX1607051_T1.fastq.gz |       | Artificial_Sputum |   1    |   0    |    reverse   |
  | AS_2   | SRX1607052_T1.fastq.gz |       | Artificial_Sputum |   2    |   0    |    reverse   |
  | AS_3   | SRX1607053_T1.fastq.gz |       | Artificial_Sputum |   3    |   0    |    reverse   |
  | MB_1   | SRX1607054_T1.fastq.gz |       | Middlebrook       |   1    |   0    |    reverse   |
  | MB_2   | SRX1607055_T1.fastq.gz |       | Middlebrook       |   2    |   0    |    reverse   |
  | MB_3   | SRX1607056_T1.fastq.gz |       | Middlebrook       |   3    |   0    |    reverse   |

**Optional Contrast Table**

**contrast_ref.tsv**
| contrast1         |	contrast2   |
| :---------------: | :---------: |
|	Artificial_Sputum | Middlebrook |

## Output

1. **fastp** directory containing adaptor-trimmed RNA-Seq files and QC results.
2. **read_counts** directory containing:
   1. `ref_gene_df.tsv`: table of genes in the annotation.
   2. `gene_counts.tsv`: raw read counts per gene.
   3. `cpm_counts.tsv`: size factor scaled counts per million (CPM).
   4. `rpkm_counts.tsv`: size factor scaled and gene length-scaled counts, expressed as reads per kilobase per million mapped reads (RPKM).
3. **PCA_samples** directory containing principal component analysis results.
4. **diff_expr** directory containing differential expression results.
5. **bigwig** directory containing BigWig files.
6. **bwa_aln** directory containing BAM files.

## Slurm Example

```bash
#!/usr/bin/bash
#SBATCH --job-name=chienlab-rnaseq-ba   # Job name
#SBATCH --partition=cpu            # Partition (queue) name
#SBATCH --ntasks=24                   # Number of CPUs
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=06:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/chienlab-rnaseq-ba_%j.log   # Standard output and error log

# Load modules

module load nextflow/24.04.3 conda/latest

# Run pipeline

nextflow run baldikacti/chienlab-rnaseq -r v0.1.0 \
    --data_dir /path/to/fastq \
    --sample_file /path/to/reference.tsv \
    --ref_genome /path/to/organism.fasta \
    --ref_ann /path/to/annotation.gff \
    --cont_tabl /path/to/contrast_ref.tsv \
    --outdir /path/to/results \
    -profile conda \
    -resume

```