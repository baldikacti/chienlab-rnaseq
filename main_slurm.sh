#!/usr/bin/bash
#SBATCH --job-name=chienlab-rnaseq-ba   # Job name
#SBATCH --partition=cpu            # Partition (queue) name
#SBATCH --ntasks=8                   # Number of CPU cores
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --mem=32gb                     # Job memory request
#SBATCH --time=06:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/chienlab-rnaseq-ba_%j.log   # Standard output and error log
date;hostname;pwd

# Load modules

module load nextflow/23.04.1 miniconda/22.11.1-1

# Run pipeline

nextflow run main_dev.nf \
    --data_dir data/test/raw \
    --sample_file data/test/reference.tsv \
    --ref_genome references/NC_011916.fasta \
    --ref_ann references/ccna.gff \
    --outdir results/test \
    -profile conda \
    -resume

date