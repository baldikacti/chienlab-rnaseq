#!/usr/bin/zsh

# Run pipeline

nextflow run main_dev.nf \
    --data_dir data/test/raw \
    --sample_file data/test/reference.tsv \
    --ref_genome references/NC_011916.fasta \
    --ref_ann references/ccna.gff \
    --outdir results/test \
    -profile conda \
    -resume