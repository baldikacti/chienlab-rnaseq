#!/usr/bin/zsh

# Run pipeline

nextflow run main.nf --data_dir data/GRL10/raw --sample_file data/GRL10/reference.tsv --cont_tabl data/GRL10/contrast_ref.tsv --ref_genome references/NC_011916.fasta --ref_ann references/ccna.gff -profile docker --outdir results/GRL10 -resume
