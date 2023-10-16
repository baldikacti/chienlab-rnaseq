#!/usr/bin/zsh

# Run pipeline

nextflow run dev.nf \
    --data_dir data/QUO1006522/raw \
    --sample_file data/QUO1006522/reference.tsv \
    --cont_tabl data/QUO1006522/contrast_ref.tsv \
    --ref_genome references/NC_011916.fasta \
    --ref_ann references/ccna.gff \
    --outdir results/test2 \
    -profile conda \
    -resume