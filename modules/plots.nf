process PCA_SAMPLES {
    tag "$tmm_counts"
    label 'process_medium'
    publishDir "${params.outdir}/PCA_samples", mode: 'copy'
    conda "conda-forge::r-base=4.1 conda-forge::r-optparse conda-forge::r-tibble conda-forge::r-rcolorbrewer conda-forge::r-ggplot2 bioconda::r-ggbiplot"

    input:
    path tmm_counts
    path meta_merged

    output:
    path '*.{rds,png}', emit: pca_out

    script:
    """
    pca.R -o ./
    """
}
