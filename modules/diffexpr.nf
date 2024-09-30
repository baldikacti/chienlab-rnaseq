process DIFF_EXPRESSION {
    tag "$gene_counts"
    label 'process_high'
    publishDir "${params.outdir}/diff_expr", mode: 'copy'
    conda "conda-forge::r-base=4.1 conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-tibble conda-forge::r-rsqlite conda-forge::r-plyr bioconda::bioconductor-deseq2 bioconda::bioconductor-enhancedvolcano"

    input:
    path gene_counts
    path meta_merged
    path cont_tabl
    val p_thresh
    val l2fc_thresh

    output:
    path '*.tsv', emit: deseq_res
    path '*.png', emit: deseq_volcano

    script:
    """
    [ ! -f contrast_table.tsv ] && ln -s $cont_tabl contrast_table.tsv
    diffexpr.R -p $p_thresh -l $l2fc_thresh -o ./ 
    """
}
