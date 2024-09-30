process MAKE_META_FILE {
    tag "$sample_file"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    conda "conda-forge::r-base=4.1 conda-forge::r-optparse"

    input:
    path sample_file

    output:
    path 'sample_metadata.tsv', emit: sample_metadata

    script:
    """
    make_meta_file.R --sample_file $sample_file --data_dir ${params.data_dir} --outf sample_metadata.tsv
    """
}
