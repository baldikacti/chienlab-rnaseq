process BAM2BIGWIG {
    tag "$bam"
    label 'process_medium'
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    conda "bioconda::deeptools=3.5.3"

    input:
    path bam
    path bai

    output:
    path "${bam.baseName}.bw"

    script:

    """
    bamCoverage -b ${bam} -o ${bam.baseName}.bw --binSize 5 -p ${task.cpus}
    """  
}
