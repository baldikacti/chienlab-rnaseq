process BAM2BIGWIG {
    tag "$bam"
    label 'process_medium'
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    conda "envs/align_map.yml"

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
