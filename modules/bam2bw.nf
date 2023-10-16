process BAM2BIGWIG {
    tag "$bam"
    label 'process_medium'
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    conda "envs/align_map.yml"

    input:
    path bam
    path bai

    output:
    path '*.bw', emit: bigwig

    script:

    """
    bamCoverage -b ${bam} -o ${bam}.bw
    """  
}
