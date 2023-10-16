process BAM2BIGWIG {
    tag "$bam"
    label 'process_medium'
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    conda "envs/align_map.yml"

    input:
    path bam
    path meta

    output:
    path '*.bw', emit: bigwig

    script:

    def name = task.ext.prefix ?: "${meta.sample_id}"

    if (meta.paired_end) {
        """
        samtools view -b -f 128 -F 16 ${name}.bam > ${name}_fwd1.bam
        samtools view -b -f 80 ${name}.bam > ${name}_fwd2.bam
        samtools merge -f ${name}_fwd.bam ${name}_fwd1.bam ${name}_fwd2.bam
        samtools index ${name}_fwd.bam
        bamCoverage -b ${name}_fwd.bam -o ${name}_fwd.bw
        samtools view -b -f 144 ${name}.bam > ${name}_rev1.bam
        samtools view -b -f 64 -F 16 ${name}.bam > ${name}_rev2.bam
        samtools merge -f ${name}_rev.bam ${name}_rev1.bam ${name}_rev2.bam
        samtools index ${name}_rev.bam
        bamCoverage -b ${name}_rev.bam -o ${name}_rev.bw
        """
    } else {
        """
        # Forward strand
        bamCoverage -b ${name}.bam -o ${name}.fwd.bw --samFlagExclude 16

        # Reverse strand
        bamCoverage -b ${name}.bam -o ${name}.rev.bw --samFlagInclude 16
        """
    }   
}
