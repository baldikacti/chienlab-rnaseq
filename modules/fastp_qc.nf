process FASTP {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy',
        saveAs: { filename ->
                      if (filename.endsWith('.html')) "fastp_report/$filename"
                      else if (filename.endsWith('failed_reads.fq.gz')) "failed_reads/$filename"
                      else params.save_trimmed ? filename : null
                }
    conda "bioconda::fastp=0.23.4"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fq.gz"),      emit: trimmed_reads
    path("*fastp.html"),                           emit: fastp_qc_reports
    tuple val(meta), path("*_failed_reads.fq.gz"), emit: failed_reads

    script:

    def name = task.ext.prefix ?: "${meta.sample_id}"

   if (meta.paired_end) {
        """
        fastp -i ${reads[0]} -I ${reads[1]} -o ${name}_1_trimmed.fq.gz -O ${name}_2_trimmed.fq.gz --failed_out ${name}_failed_reads.fq.gz \\
            --thread ${task.cpus} --html ${name}_fastp.html
        """
    } else {
        """
        fastp -i ${reads} -o ${name}_trimmed.fq.gz --failed_out ${name}_failed_reads.fq.gz --thread ${task.cpus} --html ${name}_fastp.html
        """
    }
}


process MULTIQC {
    label 'process_single'
    publishDir "${params.outdir}/fastp/multiqc", mode:'copy'

    input:
    path('*')

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc .
    """
}