process MAKE_BWA_INDEX {
    tag "$ref_genome"
    label 'process_medium'
    publishDir "${params.outdir}/bwa_idx", mode: 'copy'
    conda "bioconda::bwa=0.7.17"

    input:
    path ref_genome

    output:
    path '*', emit: bwa_idx

    script:
    """
    bwa index -p ref_idx $ref_genome
    """
}

process BWA_ALIGN {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir "${params.outdir}/bwa_aln", mode: 'copy'
    conda "bioconda::bwa=0.7.17 bioconda::samtools=1.15"

    input:
    tuple val(meta), path(trimmed_reads)
    path idx

    output:
    path '*.bam', emit: bam_files
    path '*.bai', emit: bai_files
    path '*.counts', emit: count_files

    script:

    def name = task.ext.prefix ?: "${meta.sample_id}"

    if (meta.paired_end) {
        """
        bwa mem -t ${task.cpus} ref_idx \\
            ${trimmed_reads[0]} ${trimmed_reads[1]} | \\
            samtools sort -@ ${task.cpus - 1} -O bam - > ${name}.bam
        samtools index -@ ${task.cpus} ${name}.bam
        samtools view -F 0x4 ${name}.bam | cut -f 1 | sort | uniq | wc -l > ${name}.counts
        """
    } else {
        """
        bwa mem -t ${task.cpus} ref_idx ${trimmed_reads} \\
            | samtools sort -@ ${task.cpus - 1} -O bam - > ${name}.bam
        samtools index -@ ${task.cpus} ${name}.bam
        samtools view -F 0x4 ${name}.bam | cut -f 1 | sort | uniq | wc -l > ${name}.counts
        """
    }
}

process COUNT_READS {
    tag "$gff"
    label 'process_medium'
    publishDir "${params.outdir}/read_counts", mode: 'copy'
    conda "conda-forge::r-base=4.1 conda-forge::r-optparse conda-forge::r-ape conda-forge::r-stringr conda-forge::r-ggplot2 conda-forge::r-scales conda-forge::r-rcolorbrewer conda-forge::r-reshape2 conda-forge::r-tibble conda-forge::r-rsqlite"

    input:
    path bam
    path bai
    path counts
    path meta
    path gff

    output:
    path 'gene_counts.tsv', emit: counts_df
    path 'gene_counts_pc.tsv', emit: counts_df_pc
    path 'counts_summary.tsv', emit: counts_summary
    path 'ref_gene_df.tsv', emit: ref_gene_df
    path 'library_composition.png', emit: libcomp_plot
    path 'library_composition_proportions.png', emit: libcomp_plot_prop

    script:

    """
    count_reads.R -m $meta -g $gff -t ${task.cpus}
    """
}
