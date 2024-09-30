#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
================================================================================
    Modules
================================================================================
*/

include {MAKE_META_FILE} from './modules/metadata'
include {FASTP; MULTIQC} from './modules/fastp_qc'
include {MAKE_BWA_INDEX; BWA_ALIGN; COUNT_READS} from './modules/bwa_align'
include {BAM2BIGWIG} from './modules/bam2bw'
include {TMM_NORMALISE_COUNTS} from './modules/normalisation'
include {PCA_SAMPLES} from './modules/plots'
include {DIFF_EXPRESSION} from './modules/diffexpr'

/*
================================================================================
    Validate inputs and create channels for files
================================================================================
*/

// required inputs
if (params.sample_file) {
    ch_samples = file(params.sample_file, checkIfExists: true)
} else { exit 1, 'Sample file not specified!' }

if (params.ref_genome) {
    ch_fasta_file = file(params.ref_genome, checkIfExists: true)
} else { exit 1, 'Reference genome FASTA file not specified!' }

if (params.ref_ann) {
    ch_gff_file = file(params.ref_ann, checkIfExists: true)
} else { exit 1, 'Reference genome GFF file not specified!' }

// optional inputs
if (params.cont_tabl) {
    ch_cont_file = file(params.cont_tabl, checkIfExists: true)
}

/*
================================================================================
    Functions
================================================================================
*/
def create_fastq_channel(LinkedHashMap row) {
    // create sample metadata
    def meta = [:]
    meta.sample_id    = row.sample
    meta.paired_end   = row.paired.toBoolean()
   
    // add path(s) of the fastq file(s) to the metadata
    def fastq_meta = []
    if (!file(row.file1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.file1}"
    }
    if (meta.paired_end) {
        if (!file(row.file2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.file2}"
        }
        fastq_meta = [ meta, [ file(row.file1), file(row.file2) ] ]
    } else {
        fastq_meta = [ meta, [ file(row.file1) ] ]
    }
    return fastq_meta
}


workflow {
 
    /*
     * Make metadata file linking samples with FastQ files
     */
    MAKE_META_FILE (
        ch_samples
    )
    ch_metadata = MAKE_META_FILE.out.sample_metadata


    /*
     *  Create channels for input files
     */

    ch_metadata
        .splitCsv(header: true, sep:'\t')
        .map { create_fastq_channel(it) }
        .set { ch_raw_reads_fastp }
    
   /*
     *  QC and trimming for fastq files
     */
    
    FASTP (
        ch_raw_reads_fastp
    )
    ch_trimmed_reads = FASTP.out.trimmed_reads
    ch_fastq_reports = FASTP.out.fastp_qc_reports.collect()

    //MULTIQC (
    //    ch_fastq_reports
    //)

    /*
     *  Align / pseudo-align reads
     */
    MAKE_BWA_INDEX (
        ch_fasta_file
    )
    ch_bwa_idx = MAKE_BWA_INDEX.out.bwa_idx

    BWA_ALIGN (
        ch_trimmed_reads,
        ch_bwa_idx
    )
    ch_bwa_out_bam = BWA_ALIGN.out.bam_files.collect()
    ch_bwa_out_bai = BWA_ALIGN.out.bai_files.collect()
    ch_bwa_out_count = BWA_ALIGN.out.count_files.collect()

    /*
     *  Convert BAM files to BigWig format
     */
    BAM2BIGWIG (
        ch_bwa_out_bam.flatten(),
        ch_bwa_out_bai.flatten()
    )

    /*
     *  Count reads mapped per gene; summarise library composition
     */
    COUNT_READS (
        ch_bwa_out_bam,
        ch_bwa_out_bai,
        ch_bwa_out_count,
        ch_metadata,
        ch_gff_file
    )
    ch_readcounts_df = COUNT_READS.out.counts_df
    ch_readcounts_df_pc = COUNT_READS.out.counts_df_pc
    ch_refgene_df = COUNT_READS.out.ref_gene_df


    /*
     *  Get normalised read counts per gene
     */
    TMM_NORMALISE_COUNTS (
        ch_readcounts_df,
        ch_refgene_df
    )
    ch_cpm_counts = TMM_NORMALISE_COUNTS.out.cpm_counts
    ch_rpkm_counts = TMM_NORMALISE_COUNTS.out.rpkm_counts
    // NB the resulting counts are log-transformed by default
    
    /*
     *  Principal component analysis (PCA) of samples
     */
    PCA_SAMPLES (
        ch_cpm_counts,
        ch_metadata
    )
    ch_pca_out = PCA_SAMPLES.out.pca_out
    
    /*
     *  Differential gene expression (DESeq2)
     */
    if (params.cont_tabl) {
        DIFF_EXPRESSION (
            ch_readcounts_df,
            ch_metadata,
            ch_cont_file,
            params.p_thresh,
            params.l2fc_thresh
        )
        ch_deseq_res = DIFF_EXPRESSION.out.deseq_res.collect()
    }

}