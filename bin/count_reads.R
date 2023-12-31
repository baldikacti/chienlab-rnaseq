#!/usr/bin/env Rscript
library(optparse)
library(Rsubread)
library(ape)
library(stringr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(reshape2)
library(tibble)


option_list <- list(
    make_option(c("-m", "--metadata"),
        type = "character", default = NULL,
        help = "sample metadata tsv file", metavar = "character"
    ),
    make_option(c("-g", "--gff"),
        type = "character", default = NULL,
        help = "GFF annotation file for the reference strain",
        metavar = "character"
    ),
    make_option(c("-t", "--threads"),
        type = "numeric", default = 1,
        help = "number of threads to use. default = 1.",
        metavar = "numeric"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

meta_f <- opt$metadata
gff_f <- opt$gff
threads <- opt$threads

## ------------------------------------------------------------------------------
## Read data
## ------------------------------------------------------------------------------
meta_tab <- read.table(meta_f, header = TRUE, sep = "\t")
## columns: sample	 file1   file2	group	rep_no    paired  strandedness


## cat the counts files

total_counts_list <- lapply(meta_tab$sample, function(x) {
    mapped_count <- read.table(paste0(x, ".counts"), header = FALSE)
    colnames(mapped_count) <- "mapped"
    mapped_count$sample <- x
    mapped_count
})
merged_total_counts <- as.data.frame(do.call(rbind, total_counts_list))



## ------------------------------------------------------------------------------
## Read genome annotation
## ------------------------------------------------------------------------------
ref_annot <- ape::read.gff(gff_f, na.strings = c(".", "?"), GFF3 = TRUE)

ref_annot <- subset(ref_annot, type == "gene")

gene_attr <- stringr::str_split(ref_annot$attributes, ";")
locus_tags <- unlist(lapply(gene_attr, function(x) {
    x[grepl("locus_tag", x)]
}))
gene_biotypes <- unlist(lapply(gene_attr, function(x) {
    x[grepl("gene_biotype", x)]
}))
common_gene_names <- unlist(lapply(gene_attr, function(x) {
    x <- x[grepl("gene=", x)]
    x[identical(x, character(0))] <- ""
    x
}))
gene_lengths <- (ref_annot$end - ref_annot$start) + 1

ref_gene_df <- data.frame(
    locus_tag = locus_tags,
    biotype = gene_biotypes,
    gene_name = common_gene_names,
    gene_length = gene_lengths
)
ref_gene_df$locus_tag <- gsub("locus_tag=", "", ref_gene_df$locus_tag)
ref_gene_df$biotype <- gsub("gene_biotype=", "", ref_gene_df$biotype)
ref_gene_df$gene_name <- gsub("gene=", "", ref_gene_df$gene_name)

write.table(
    ref_gene_df, "ref_gene_df.tsv",
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)


## ------------------------------------------------------------------------------
## Count reads mapping to genes
## ------------------------------------------------------------------------------
bamfilesCount <- paste0(meta_tab$sample, ".bam")
ispaired <- as.logical(meta_tab$paired)
strand <- sapply(meta_tab$strandedness, function(x) {
    switch(as.character(x),
        "unstranded" = 0,
        "forward" = 1,
        "reverse" = 2,
        stop("Invalid input")
    )
}, USE.NAMES = FALSE)
# 0 (unstranded), 1 (stranded) and 2 (reversely stranded)


gene_counts <- Rsubread::featureCounts(
    bamfilesCount,
    annot.ext = gff_f,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "gene",
    GTF.attrType = "locus_tag",
    nthreads = threads,
    countMultiMappingReads = TRUE,
    fraction = TRUE, ## assign fractional counts to multimappers
    isPairedEnd = ispaired,
    strandSpecific = strand
)
colnames(gene_counts$counts) <- gsub(".bam", "", colnames(gene_counts$counts))
colnames(gene_counts$counts) <- gsub("\\.", "_", colnames(gene_counts$counts))


counts_mat <- gene_counts$counts
counts_mat <- tibble::rownames_to_column(as.data.frame(counts_mat), "feature_id")


write.table(
    counts_mat, "gene_counts.tsv",
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

## protein-coding genes only
gene_counts_pc <- counts_mat[ref_gene_df$biotype == "protein_coding", ]
write.table(
    gene_counts_pc, "gene_counts_pc.tsv",
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)


## ------------------------------------------------------------------------------
## Plot library composition
## ------------------------------------------------------------------------------
## set up plots
brewer_pallette1 <- brewer.pal(9, "Set1")
brewer_pallette3 <- brewer.pal(8, "Dark2")

gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
ggColsDefault <- (gg_color_hue(4))
ggCols <- brewer_pallette1[c(1, 3, 4, 5, 2, 7, 8)]

## summarise counts per sample

all_biotypes <- unique(ref_gene_df$biotype)

biotype_counts <- data.frame(do.call(
    cbind,
    lapply(all_biotypes, function(biotype) {
        colSums(gene_counts$counts[ref_gene_df$biotype == biotype, , drop = FALSE])
    })
))
colnames(biotype_counts) <- all_biotypes



counts_summary <- data.frame(
    sample = meta_tab$"sample",
    group = meta_tab$"group",
    rep = meta_tab$"rep_no")

counts_summary$rRNA <- biotype_counts$rRNA

## add total mapped counts
counts_summary <- merge(counts_summary, merged_total_counts, by = "sample")

counts_summary$other <- counts_summary$mapped - counts_summary$rRNA

write.table(
    counts_summary, "counts_summary.tsv",
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)



counts_summary <- counts_summary[rev(order(counts_summary$sample)), ]


cc1 <- 12

all_biotypes <- c(all_biotypes, "other")
non_rRNA_btypes <- all_biotypes[!all_biotypes == "rRNA"]


#############################
## raw counts plot
#############################
counts_melt <- reshape2::melt(
    counts_summary,
    id.vars = c("sample"),
    measure.vars = c(
        "other",
        "rRNA"
        # all_biotypes
    )
)
counts_melt$sample <- factor(
    counts_melt$sample,
    levels = rev(unique(sort(counts_melt$sample)))
)
counts_melt$variable <- factor(counts_melt$variable, levels = c(
    "rRNA",
    "other"
    # non_rRNA_btypes
))


ylabel <- ifelse(isTRUE(ispaired), "Million read pairs", "Million reads")

p1 <- ggplot(
    counts_melt,
    aes(x = sample, colour = variable, fill = variable, y = value)
) +
    geom_bar(position = "stack", stat = "identity", width = 0.7) +
    coord_flip() +
    xlab("Sample") +
    ylab(ylabel) +
    scale_fill_manual(
        "",
        values = ggCols,
        guide = guide_legend(reverse = TRUE)
    ) +
    scale_colour_manual(values = ggCols, guide = FALSE) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6)) +
    theme_bw(base_size = cc1 * 1.3) +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = cc1),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
    )

nsamps <- ncol(gene_counts$counts)

ggsave(
    p1,
    file = paste0("library_composition.png"),
    device = "png",
    width = 8, height = (nsamps / 2.2),
    dpi = 300
)


#############################
## proportions plot
#############################
## get the proportions of reads per library

propCols <- counts_summary[c(
    "other",
    "rRNA"
    # all_biotypes
)] / counts_summary$mapped

propCols$sample <- counts_summary$sample

# rowSums(propCols) ## each row should sum to 1
prop_melt <- melt(
    propCols,
    id.vars = c("sample"),
    measure.vars = c(
        "other",
        "rRNA"
        # all_biotypes
    )
)
prop_melt$sample <- factor(
    prop_melt$sample,
    levels = rev(unique(sort(prop_melt$sample)))
)
prop_melt$variable <- factor(prop_melt$variable, levels = c(
    "rRNA",
    "other" # ,
    # non_rRNA_btypes
))

p2 <- ggplot(
    prop_melt,
    aes(x = sample, colour = variable, fill = variable, y = value)
) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    xlab("Sample") +
    ylab("Proportion of reads") +
    scale_fill_manual(
        "",
        values = ggCols,
        guide = guide_legend(reverse = TRUE)
    ) +
    scale_colour_manual(values = ggCols, guide = FALSE) +
    scale_y_continuous(labels = comma) +
    theme_bw(base_size = cc1 * 1.3) +
    theme(
        legend.position = "top",
        legend.text = element_text(size = cc1),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
    )

ggsave(
    p2,
    file = "library_composition_proportions.png",
    device = "png",
    width = 8, height = (nsamps / 2.2),
    dpi = 300
)
