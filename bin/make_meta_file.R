#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("--sample_file"),
              type = "character", default = NULL,
              help = "sample metadata tsv file containing following columns.\n Headers: sample, file1, file2, group, rep_no"
  ),
  make_option(c("--data_dir"),
              type = "character", default = NULL,
              help = "Directory containing RNA-seq data"
  ),
  make_option(c("--outf"),
              type = "character",
              help = "File for results"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read data and add data_dir to file paths in file1 and file2 columns
sample_dat <- suppressWarnings(read.delim(opt$sample_file, header = TRUE, sep = "\t", na.strings = c("", "NA")))
sample_dat <- sample_dat[!is.na(sample_dat$sample),]
sample_dat$file1 <- paste0(opt$data_dir, "/", sample_dat$file1)
sample_dat$file2 <- ifelse(sample_dat$paired == 1, paste0(opt$data_dir, "/", sample_dat$file2), "")
# Export sample_file to outf location
write.table(sample_dat, file = opt$outf, sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)