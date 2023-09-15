library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)
#library(BSgenome.Mmusculus.UCSC.mm10)

dir_path <- "results/QUO1006522/bwa_aln"

# create a function to read in data
# extend reads, get coverage, normalize to RPM
# and export to bw
bam2bw <- function(bamfile){
  
  extLen <- 150
  cat("opening:", bamfile, sep="\n")
  bd <- readGAlignments(bamfile)
  
  cat("convert to GRanges\n")
  mygr <- as(bd,"GRanges")
  cat("extending reads\n")
  mygr <- resize(mygr, extLen)
  mygr <- trim(mygr)
  
  totalReads <- length(mygr)
  
  cat("getting coverage\n")
  # get coverage                                                                                                      
  cov <- coverage(mygr)
  
  rpm <- lapply(cov, function(x) signif(10^6 * x/totalReads,3))
  rpm <- as(rpm,"SimpleRleList")
  
  # export rpm to bigWig                                                                                              
  outfile <- sub(".bam", "", bamfile)
  outfile <- paste(outfile, ".bw", sep="")
  cat(paste("exporting to bigwig", outfile, "\n", sep="\t"))
  export.bw(rpm, outfile)
}

# read the bam file names into a vector
bamfiles <- dir(path = dir_path, pattern = ".bam$", full.names = TRUE)

# for each element of our vector, call the bam2bw function
sapply(bamfiles, bam2bw)