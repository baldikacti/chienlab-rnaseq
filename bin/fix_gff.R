library(microseq)

ref_annot <- readGFF("references/genomic.gff")
ref_annot$Seqid <- "gi|221232939|ref|NC_011916.1|"
writeGFF(ref_annot, "references/ccna.gff")