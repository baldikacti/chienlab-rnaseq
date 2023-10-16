library(tidyverse)

ref_annot <- microseq::readGFF("references/genomic.gff")
ref_annot$Seqid <- "gi|221232939|ref|NC_011916.1|"
microseq::writeGFF(ref_annot, "references/ccna.gff")

ref_annot <- microseq::readGFF("references/genomic.gtf")
ref_annot$Seqid <- "gi|221232939|ref|NC_011916.1|"
microseq::writeGFF(ref_annot, "references/ccna.gtf")