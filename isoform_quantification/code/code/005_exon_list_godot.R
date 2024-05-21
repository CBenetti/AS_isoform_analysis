###library
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("AnnotationHub", quietly = TRUE)) BiocManager::install(AnnotationHub)
if (!require("ensembldb", quietly = TRUE)) BiocManager::install(ensembldb)
library(AnnotationHub)
library(ensembldb)

###code
	ah <- AnnotationHub()
	ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 110))
	ahEdb <- ahDb[[1]]
	exon <- exonicParts(ahEdb)
	str <- as.vector(strand(exon))
        chr <- paste("chr",as.vector(seqnames(exon)),sep="")
        gene <- mcols(exon)$gene_id  
	start <- start(ranges(exon))
	end <- end(ranges(exon))
	pot_exons <- paste(chr,start,end,str,sep=":")

###save
	save(pot_exons,file="data/potential_exons_ensembl_110.rda")
