###library
	library(GenomicRanges)
	library(ensembldb)
	library(AnnotationHub)
	library(Biobase)
###loading
        ah <- AnnotationHub()
        ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 110))
        ahEdb <- ahDb[[1]]

###code
        dirl <- system("ls out/stringtie/out/", intern=T)
        dirl <- dirl[grep(".STAR",dirl)]
	cov <- matrix(0,length(transcripts(ahEdb)),length(dirl))
	colnames(cov) <- dirl
	rownames(cov) <- transcripts(ahEdb)$tx_id_version
	FPKM <- TPM <- cov
	gene_name <- genes(ahEdb)[transcripts(ahEdb)$gene_id]$gene_name
	biotype <- transcripts(ahEdb)$tx_biotype
	t_length <- width(transcripts(ahEdb))
	for(d in dirl){
		print(d)
		load(paste("out/stringtie/out/",d,"/stringtie_transcript.rda",sep=""))
		anno <- transcript[which(transcript$ENS_ref_ID!="not_annotated")]
		cov[anno$ENS_ref_ID,d] <- as.numeric(anno$cov)
		TPM[anno$ENS_ref_ID,d] <- as.numeric(anno$TPM)
		FPKM[anno$ENS_ref_ID,d] <- as.numeric(anno$FPKM)
	}
	non_null <- apply(cov,1,sum)
	coverage <- ExpressionSet(assayData=cov[which(non_null>0),])
	TPM <- ExpressionSet(assayData=TPM[which(non_null>0),])
	FPKM <- ExpressionSet(assayData=FPKM[which(non_null>0),])
	fData(coverage) <- fData(TPM) <- fData(FPKM) <- data.frame(Geneid=gene_name[which(non_null>0)], biotype=biotype[which(non_null>0)], length=t_length[which(non_null>0)])
	isoforms_cell_bank <- list(cov=coverage,TPM=TPM,FPKM=FPKM)
	save(isoforms_cell_bank, file="out/stringtie/isoforms_cell_bank.rda")
