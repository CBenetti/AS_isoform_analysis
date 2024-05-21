###library
	library(AnnotationHub)
	library(ensembldb)

###loading
	ah <- AnnotationHub()
        ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 110))
        ahEdb <- ahDb[[1]]
	dirl <- system("ls out/stringtie/out", intern=TRUE)
	dirl <- dirl[grep(".STAR",dirl)]
	transcripts_info <- as.list(dirl)
	names(transcripts_info) <- dirl
	transcripts_biotype_tpm <- transcripts_info
	for (d in dirl) {
                load(paste("out/stringtie/out/",d,"/stringtie_transcript.rda",sep=""))
                anno <- transcript[which(transcript$ENS_ref_ID!="not_annotated")]
                index <- match(anno$ENS_ref_ID,transcripts(ahEdb)$tx_id_version)
                bio <- transcripts(ahEdb)[index]$tx_biotype
		transcripts_info[[d]] <- table(bio)
		indList <- split(seq_along(bio), bio)
		transcripts_biotype_tpm[[d]] <- lapply(indList,function(x){mean(as.numeric(anno[x]$TPM))})
			
	}
	save(transcripts_info, file="out/stringtie/transcripts_info.rda")
	save(transcripts_biotype_tpm, file="out/stringtie/transcripts_biotype_tpm.rda")
	
