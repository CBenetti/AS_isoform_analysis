###library
	library("GenomicRanges")
	library("AnnotationHub")
	library("ensembldb")
###loading
        ah <- AnnotationHub()
        ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 110))
        ahEdb <- ahDb[[1]]
	
writeLines(paste("sample","total_transcripts","total_annotated","total_unannotated","assigned_transcript","assigned_annotated","assigned_unannotated","unannotated_IR_containing_transcripts",sep="\t"),"out/stringtie/summary.txt")

###code
	x <- split(seq_along(transcripts(ahEdb)$gene_id),transcripts(ahEdb)$gene_id )
	x<-x[genes(ahEdb)$gene_id]
	nam <- transcripts(ahEdb)$tx_id_version
	p <- CharacterList(lapply(x,function(z){nam[z]}))
	ens_gene <-GRanges(seqnames=paste("chr",seqnames(genes(ahEdb)),sep=""),ranges=ranges(genes(ahEdb)),strand=strand(genes(ahEdb)), gene_id=genes(ahEdb)$gene_id, transcripts=p, gene_name=genes(ahEdb)$gene_name) 
	transcript_association_ens <- data.frame(tx_id=transcripts(ahEdb)$tx_id_version,gene_id=transcripts(ahEdb)$gene_id)
	rownames(transcript_association_ens) <- transcript_association_ens$tx_id
	dirl <- system("ls out/stringtie/out/", intern=T)
	dirl <- dirl[grep(".STAR",dirl)]
	sample_list_tot <- as.list(dirl)
	names(sample_list_tot) <- dirl
	sample_list_anno <- as.list(dirl)
	names(sample_list_anno) <- dirl
	sample_list_not_anno <- as.list(dirl)
	names(sample_list_not_anno) <- dirl
	for(d in dirl){
   
   ##loading file and subsetting annotated and not - annotated transcripts
		print(d)
		load(paste("out/stringtie/out/",d,"/stringtie_transcript.rda",sep=""))
		load(paste("out/stringtie/out/",d,"/stringtie_overlap_clean_cutoff_0.05.rda",sep=""))
		load(paste("out/stringtie/out/",d,"/IRF_overlap_clean_cutoff_0.05.rda",sep=""))
		anno <- transcript[which(transcript$ENS_ref_ID!="not_annotated")]
   
     #retrieving gene names for annotated transcripts
		anno$gene_name <- ens_gene[transcript_association_ens[anno$ENS_ref_ID,"gene_id"]]$gene_name
		not_anno <- transcript[which(transcript$ENS_ref_ID=="not_annotated")]

   ##associating references to univocally assigned stringtie genes
		coupling_anno <- unique(data.frame(gene_id=anno$gene_id, ref_id=transcript_association_ens[anno$ENS_ref_ID,"gene_id"]))
		reps <- split(seq_along(coupling_anno[,1]),coupling_anno[,1])
		le <- lapply(reps,length)
		coupling_anno <- coupling_anno[unlist(reps[which(le==1)]),] 
		rownames(coupling_anno) <- coupling_anno$gene_id

   ##using assigned references for unassigned transcript which match stringtie genes of assigned ones
		associated <- not_anno[which(not_anno$gene_id%in%rownames(coupling_anno)==TRUE)]
		associated$gene_name <- ens_gene[coupling_anno[associated$gene_id,"ref_id"]]$gene_name		
   
   ##mapping unassigned transcripts on ensembl
		not_anno_sub <- not_anno[which(not_anno$gene_id%in%coupling_anno$gene_id==FALSE)]
		o <- findOverlaps(ens_gene,not_anno_sub)
		pintersect(ens_gene[queryHits(o)],not_anno_sub[subjectHits(o)]) -> inters
		Perc<- width(inters)/width(not_anno_sub[subjectHits(o)])
		selection <- o[which(Perc>0.9)]
		perc_sel <- Perc[which(Perc>0.9)]
		sub <- cbind(ens_gene$gene_name[queryHits(selection)],not_anno_sub$transcript_id[subjectHits(selection)])
		reps <- split(seq_along(sub[,2]),sub[,2])
		reps <-lapply(reps, function(x) x[!(x %in% which(duplicated(sub)==TRUE))])
		le <- lapply(reps,length)
		tmp <- lapply(reps[which(le>1)],function(x){x[which(perc_sel[x]==max(perc_sel[x]))]})
		reps[names(tmp)] <- tmp
		reps <-lapply(reps, function(x) x[!(x %in% which(sub[,1]==""))])
		le <- lapply(reps,length)
		reps_found <- reps[which(le==1)]
		reps_found <- sort(unlist(reps_found))
		names(not_anno_sub) <- not_anno_sub$transcript_id
   ##selecting all assigned transcripts
		assigned <- not_anno_sub[names(reps_found)]
		assigned$gene_name <- ens_gene[queryHits(selection[reps_found])]$gene_name
		assigned <- c(assigned,associated)
			names(stringtie_overlap_clean$transcript_id) <-IRF_overlap_clean$gene
			li <- unlist(stringtie_overlap_clean$transcript_id)
			li<- setNames(names(li), li)
		if(length(which(unlist(stringtie_overlap_clean$transcript_id)%in%assigned$transcript_id==FALSE))>0){
			li_na <- li[which(names(li)%in%assigned$transcript_id==FALSE)]
			IR_iso <- not_anno[which(not_anno$transcript_id%in%names(li_na)==TRUE)]
			IR_iso$gene_name <- li_na[IR_iso$transcript_id]
			assigned <- c(assigned,IR_iso)			
		}
		tot <- c(anno,assigned)
		tot <- tot[-which(tot$gene_name=="")]
		untranslated_regions <- tot[which(tot$gene_name=="")]
		write(paste(d,length(transcript),length(anno),length(not_anno),length(tot),length(anno[which(anno$gene_name!="")]),length(assigned[which(assigned$gene_name!="")]),length(which(not_anno$transcript_id%in%names(li)==TRUE)),sep="\t"),file = "out/stringtie/summary.txt" , append = TRUE)
		sample_list_tot[[d]] <- v_tot <- table(tot$gene_name)
		sample_list_anno[[d]] <- a_tot <- table(anno[which(anno$gene_name!="")]$gene_name)
		sample_list_not_anno[[d]] <- u_tot <- table(assigned[which(assigned$gene_name!="")]$gene_name)
		isoforms <- data.frame(tot=rep(0,length(v_tot)),annotated=rep(0,length(v_tot)),unannotated=rep(0,length(v_tot)))
		rownames(isoforms) <- names(v_tot)		
		isoforms[names(a_tot),"annotated"] <- a_tot 
		isoforms$tot <- v_tot 
		isoforms[names(u_tot),"unannotated"] <- u_tot
		write.table(isoforms,paste("out/stringtie/out/",d,"/n_of_variants.txt",sep=""),sep="\t",row.names=TRUE)
		save(tot,file=paste("out/stringtie/out/",d,"/total_assigned_transcripts.rda",sep=""))
	}
###out
	
	total_isoforms <- matrix(0,ncol=length(dirl),nrow=length(unique(unlist(lapply(sample_list_tot,names)))))
	annotated_isoforms <- matrix(0,ncol=length(dirl),nrow=length(unique(unlist(lapply(sample_list_tot,names)))))
	not_annotated_isoforms <- matrix(0,ncol=length(dirl),nrow=length(unique(unlist(lapply(sample_list_tot,names)))))
	rownames(total_isoforms) <- rownames(annotated_isoforms) <- rownames(not_annotated_isoforms) <- unique(unlist(lapply(sample_list_tot,names)))
	colnames(total_isoforms) <- colnames(annotated_isoforms) <- colnames(not_annotated_isoforms) <- dirl
	for(d in dirl){
		total_isoforms[names(sample_list_tot[[d]]),d] <- sample_list_tot[[d]]		
		annotated_isoforms[names(sample_list_anno[[d]]),d] <- sample_list_anno[[d]]		
		not_annotated_isoforms[names(sample_list_not_anno[[d]]),d] <- sample_list_not_anno[[d]]		
	}
	save(total_isoforms,file="out/stringtie/total_isoforms.rda")
	save(annotated_isoforms,file="out/stringtie/annotated_isoforms.rda")
	save(not_annotated_isoforms,file="out/stringtie/not_annotated_isoforms.rda")

