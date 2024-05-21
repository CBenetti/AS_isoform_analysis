###library and loading
	load("out/IRFinder/IR_ratio.rda")
	load("data/potential_exons_ensembl_110.rda")
	if(!require("GenomicRanges", quietly=TRUE)) install.packages("GenomicRanges")
	library(GenomicRanges)

###code
##building IRFinder GRobject
	IRFinder_list <- paste("chr",unlist(lapply(rownames(IR_ratio),function(x){paste(c(strsplit(x,":")[[1]][c(1,2,3,7)]),collapse=":")})),sep="")
	IRFinder_type <- unlist(lapply(rownames(IR_ratio),function(x){strsplit(x,":")[[1]][6]}))
	IRF <- GRanges(seqnames = Rle(unlist(lapply(IRFinder_list,function(x){strsplit(x,":")[[1]][1]}))),ranges =
IRanges(start=as.numeric(unlist(lapply(IRFinder_list,function(x){strsplit(x,":")[[1]][2]}))), end
=as.numeric(unlist(lapply(IRFinder_list,function(x){strsplit(x,":")[[1]][3]})))),  strand =
Rle(gsub("[.]","*",unlist(lapply(IRFinder_list,function(x){strsplit(x,":")[[1]][4]})))),type=unlist(lapply(rownames(IR_ratio),function(x){strsplit(x,":")[[1]][6]})),gene=unlist(lapply(rownames(IR_ratio),
function(x){strsplit(x,":")[[1]][4]})))
	save(IRF,file="out/IRFinder/IRFinder_GR_object.rda")
##for each sample:
dirl <- system("ls out/stringtie/out/", intern=T)
dirl <- dirl[which(dirl%in%colnames(IR_ratio==TRUE))]
hit_list <- as.list(dirl)
names(hit_list) <- dirl
cutoffs <- c(0.05,0.10,0.2,0.3)
names(cutoffs) <- paste("above",cutoffs,sep="")
hit_list <- lapply(hit_list,function(x){as.list(cutoffs)})
for(d in dirl){
	stringtie_tot <- readLines(paste("out/stringtie/out/",d,"/output.gtf",sep=""))[-c(1:2)]

#building sample stringtie object
	#stringtie_tot <- lapply(stringtie_tot,function(x){strsplit(x,"\t")})
	stringtie_exons <- stringtie_tot[grep("exon",stringtie_tot)]
	stringtie_exons <- strsplit(stringtie_exons,"\t")
	stringtie_exons_name <- paste(sapply(stringtie_exons,"[[",1),sapply(stringtie_exons,"[[",4),sapply(stringtie_exons,"[[",5),sapply(stringtie_exons,"[[",7),sep=":")
	stringtie_exons_unique <- stringtie_exons[-which(duplicated(stringtie_exons_name))]
	indList <- split(seq_along(stringtie_exons_name), stringtie_exons_name)
	stringtie_exon_name_unique <- paste(sapply(stringtie_exons_unique,"[[",1),sapply(stringtie_exons_unique,"[[",4),sapply(stringtie_exons_unique,"[[",5),sapply(stringtie_exons_unique,"[[",7),sep=":")
	indList <- indList[stringtie_exon_name_unique]
	stringtie_exons_info <- sapply(stringtie_exons,"[[",9)
	stringtie_exons_info <- gsub('["]',"",stringtie_exons_info)
	stringtie_exons_info <- gsub(" ",";",stringtie_exons_info)
	stringtie_exons_info <- gsub(";;",";",stringtie_exons_info)
	stringtie_exons_info <- strsplit(stringtie_exons_info,";")
	stringtie_exons_info <- lapply(stringtie_exons_info,function(x){a=x[seq(1,length(x),2)];b=x[seq(2,length(x),2)];names(b)=a;b})
	gene_id <- sapply(stringtie_exons_info,"[[","gene_id")
	exon_number <- sapply(stringtie_exons_info,"[[","exon_number")
	transcript_id <- sapply(stringtie_exons_info,"[[","transcript_id")
	cov <- sapply(stringtie_exons_info,"[[","cov")
	annotation <- lapply(stringtie_exons_info,function(x){if(length(x)>5){x["reference_id"]}else{"not-annotated"}})
	gene_id <- CharacterList(lapply(indList,function(x){gene_id[x]}))
	exon_number <- CharacterList(lapply(indList,function(x){exon_number[x]}))
	transcript_id <- CharacterList(lapply(indList,function(x){transcript_id[x]}))
	cov <- CharacterList(lapply(indList,function(x){cov[x]}))
	annotation_state=CharacterList(lapply(indList,function(x){unique(annotation[x])}))

	stringtie <- GRanges(seqnames = Rle(sapply(stringtie_exons_unique,"[[",1)),ranges =
IRanges(start=as.numeric(sapply(stringtie_exons_unique,"[[",4)), end =as.numeric(sapply(stringtie_exons_unique,"[[",5))),
strand = Rle(gsub("[.]","*",sapply(stringtie_exons_unique,"[[",7))),gene_id=gene_id,exon_number=exon_number,transcript_id=transcript_id,cov=cov,annotation_state=annotation_state)
 

	#stringtie_exons_nu <- unlist(lapply(stringtie_tot,function(x){if(unlist(x)[3]=="exon"){paste(unlist(x)[1],unlist(x)[4],unlist(x)[5],unlist(x)[7],sep=":")}}))
	#stringtie_info_nu <- unlist(lapply(stringtie_tot,function(x){if(unlist(x)[3]=="exon"){unlist(x)[9]}}))
	#stringtie_info_nu <- gsub("gene_id ","",stringtie_info_nu)
	#stringtie_info_nu <- gsub("transcript_id ","",stringtie_info_nu)
	#stringtie_info_nu <- gsub("exon_number ","",stringtie_info_nu)
	#exon_id <- unlist(lapply(strsplit(stringtie_info_nu,"[\"]"),function(x){x[2]}))
	#transcript_id <- unlist(lapply(strsplit(stringtie_info_nu,"[\"]"),function(x){x[4]}))
	#exon_number  <- unlist(lapply(strsplit(stringtie_info_nu,"[\"]"),function(x){x[6]}))
	#exon_id <- lapply(stringtie_exons,function(x){exon_id[which(stringtie_exons_nu ==x)]})
	#transcript_id <- lapply(stringtie_exons,function(x){transcript_id[which(stringtie_exons_nu ==x)]})
	#exon_number <- lapply(stringtie_exons,function(x){exon_number[which(stringtie_exons_nu ==x)]})
	#unique(stringtie_exons%in%pot_exons)
	#unique(stringtie_exons%in%IRFinder_list)
	#annotation_state <- rep("annotated", length(stringtie_exons))
	#annotation_state[which(stringtie_exons%in%pot_exons==FALSE)] <- "not-annotated"
	#not_annotated_exons <- stringtie_exons[which(stringtie_exons%in%pot_exons==FALSE)]
	 #stringtie <- GRanges(seqnames = Rle(unlist(lapply(stringtie_exons,function(x){strsplit(x,":")[[1]][1]}))),ranges = 
#IRanges(start=as.numeric(unlist(lapply(stringtie_exons,function(x){
#strsplit(x,":")[[1]][2]}))), end =as.numeric(unlist(lapply(stringtie_exons,function(x){strsplit(x,":")[[1]][3]})))),  
#strand = Rle(gsub("[.]","*",unlist(lapply(stringtie_exons,function(x){strsplit(x,":")[[1]][4]})))),transcript_id=CharacterList(transcript_id),exon_number=CharacterList(exon_number),annotation_state=annotation_state)

#building transcript GRobject
	stringtie_transcript <- stringtie_tot[grep("\ttranscript",stringtie_tot)]
	stringtie_transcript <- strsplit(stringtie_transcript,"\t")
	stringtie_transcript_info <- sapply(stringtie_transcript,"[[",9)
	stringtie_transcript_info <- gsub('["]',"",stringtie_transcript_info)
	stringtie_transcript_info <- gsub(" ",";",stringtie_transcript_info)
	stringtie_transcript_info <- gsub(";;",";",stringtie_transcript_info)
	stringtie_transcript_info <- strsplit(stringtie_transcript_info,";")
	stringtie_transcript_info <- lapply(stringtie_transcript_info,function(x){a=x[seq(1,length(x),2)];b=x[seq(2,length(x),2)];names(b)=a;b})
	transcript <- GRanges(seqnames = Rle(sapply(stringtie_transcript,"[[",1)),ranges = 
IRanges(start=as.numeric(sapply(stringtie_transcript,"[[",4)), end =as.numeric(sapply(stringtie_transcript,"[[",5))),
strand = Rle(gsub("[.]","*",sapply(stringtie_transcript,"[[",7))), gene_id=unlist(lapply(stringtie_transcript_info,function(x){x["gene_id"]})),
transcript_id=unlist(lapply(stringtie_transcript_info,function(x){x["transcript_id"]})),cov=unlist(lapply(stringtie_transcript_info,function(x){x["cov"]})),
TPM=unlist(lapply(stringtie_transcript_info,function(x){x["TPM"]})),FPKM=unlist(lapply(stringtie_transcript_info,function(x){x["FPKM"]})),
ENS_ref_ID=unlist(lapply(stringtie_transcript_info,function(x){if("reference_id"%in%names(x)){x["reference_id"]}else{"not_annotated"}})))

#distribution of length: quality control check
	not_anno <-  stringtie[grep("not-annotated",as.list(stringtie$annotation_state)),]
	na_length <- as.numeric(end(not_anno))-as.numeric(start(not_anno))
	stringtie_length <- as.numeric(end(stringtie))-as.numeric(start(stringtie))
	ns_length <- unlist(lapply(pot_exons,function(x){as.numeric(strsplit(x,":")[[1]][3])-as.numeric(strsplit(x,":")[[1]][2])}))
	pdf(paste(c("out/stringtie/out/",d,"/distribitions.pdf"),collapse=""))
	xl <- c(1000,5000,max(na_length))
	for(i in 1:3){
		plot(names(table(ns_length)),100*table(ns_length)/length(ns_length),xlim=c(0,xl[i]),pch=".",xlab="length",ylab="%")
		points(names(table(stringtie_length)),100*table(stringtie_length)/length(stringtie_length),xlim=c(0,xl[i]),pch=".",col="red")
		points(names(table(na_length)),100*table(na_length)/length(na_length),xlim=c(0,xl[i]),pch=".",col="green")
		legend("topright",legend=c("ensemble","stringtie","non-annotated in stringtie"),pch=20,col=c("black","red","green"))
	}
	dev.off()
	for(c in cutoffs){ 
##overlaps between introns form IRFinder and stringtie not annotated
		IRF_sample <- IRF[which(IR_ratio[,d]>c)]
		overlap <- findOverlaps(IRF_sample, not_anno)
		if(length(overlap)>0){
			pintersect(IRF_sample[queryHits(overlap),],not_anno[subjectHits(overlap),]) -> inters
			inters_value <- vector()
			for(i in 1:length(inters$hit)){inters_value[i] <- width(inters[i,])}
			overlap_good <- overlap[which(inters_value>1),]
##subset GRobject
			IRF_overlap <- IRF_sample[queryHits(overlap_good),]
			stringtie_overlap <- not_anno[subjectHits(overlap_good)]
			IRF_overlap_clean <- IRF_overlap[which(IRF_overlap$type=="clean")]
			stringtie_overlap_clean <- stringtie_overlap[which(IRF_overlap$type=="clean")]
			dup <- unique(IRF_overlap_clean[which(duplicated(IRF_overlap_clean))])
	
			counts <- unlist(lapply(stringtie_overlap_clean$transcript_id,function(x){if(length(x)==1){transcript$TPM[which(transcript$transcript_id==x)]}else{a=match(x,transcript$transcript_id,nomatch=0);
sum(as.numeric(transcript$TPM[a]))/length(a)}})) 
			if(length(dup)>0){for(i in 1:length(dup)){a=match(unlist(stringtie_overlap_clean$transcript_id[which(IRF_overlap_clean==dup[i])]),transcript$transcript_id);
 counts[which(IRF_overlap_clean==dup[i])]=sum(as.numeric(transcript$TPM[a]))/length(a)}}
			names_counts  <- paste(as.vector(seqnames(IRF_overlap_clean)),start(IRF_overlap_clean),end(IRF_overlap_clean),strand(IRF_overlap_clean),IRF_overlap_clean$gene,sep=":")
			tmp <- unique(cbind(names_counts,counts))
			counts_uniqued <- as.numeric(tmp[,2])
			names(counts_uniqued) <- tmp[,1]
			hit_list[[d]][[paste("above",c,sep="")]] <- counts_uniqued
			save(stringtie_overlap_clean,file=paste("out/stringtie/out/",d,"/stringtie_overlap_clean_cutoff_",c,".rda",sep=""))
			save(IRF_overlap_clean,file=paste("out/stringtie/out/",d,"/IRF_overlap_clean_cutoff_",c,".rda",sep=""))
			}
		}
	
#save
	save(transcript,file=paste("out/stringtie/out/",d,"/stringtie_transcript.rda",sep=""))
	save(stringtie,file=paste("out/stringtie/out/",d,"/stringtie.rda",sep=""))
	save(not_anno,file=paste("out/stringtie/out/",d,"/not_annotated.rda",sep=""))
	print(d)
	}
	save(hit_list, file="out/stringtie/out/list_count_results.rda")
	for(c in cutoffs){
		nam <- unique(unlist(lapply(hit_list,function(x){names(x[[paste("above",c,sep="")]])})))
		raw_tab <- matrix(0,ncol=length(hit_list),nrow=length(nam))
		rownames(raw_tab) <- nam
		colnames(raw_tab) <- names(hit_list)
		for(i in names(hit_list)){raw_tab[names(hit_list[[i]][[paste("above",c,sep="")]]),i] <- hit_list[[i]][[paste("above",c,sep="")]]}
		save(raw_tab, file=paste("out/stringtie/out/count_raw_results_cutoff_",c,".rda",sep=""))
	}
