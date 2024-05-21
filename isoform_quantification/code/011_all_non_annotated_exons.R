###library	
	library(ensembldb)
	library(GenomicRanges)
	library(AnnotationHub)
###loading
        ah <- AnnotationHub()
        ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 110))
        ahEdb <- ahDb[[1]]
        dirl <- system("ls out/stringtie/out", intern=TRUE)
        dirl <- dirl[grep(".STAR",dirl)]
	chrom <- c(as.character(seq(1,22,1)),"X","Y")
###code
	print(dirl[1])
	load(paste("out/stringtie/out/",dirl[1],"/stringtie.rda",sep=""))	
	seq <- sub("chr","",seqnames(stringtie))
	seqlevels(stringtie) <- levels(seq)
	seqnames(stringtie) <- seq
	o <- findOverlaps(stringtie,exons(ahEdb),type="equal")	
	novel_exons <- stringtie[-unique(queryHits(o))]
	load(paste("out/stringtie/out/",dirl[1],"/total_assigned_transcripts.rda",sep=""))
	annotation <- tot$gene_name
	names(annotation) <- tot$gene_id
	gene_id <- lapply(novel_exons$gene_id,function(x){if(length(unique(unique(unlist(x))%in%names(annotation)))<2){unique(annotation[x])}else{return("")}})
	gene_id <- unlist(lapply(gene_id,function(x) if(length(x)>1){""}else{x}))	
        chr <- as.character(seqnames(novel_exons))
        chr[which(chr%in%chrom==TRUE)] <- paste("chr",chr[which(chr%in%chrom==TRUE)],sep="")
	exon_list <- cbind(chr,rep("STRINGTIE",length(novel_exons)),rep("exon",length(novel_exons)),as.numeric(start(novel_exons)),as.numeric(end(novel_exons)),rep(".",length(novel_exons)),as.character(strand(novel_exons)),paste("gene_ID=",gene_id,";",sep=""))
	write.table(exon_list,"out/stringtie/stringtie_novel_exons.gtf",sep="\t",col.names=FALSE,row.names=F)	

	for(d in dirl[-1]){
		print(d)
		load(paste("out/stringtie/out/",d,"/stringtie.rda",sep=""))
		seq <- sub("chr","",seqnames(stringtie))
        	seqlevels(stringtie) <- levels(seq)
        	seqnames(stringtie) <- seq
        	o <- findOverlaps(stringtie,exons(ahEdb),type="equal")
        	novel_exons <- stringtie[-unique(queryHits(o))]
		#o <- findOverlaps(stringtie,novel_exons,type="equal")
		#new <- stringtie[-c(unique(queryHits(o)))]
		load(paste("out/stringtie/out/",d,"/total_assigned_transcripts.rda",sep=""))
		annotations <- tot$gene_name
		names(annotations) <- tot$gene_id
		gene_id <- lapply(novel_exons$gene_id,function(x){if(length(unique(unique(unlist(x))%in%names(annotation)))<2){unique(annotation[x])}else{return("")}})
		gene_id <- unlist(lapply(gene_id,function(x) if(length(x)>1){""}else{x}))	
		#novel_exons <- c(novel_exons,new)
		#gene_id <- c(gene_id,gene_id_s)
		chr <- as.character(seqnames(novel_exons))
		chr[which(chr%in%chrom==TRUE)] <- paste("chr",chr[which(chr%in%chrom==TRUE)],sep="")
		exon_list <- cbind(chr,rep("STRINGTIE",length(novel_exons)),rep("exon",length(novel_exons)),as.numeric(start(novel_exons)),as.numeric(end(novel_exons)),rep(".",length(novel_exons)),as.character(strand(novel_exons)),paste("gene_ID=",gene_id,";",sep=""))
		write.table(exon_list,"out/stringtie/stringtie_novel_exons.gtf",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)	
		system("awk '!visited[$0]++' out/stringtie/stringtie_novel_exons.gtf  > out/stringtie/tmp.gtf",intern=T)
		system("mv out/stringtie/tmp.gtf out/stringtie/stringtie_novel_exons.gtf",intern=T)
	}

###out
