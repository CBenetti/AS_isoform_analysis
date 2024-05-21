###library
	library(ensembldb)
	library(AnnotationHub)
	library(ggplot2)
###loading
	load("out/annotated_isoforms_cell_bank.rda")
	load("out/annotated_genes_counts_cell_bank.rda")
###code
	ah <- AnnotationHub()
        ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 110))
        ahEdb <- ahDb[[1]]
	transcr_data <- transcripts(ahEdb)
	undefined <- rownames(exprs(isoforms_cell_bank$TPM))[which(fData(isoforms_cell_bank$TPM)$Geneid=="")]
	fData(isoforms_cell_bank$cov)$Geneid[which(fData(isoforms_cell_bank$cov)$Geneid=="")] <- transcr_data[match(undefined,transcr_data$tx_id_version,nomatch=0)]$gene_id
	transcriptLengths(txdb=ahEdb,with.cds_len=TRUE) -> tl
	fData(isoforms_cell_bank$cov)$transcript_length <- tl$tx_len[which(match(transcr_data$tx_id_version,rownames(exprs(isoforms_cell_bank$cov)),nomatch=0)>0)]
	count_matrix_75 <- 
apply(exprs(isoforms_cell_bank$cov),2,function(x){ceiling(x*fData(isoforms_cell_bank$cov)$transcript_length/75)})
	count_matrix_150 <- 
apply(exprs(isoforms_cell_bank$cov),2,function(x){ceiling(x*fData(isoforms_cell_bank$cov)$transcript_length/150)})
	isoforms_cell_bank$counts_75 <- ExpressionSet(assayData=as.matrix(count_matrix_75))
	isoforms_cell_bank$counts_150 <- ExpressionSet(assayData=as.matrix(count_matrix_150))
	fData(isoforms_cell_bank$counts_75) <- fData(isoforms_cell_bank$counts_150) <- fData(isoforms_cell_bank$cov)
	pData(isoforms_cell_bank$counts_75) <- pData(isoforms_cell_bank$counts_150) <- pData(isoforms_cell_bank$cov)
	sum_150 <- apply(exprs(isoforms_cell_bank$counts_150),2,sum)
	sum_75 <- apply(exprs(isoforms_cell_bank$counts_75),2,sum)
	s_genes <- apply(exprs(annotated_genes_counts_cell_bank),2,sum)
	cpm <- apply(exprs(isoforms_cell_bank$counts_75),2,function(x){x*1000000/sum(x)})
	transcr <- apply(cpm,2,function(x){length(which(x>1))})
	genes <- apply(cpm,2,function(x){length(unique(fData(isoforms_cell_bank$counts_75)$Geneid[which(x>1)]))})
	gps <- lapply(unique(fData(isoforms_cell_bank$counts_75)$Geneid),function(x){
		tab=cpm[which(fData(isoforms_cell_bank$counts_75)$Geneid==x),]
		if(is.null(dim(tab))!=TRUE){
			return(apply(tab,2,function(y){length(which(y>0))}))
		}else{
			return(rep(1,length(tab)))
		}
	})
	pdf("out/represented_isoforms.pdf")
	plot(ggplot(data.frame(values=unlist(gps)), aes(x=values)) +
    	geom_histogram(aes(y=..density..),fill="#00AFBB", alpha=0.8,binwidth=1) +
    	ggtitle("Number of represented isoforms per gene")+ 
  	stat_function(fun = dpois, args = list(lambda=mean(as.numeric(unlist(gps)))),geom="area",fill="#00AFBB",alpha=0.25)+ stat_function(fun = dpois, args = 
	list(lambda=mean(as.numeric(unlist(gps))))))	
	dev.off()	
###out
	pdf("out/counts_conversion_results.pdf")
	plot(s_genes,sum_75,pch=20,xlab="lib size counts\n for genes",ylab="lib size 
converted counts\nfor transcripts",main=paste("Converted counts with read 
length = 75\n cor=",cor(s_genes,sum_75),sep=" "))
	plot(s_genes,sum_150,pch=20,xlab="lib size genes",ylab="lib size 
converted counts\nfor transcripts",main=paste("Converted counts with read 
length = 150\n cor=",cor(s_genes,sum_150),sep=" "))
	cols <- colnames(isoforms_cell_bank$counts_75)[order(transcr,decreasing=F)]
	 db <- data.frame(cells=factor(rep(cols,2),levels=cols),type=factor(sort(rep(c("transcript","genes"),length(cols))),levels=c("transcript","genes")),number=c(genes[cols],transcr[cols]-genes[cols]))
	plot(ggplot(data=db, aes(x=cells, y=number, fill=type)) + geom_bar(stat="identity")+ggtitle("Transcripts and genes with more than 1 
CPM")+theme(axis.text.x=element_blank()))
	dev.off()
	save(isoforms_cell_bank,file="out/annotated_isoforms_cell_bank_counts.rda")
