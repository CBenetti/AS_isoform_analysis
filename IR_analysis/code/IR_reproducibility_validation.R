###library
	library(Biobase)
###loading
	load("data/count_raw_results_cutoff_0.3.rda")
	read.table("data/fastq.featurecounts.ribo.ex.count",sep="\t",header=T)-> exp
	
###code
annotated_genes_counts_cell_bank <- ExpressionSet(assayData=as.matrix(exp[,-c(1:6)]))
colnames(annotated_genes_counts_cell_bank)<- paste(colnames(annotated_genes_counts_cell_bank),"STAR",sep=".")
fData(annotated_genes_counts_cell_bank) <- exp[,c(1:6)]	
pdf("out/IR_reproducibility_validation.pdf")	
	samples <- c("C10.STAR","DLD1.STAR","LS180.STAR","SNU1040.STAR","LIM1215.STAR","HROC277.STAR","HROC24.STAR","SW620.STAR")
	replicates <- c("C10S1.STAR","DLD1S12.STAR","LS180S23.STAR","SNU1040S11.STAR","LIM1215S22.STAR","HROC277MET2.STAR","HROC24S4.STAR","SW620S14.STAR")
	names(replicates) <- samples
	for (i in samples){
		hits_sample <- which(raw_tab[,i]>0)
		hits_replicates <- which(raw_tab[,replicates[i]]>0)
		intersection <- intersect(hits_sample,hits_replicates)
		hits_sample_only <- setdiff(hits_sample,intersection)
		hits_replicates_only <- setdiff(hits_replicates,intersection)
		s <- strsplit(rownames(raw_tab)[hits_sample_only],":")
		r <- strsplit(rownames(raw_tab)[hits_replicates_only],":")
		b <- strsplit(rownames(raw_tab)[intersection],":")
		s <- sapply(s[which(lapply(s,length)>4)],"[[",5)
		r <- sapply(r[which(lapply(r,length)>4)],"[[",5)
		b <- sapply(b[which(lapply(b,length)>4)],"[[",5)
		ord_sample <- order(exprs(annotated_genes_counts_cell_bank)[,i])
		ord_rep <- order(exprs(annotated_genes_counts_cell_bank)[,i])
		col <- ifelse(fData(annotated_genes_counts_cell_bank)$Geneid%in%b==TRUE,"shared IR","no IR")
		col[which(fData(annotated_genes_counts_cell_bank)$Geneid%in%s==TRUE)] <- "IRsample"
		col[which(fData(annotated_genes_counts_cell_bank)$Geneid%in%r==TRUE)] <-"IRreplicate"
		par(mfrow=c(1,2))
		
boxplot(sort(log(exprs(annotated_genes_counts_cell_bank)[,replicates[i]]))~col[ord_rep],main=replicates[i],ylab="log(counts)",xlab="",las=2, 
col=c("blue","red","black","purple"))
		
boxplot(sort(log(exprs(annotated_genes_counts_cell_bank)[,i]))~col[ord_sample],main=i,ylab="log(counts)",xlab="",las=2,col=c("blue","red","black","purple"))
	}
	dev.off()
