###library
	library(Biobase)
	library(circlize)
	library(ComplexHeatmap)
	library(xlsx)
	library(ggplot2)
###loading
	load("out/annotated_isoforms_cell_bank.rda")
	gene_counts <- read.table("data/fastq.featurecounts.ribo.ex.count",sep="\t",header=T)	
	s_r <- read.xlsx("data/Sample_replicates.xlsx",sheetIndex=1,header=F)
###code
	mat <- gene_counts[,-c(1:6)]
	colnames(mat) <- paste(colnames(mat),".STAR",sep="")
	annotated_genes_counts_cell_bank <- ExpressionSet(assayData=as.matrix(mat[,colnames(exprs(isoforms_cell_bank$TPM))]))
	fData(annotated_genes_counts_cell_bank) <- gene_counts[,c(1:6)]
	pData(annotated_genes_counts_cell_bank) <- pData(isoforms_cell_bank$TPM)
	cpm <- apply(mat,2,function(x){x*1000000/sum(x)})
	genes_expressed_samples <- apply(cpm,2,function(x){length(which(x>1))})
	genes_expressed_samples_oc <- genes_expressed_samples[colnames(annotated_genes_counts_cell_bank)]
	pdf("out/expressed_genes.pdf")
plot(ggplot(data=data.frame(cells=factor(names(genes_expressed_samples),levels=names(genes_expressed_samples)[order(genes_expressed_samples,decreasing=F)]),genes=genes_expressed_samples), 
aes(x=cells, y=genes)) +
  geom_bar(stat="identity",fill="#00AFBB")+theme(axis.text.x=element_blank())+ggtitle("Genes with al least 1 CPM"))
plot(ggplot(data=data.frame(cells=factor(names(genes_expressed_samples_oc),levels=names(genes_expressed_samples_oc)[order(genes_expressed_samples_oc,decreasing=F)]),genes=genes_expressed_samples_oc), 
aes(x=cells, y=genes)) +
  geom_bar(stat="identity",fill="#00AFBB")+theme(axis.text.x=element_blank())+ggtitle("Genes with al least 1 CPM"))
	 dev.off()
	cpm_log <- apply(cpm,1,function(x){log2(x+0.1)})
	ltr <- t(apply(cpm_log,2,function(x){x-mean(x)}))
	rownames(ltr) <-fData(annotated_genes_counts_cell_bank)$Geneid
	v <- apply(cpm,1,var)
	ltr_ord <- ltr[order(v,decreasing=T),]
	genes_to_filter_out <- apply(cpm,1,function(x){if(quantile(x,0.8)>10){0}else{1}})
	genes_to_filter_out_names <- fData(annotated_genes_counts_cell_bank)$Geneid[which(genes_to_filter_out==1)]
	ltr_ord_filt <- ltr_ord[setdiff(rownames(ltr_ord),genes_to_filter_out_names),]
	ltr_ord_filt <- ltr_ord_filt[,setdiff(colnames(ltr_ord_filt),c("COGA122.STAR","COLO94H2.STAR","COLO3202.STAR"))]
	col2 <- rep("black",length(colnames(ltr_ord_filt)))
	names(col2) <- colnames(ltr_ord_filt)
	col2[s_r[,1]] <- "orange"
	col2[s_r[,2]] <- "darkred"
	pdf("out/cell_similarity_genes.pdf",width=15)
	draw(Heatmap(ltr_ord_filt,cluster_rows=F,show_row_names=F,col=colorRamp2(c(-2,0,2),c("blue","white","red")),column_names_gp=gpar(col=col2,fontsize=5),column_names_side = "top"))
	draw(Heatmap(ltr_ord_filt[,colnames(exprs(isoforms_cell_bank$TPM))],cluster_rows=F,show_row_names=F,col=colorRamp2(c(-2,0,2),c("blue","white","red")),column_names_gp=gpar(col=col2[colnames(exprs(isoforms_cell_bank$TPM))],fontsize=5),column_names_side = "top"))	
	dev.off()
###out

	save(annotated_genes_counts_cell_bank,file="out/annotated_genes_counts_cell_bank.rda")
