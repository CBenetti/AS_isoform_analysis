###library
	library(ggplot2)
	library(Rtsne)
	library(RColorBrewer)
	library(xlsx)
	library(Biobase)

###loading
	load("out/DGE_results.rda")
	load("out/annotated_genes_counts_cell_bank.rda")
	load("out/annotated_isoforms_cell_bank_proportions.rda")
	load("out/DIE_results_DRIMSeq_withL2R.rda")
###code
	cells <- colnames(exprs(annotated_genes_counts_cell_bank))[grep("only",pData(annotated_genes_counts_cell_bank)$CRIS_class)]
	group <- pData(annotated_genes_counts_cell_bank)$CRIS_class[grep("only",pData(annotated_genes_counts_cell_bank)$CRIS_class)]
	names(group) <- cells
	group_colors <- c("red","darkgreen","lightblue","orange","blue")
	names(group_colors) <- unique(group)
##finding variable genes
	means <- rowMeans(exprs(annotated_genes_counts_cell_bank))
	vars <- apply(exprs(annotated_genes_counts_cell_bank),1,var)
	cv2 <- vars/means^2
	var_genes <-fData(annotated_genes_counts_cell_bank)$Geneid[order(cv2,decreasing=T)]
##selecting significative genes (both positive and negative)
	DGE_genes <- lapply(DGE[["CRIS_adjusted"]],function(x){rownames(x$table)[which(x$table$FDR < 0.01 & abs(x$table$logFC) > 2)]})
	DTU_isoforms <- lapply(DIE[["CRIS_adjusted"]],function(x){x$rownam[which(x$FDR < 0.01 & x$proportionSD > 0.1)]})
	 
##selecting significative genes (only positive)
	DGE_genes_positive <- lapply(DGE[["CRIS_adjusted"]],function(x){rownames(x$table)[which(x$table$FDR < 0.01 & x$table$logFC > 2)]})
	dup <- unlist(DGE_genes_positive)[which(duplicated(unlist(DGE_genes_positive))==TRUE)]
        lapply(DGE_genes_positive,function(x){setdiff(x,dup)})-> DGE_genes_positive
	DTU_isoforms_positive <- lapply(DIE[["CRIS_adjusted"]],function(x){x$rownam[which(x$FDR < 0.01 & x$logFC >0 & x$proportionSD > 0.1)]})
#biotype correlation
	rownames(fData(isoforms_cell_bank$isoform_proportions)) <- rownames(exprs(isoforms_cell_bank$isoform_proportions))
	pdf("out/biotype_positively_correlated_isoforms.pdf")
        biot <- lapply(DTU_isoforms_positive,function(x){table(fData(isoforms_cell_bank$isoform_proportions)[x,"biotype"])})
        types <-  lapply(biot,names)                
        col <- brewer.pal(n = length(unique(unlist(types))), name = "Paired")
        names(col) <- unique(unlist(types))
        lapply(biot,function(x){ggplot(data=data.frame(events=as.numeric(x),biotype=names(x)),aes(x=biotype,y=events,fill=biotype))+geom_bar(stat="identity")+scale_fill_manual(values=col[names(x)])+                                                           
        theme_classic()+theme(axis.text.x=element_blank())})
        dev.off()

##filtering for expression
	m_h_expressed_genes <- fData(annotated_genes_counts_cell_bank)$Geneid[which(means > 100)]
	DGE_genes_positive_filt <- lapply(DGE_genes_positive,function(x){x[which(x%in%m_h_expressed_genes==TRUE)]})
	DTU_isoforms_positive_filt <- lapply(DTU_isoforms_positive,function(x){x[which(fData(isoforms_cell_bank$isoform_proportions)[x,"Geneid"]%in%m_h_expressed_genes==TRUE)]})
	DGE_genes_filt <- lapply(DGE_genes,function(x){x[which(x%in%m_h_expressed_genes==TRUE)]})
	DTU_isoforms_filt <- lapply(DTU_isoforms,function(x){x[which(fData(isoforms_cell_bank$isoform_proportions)[x,"Geneid"]%in%m_h_expressed_genes==TRUE)]})
	var_genes_filt <- var_genes[which(var_genes%in%m_h_expressed_genes==TRUE)]
	var_genes_filt <- var_genes_filt[c(1:100)]

##creating lists	
	x_list <- list(
		var_genes = exprs(annotated_genes_counts_cell_bank)[which(var_genes%in%fData(annotated_genes_counts_cell_bank)$Geneid==TRUE),cells],
		DGE_all = exprs(annotated_genes_counts_cell_bank)[which(unique(unlist(DGE_genes))%in%fData(annotated_genes_counts_cell_bank)$Geneid==TRUE),cells],
		DGE_pos = exprs(annotated_genes_counts_cell_bank)[which(unlist(DGE_genes_positive)%in%fData(annotated_genes_counts_cell_bank)$Geneid==TRUE),cells],
		DTU_all = exprs(isoforms_cell_bank$isoform_proportions)[unlist(unique(DTU_isoforms)),cells],
		DTU_all_filt =exprs(isoforms_cell_bank$isoform_proportions)[unlist(unique(DTU_isoforms)),setdiff(cells,c("COLO94H.STAR","SW1463.STAR"))],
		DTU_pos = exprs(isoforms_cell_bank$isoform_proportions)[unlist(DTU_isoforms_positive),cells],
		DTU_pos_filt = exprs(isoforms_cell_bank$isoform_proportions)[unlist(DTU_isoforms_positive),setdiff(cells,c("COLO94H.STAR","SW1463.STAR"))]
	)

	
##PCA
	PCA_res <- lapply(x_list,function(x){prcomp(t(x))})
	names(PCA_res) <- names(x_list)
	for(i in names(PCA_res)){PCA_res[[i]]$var_explained = PCA_res[[i]]$sdev^2 / sum(PCA_res[[i]]$sdev^2)}
##tsne
	#tsne_res <- lapply(x_list,function(x){tsne(x,perplexity=20,check_duplicates=FALSE)})
	#names(tsne_res) <- names(x_list)
###out 
##scree plot

##PCA plot
	pdf("out/PCA_plots.pdf")
		lapply(names(x_list),function(x){
		barplot(PCA_res[[x]]$var_explained,xlab="PC",ylab="Percentage of\nexplained variance",main=paste("Scree plot",x,sep=" "),col="cornflowerblue")
		points(PCA_res[[x]]$var_explained,col="black",pch=20)
		lines(PCA_res[[x]]$var_explained,col="black",lwd=2)
		plot(PCA_res[[x]]$x[,1],PCA_res[[x]]$x[,2],xlab="PC1",ylab="PC2",pch=20,col=group_colors[group[rownames(PCA_res[[x]]$x)]],main=paste("PCA",x,sep=" "))	
	})
	
	dev.off()
##tsne plot
  #      pdf("out/tsne_plots.pdf")
   #     lapply(names(x_list),function(x){
#		plot(tsne_res[[x]]$Y,col=group_colors[group[cells]],xlab="tsne1",ylab="tsne2",main=paste("tsne",x,sep=" "))
#	})
#	dev.off()
	lapply(names(DTU_isoforms),function(x){
		genes <- unique(fData(isoforms_cell_bank$isoform_proportions)[DTU_isoforms[[x]],"Geneid"])
		iso <- rownames(fData(isoforms_cell_bank$isoform_proportions))[which(fData(isoforms_cell_bank$isoform_proportions)$Geneid%in%genes==TRUE)]
		tab <- 
data.frame(isoform=iso,logFC=DIE[["CRIS_adjusted"]][[x]][match(iso,DIE[["CRIS_adjusted"]][[x]]$rownam),"logFC"],
FDR=DIE[["CRIS_adjusted"]][[x]][match(iso,DIE[["CRIS_adjusted"]][[x]]$rownam),"FDR"],
SD=DIE[["CRIS_adjusted"]][[x]][match(iso,DIE[["CRIS_adjusted"]][[x]]$rownam),"proportionSD"],fData(isoforms_cell_bank$isoform_proportions)[iso,])
#		write.xlsx(tab,file="out/biotype_annotatedDTU_results.xlsx",sheetName=x,append=TRUE)
	})
	save(PCA_res, file="out/PCA_results.R")
#	save(tsne_res, file="out/tsne_results.R")

