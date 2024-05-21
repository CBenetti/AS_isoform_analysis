###library
	library(Biobase)
###loading
	load("out/annotated_isoforms_cell_bank_proportions.rda")
	load("out/annotated_genes_counts_cell_bank.rda")
	load("out/annotated_genes_counts_cell_bank.rda")
	load("out/DIE_results_DRIMSeq_withL2R.rda")	
	load("out/DGE_results.rda")
	read.xlsx("data/CRIS_genes.xlsx",sheetIndex=1)-> CRIS_genes
###code
	DGE_genes_positive <- lapply(DGE[["CRIS_adjusted"]],function(x){rownames(x$table)[which(x$table$FDR < 0.01 & x$table$logFC > 2)]})
	dup <- unlist(DGE_genes_positive)[which(duplicated(unlist(DGE_genes_positive))==TRUE)]
	lapply(DGE_genes_positive,function(x){setdiff(x,dup)})-> DGE_genes_positive
        DTU_isoforms_positive <- lapply(DIE[["CRIS_adjusted"]],function(x){x$rownam[which(x$FDR < 0.01 & x$logFC >0 & x$proportionSD > 0.1)]})
	DGE_length <- lapply(DGE_genes_positive,length)
	DTU_length <- lapply(DTU_isoforms_positive,length)
	class_specificity_genes <- unlist(lapply(names(DGE_length),function(x){rep(x,DGE_length[[x]])}))
	class_specificity_isoforms <- unlist(lapply(names(DTU_length),function(x){rep(x,DTU_length[[x]])}))	
	rownames(fData(isoforms_cell_bank$isoform_proportions)) <- rownames(exprs(isoforms_cell_bank$isoform_proportions))
	exprs(annotated_genes_counts_cell_bank) <- apply(exprs(annotated_genes_counts_cell_bank),2,function(x){(x*1000000)/sum(x)})
	CRIS_genes <- CRIS_genes[which(CRIS_genes$Gene.ID%in%fData(annotated_genes_counts_cell_bank)$Geneid==TRUE),]
	signature_genes <- 
annotated_genes_counts_cell_bank[which(fData(annotated_genes_counts_cell_bank)$Geneid%in%unlist(CRIS_genes$Gene.ID)==TRUE),grep("onlyCRIS",pData(annotated_genes_counts_cell_bank)$CRIS_class)]
	predictive_genes <- 
annotated_genes_counts_cell_bank[which(fData(annotated_genes_counts_cell_bank)$Geneid%in%unlist(DGE_genes_positive)==TRUE),grep("onlyCRIS",pData(annotated_genes_counts_cell_bank)$CRIS_class)]	
	predictive_isoform_proportions <- isoforms_cell_bank$isoform_proportions[unlist(DTU_isoforms_positive),grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)]
	fData(predictive_genes)$class_specificity <- class_specificity_genes
	fData(signature_genes)$class_specificity <- paste("only",gsub("-","",CRIS_genes$CRIS.Class),sep="")
	fData(predictive_isoform_proportions)$class_specificity <- class_specificity_isoforms

###out
	save(signature_genes,file="out/signature_genes.rda")
	save(predictive_genes,file="out/predictive_genes.rda")
	save(predictive_isoform_proportions,file="out/predictive_isoform_proportions.rda")

