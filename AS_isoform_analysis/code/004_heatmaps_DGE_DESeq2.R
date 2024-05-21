###library
	library(ComplexHeatmap)
	library(circlize)
	library(Biobase)

###loading
	load("out/annotated_genes_counts_cell_bank.rda")
	load("out/DGE_results_DESeq2.rda")
###function
	
	DGEHeatmap <- function(cells,group,DGE_res,title){
		log2cpm <- t(apply(cpm[,cells],1,function(x){log2((x/mean(x))+0.01)}))	
		if(is.data.frame(DGE_res)==TRUE){
			genes <- rownames(DGE_res)[which(DGE_res$FDR < 0.01 & abs(DGE_res$logFC)>1.5)]
			genes <- genes[order(DGE_res[genes,]$logFC,decreasing=T)]	
			row_split <- ifelse(DGE_res[genes,]$logFC>1,0,1)
			columns_split <- group
			row_anno <- rowAnnotation(gene_modules=anno_block(labels=paste(table(row_split),"genes",sep=" "),gp=gpar(fill=c("red","blue"))))
			col_anno <- HeatmapAnnotation(classes=anno_block(labels=levels(group),gp=gpar(fill=c("darkslategray3","darkgoldenrod1"))))
			return(Heatmap(log2cpm[genes,],cluster_rows=FALSE, cluster_columns=FALSE,name="L2R\ngene\nexpression",col=colorRamp2(c(2,0,-2),c("red","white","blue")),
top_annotation=col_anno,left_annotation=row_anno,column_title=title,row_title=NULL,show_row_names=FALSE,show_column_names=FALSE,row_split=row_split,column_split=columns_split))
		}else{
			log2cpm <- t(apply(cpm[,cells],1,function(x){log2((x/mean(x))+0.01)}))
			genes <- lapply(DGE_res,function(x){rownames(x)[which(x$FDR < 0.01 & x$logFC>2)]})
			names(genes) <- names(DGE_res)
			le <- lapply(genes,length)
			row_split <- lapply(names(DGE_res),function(x){rep(x,length(genes[[x]]))})
			row_split <- row_split[which(unlist(le)!=0)]
			columns_split <- group
			row_anno <- rowAnnotation(gene_modules=anno_block(labels=paste(le[which(unlist(le)!=0)],sep=" "),gp=gpar(fill=c("red","darkgreen","lightblue","orange","blue")[1:length(levels(group))])))
			col_anno <- 
HeatmapAnnotation(classes=anno_block(labels=sub("only","",levels(group)),labels_gp=gpar(col="white",fontsize=10),gp=gpar(fill=c("red","darkgreen","lightblue","orange","blue")[which(unlist(le)!=0)])))
			return(Heatmap(log2cpm[unlist(genes),],cluster_rows=FALSE, 
cluster_columns=FALSE,name="L2R 
\ngene\nexpression",col=colorRamp2(c(2,0,-2),c("red","white","blue")),show_column_names=FALSE,row_split=unlist(row_split),column_split=columns_split,top_annotation=col_anno,column_title=title,row_title=NULL,
show_row_names=FALSE,left_annotation=row_anno))

		}
	}

###code
	cpm <- apply(exprs(annotated_genes_counts_cell_bank),2,function(x){x*1000000/sum(x)})
	rownames(cpm) <- fData(annotated_genes_counts_cell_bank)$Geneid
	pdf("out/Heatmaps_DGE_DESeq2.pdf")

##MS status
	cells <- colnames(annotated_genes_counts_cell_bank)[which(pData(annotated_genes_counts_cell_bank)$MS_status!="NA")]
	group <- factor(pData(annotated_genes_counts_cell_bank)$MS_status[which(pData(annotated_genes_counts_cell_bank)$MS_status!="NA")])
	DGE_res <- DGE$MS_status_adjusted
	DGEHeatmap(cells,group,DGE_res,"MSI/MSS") -> p
	draw(p)

##CRIS

	group <- factor(pData(annotated_genes_counts_cell_bank)$CRIS_class[grep("only",pData(annotated_genes_counts_cell_bank)$CRIS_class)])
	cells <- colnames(annotated_genes_counts_cell_bank)[grep("only",pData(annotated_genes_counts_cell_bank)$CRIS_class)]
	DGE_res <- DGE$CRIS_adjusted
	DGEHeatmap(cells,group,DGE_res,"CRIS") -> p
	draw(p)

##CMS

	group <- factor(pData(annotated_genes_counts_cell_bank)$CMS_class[which(pData(annotated_genes_counts_cell_bank)$CMS_class!="NA")])
	cells <- colnames(annotated_genes_counts_cell_bank)[which(pData(annotated_genes_counts_cell_bank)$CMS_class!="NA")]
	DGE_res <- DGE$CMS_adjusted
	DGEHeatmap(cells,group,DGE_res,"CMS") -> p
	draw(p)

## cetuximab DGE

	group <- factor(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity[which(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
	cells <- colnames(annotated_genes_counts_cell_bank)[which(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity%in%c(0,1)==TRUE)]
	DGE_res <- DGE$Cetuximab_adjusted
	DGEHeatmap(cells,group,DGE_res,"Cetuximab\nsensitivity") -> p
	draw(p)
## mutations 

	pData(annotated_genes_counts_cell_bank)$mutations[which(pData(annotated_genes_counts_cell_bank)$mutations=="NRAS")] <- NA
        group <- factor(pData(annotated_genes_counts_cell_bank)$mutations[which(pData(annotated_genes_counts_cell_bank)$mutations!="NA")])
	cells <- colnames(annotated_genes_counts_cell_bank)[which(pData(annotated_genes_counts_cell_bank)$mutations!="NA")]
	DGE_res <- DGE$mutations_adjusted
	DGEHeatmap(cells,group,DGE_res,"Mutations") -> p
	draw(p)

## CIMP

        sel <- pData(annotated_genes_counts_cell_bank)$CIMP_status[which(pData(annotated_genes_counts_cell_bank)$CIMP_status!="0")]
        group <- factor(ifelse(sel=="CIMP-H","CIMP-H","rest"))
	cells <- colnames(annotated_genes_counts_cell_bank)[which(pData(annotated_genes_counts_cell_bank)$CIMP_status!="0")]
	DGE_res <- DGE$CIMP_adjusted
	DGEHeatmap(cells,group,DGE_res,"CIMP-H/CIMP-L") -> p
	draw(p)

dev.off()
