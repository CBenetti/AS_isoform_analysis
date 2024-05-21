###library
	library(ComplexHeatmap)
	library(circlize)
	library(Biobase)

###loading
	load("out/annotated_isoforms_cell_bank_counts.rda")
	load("out/DGE_signif_res.rda")
	load("out/DTE_results.rda")
###function
	
	DTEHeatmap <- function(cells,group,DTE_res,title){
	##our normal way
		log2cpm <- t(apply(cpm[,cells],1,function(x){log2((x/mean(x))+0.01)}))
		
	##What edgeR suggests
		#x <- exprs(isoforms_cell_bank$counts_75)[,cells]
                #rownames(x)<- fData(isoforms_cell_bank$counts_75)$Geneid
                #data <- DGEList(counts=x,group=group)
                #keep <- filterByExpr(data)
                #filt_data <- data[keep,,keep.lib.sizes=FALSE]
                #filt_data <- normLibSizes(filt_data)
		#logcpm <- cpm(filt_data, log=TRUE)
		#log2cpm <- t(apply(logcpm,1,function(x){x-mean(x)}))	
		if(is.data.frame(DTE_res)==TRUE){
			genes <- rownames(DTE_res)[which(DTE_res$FDR < 0.01 & abs(DTE_res$logFC)>2)]
			genes <- genes[order(DTE_res[genes,]$logFC,decreasing=T)]	
			row_split <- ifelse(DTE_res[genes,]$logFC>1,0,1)
			columns_split <- group
			row_anno <- rowAnnotation(gene_modules=anno_block(labels=paste(table(row_split),"isoforms",sep=" "),gp=gpar(fill=c("red","blue"))))
			col_anno <- HeatmapAnnotation(classes=anno_block(labels=levels(group),gp=gpar(fill=c("darkslategray3","darkgoldenrod1"))))
			right <- rowAnnotation(DGE_results=as.character(DTE_res[genes,"DGE_res"]),col=list(DGE_results=c("-1"="blue","0"="white","1"="red")))			
return(Heatmap(log2cpm[genes,],cluster_rows=FALSE,cluster_columns=FALSE,name="L2R\nisoform\nexpression",col=colorRamp2(c(3,0,-3),c("red","white","blue")),
top_annotation=col_anno,left_annotation=row_anno,right_annotation=right,column_title=title,row_title=NULL,show_row_names=FALSE,show_column_names=FALSE,row_split=row_split,column_split=columns_split))
		}else{
			genes <- lapply(DTE_res,function(x){rownames(x$table)[which(x$table$FDR < 0.05 & x$table$logFC>2)]})
			names(genes) <- names(DTE_res)
			le <- lapply(genes,length)
			#if(identical(unique(unlist(le)),0)==FALSE){
			row_split <- lapply(names(DTE_res),function(x){rep(x,length(genes[[x]]))})
			row_split <- row_split[which(unlist(le)!=0)]
			columns_split <- group
			row_anno <- rowAnnotation(gene_modules=anno_block(labels=paste(le[which(unlist(le)!=0)],sep=" "),gp=gpar(fill=c("red","darkgreen","lightblue","orange","blue")[which(unlist(le)!=0)])))
			col_anno <- 
HeatmapAnnotation(classes=anno_block(labels=sub("only","",levels(group)),labels_gp=gpar(col="white",fontsize=10),gp=gpar(fill=c("red","darkgreen","lightblue","orange","blue")[1:length(levels(group))])))
			DGE_res <- unlist(lapply(names(DTE_res),function(x){DTE_res[[x]]$table[genes[[x]],"DGE_res"]}))
			right <- rowAnnotation(DGE_results=as.character(DGE_res),col=list(DGE_results=c("-1"="blue","0"="white","1"="red")))
			return(Heatmap(log2cpm[unlist(genes),],cluster_rows=FALSE, 
cluster_columns=FALSE,name="L2R 
\nisoform\nexpression",col=colorRamp2(c(3,0,-3),c("red","white","blue")),show_column_names=FALSE,row_split=unlist(row_split),column_split=columns_split,top_annotation=col_anno,column_title=title,row_title=NULL,
show_row_names=FALSE,left_annotation=row_anno,right_annotation=right))
#		}
		}
	}

###code
	cpm <- apply(exprs(isoforms_cell_bank$counts_75),2,function(x){x*1000000/sum(x)})
	rownames(cpm) <- rownames(exprs(isoforms_cell_bank$counts_75))
	pdf("out/Heatmaps_DTE.pdf")

##MS status
	cells <- colnames(isoforms_cell_bank$counts_75)[which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")]
	group <- factor(pData(isoforms_cell_bank$counts_75)$MS_status[which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")])
	DTE_res <- DTE$MS_status_adjusted$table
	DTEHeatmap(cells,group,DTE_res,"MSI/MSS") -> p
	draw(p)

##CRIS

	group <- factor(pData(isoforms_cell_bank$counts_75)$CRIS_class[grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)])
	cells <- colnames(isoforms_cell_bank$counts_75)[grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)]
	DTE_res <- DTE$CRIS_adjusted
	DTEHeatmap(cells,group,DTE_res,"CRIS") -> p
	draw(p)

##CMS

	group <- factor(pData(isoforms_cell_bank$counts_75)$CMS_class[which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")])
	cells <- colnames(isoforms_cell_bank$counts_75)[which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")]
	DTE_res <- DTE$CMS_adjusted
#	DTEHeatmap(cells,group,DTE_res,"CMS") -> p
#	draw(p)

## cetuximab DTE

	group <- factor(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity[which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
	cells <- colnames(isoforms_cell_bank$counts_75)[which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)]
	DTE_res <- DTE$Cetuximab_adjusted$table
	DTEHeatmap(cells,group,DTE_res,"Cetuximab\nsensitivity") -> p
	draw(p)
## mutations 

	pData(isoforms_cell_bank$counts_75)$mutations[which(pData(isoforms_cell_bank$counts_75)$mutations=="NRAS")] <- NA
        group <- factor(pData(isoforms_cell_bank$counts_75)$mutations[which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")])
	cells <- colnames(isoforms_cell_bank$counts_75)[which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")]
	DTE_res <- DTE$mutations_adjusted
	DTEHeatmap(cells,group,DTE_res,"Mutations") -> p
	draw(p)

## CIMP

        sel <- pData(isoforms_cell_bank$counts_75)$CIMP_status[which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")]
        group <- factor(ifelse(sel=="CIMP-H","CIMP-H","rest"))
	cells <- colnames(isoforms_cell_bank$counts_75)[which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")]
	DTE_res <- DTE$CIMP_adjusted$table
	DTEHeatmap(cells,group,DTE_res,"CIMP-H/CIMP-L") -> p
	draw(p)

dev.off()
