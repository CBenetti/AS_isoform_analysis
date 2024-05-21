###library
	library(ComplexHeatmap)
	library(circlize)
	library(Biobase)

###loading
	load("out/annotated_isoforms_cell_bank_proportions.rda")
	DTU_table <- as.data.frame(read.table("out/DTU_table.txt",sep="\t",header=T,row.names=1))
	load("out/DGE_signif_res.rda")
###function
	
	DTUHeatmap <- function(cells,group,cla,title){
		log2proportion <- t(apply(exprs(isoforms_cell_bank$isoform_proportions)[,cells],1,function(x){log2((x/mean(x))+0.001)}))	
		if(length(levels(group))<3){
			transcripts <- rownames(DTU_table)[which(DTU_table[,cla]!=0)]
			transcripts <- transcripts[order(DTU_table[transcripts,cla],decreasing=T)]	
			row_split <- ifelse(DTU_table[transcripts,cla]>0,0,1)
			columns_split <- group
			is_signif_gene <- DGE_signif_res[transcripts,cla]
			row_anno_left <- rowAnnotation(transcript_n=anno_block(labels=paste(table(row_split),"isoforms",sep=" "),gp=gpar(fill=c("red","blue"))))
			row_anno_right <- rowAnnotation(DGE_results=as.character(is_signif_gene),col=list(DGE_results=c("-1"="blue","0"="grey","1"="red")))
			col_anno <- HeatmapAnnotation(classes=anno_block(labels=levels(group),gp=gpar(fill=c("darkslategray3","darkgoldenrod1"))))
			return(Heatmap(log2proportion[transcripts,],cluster_rows=FALSE, 
cluster_columns=FALSE,name="L2R\nisoform\nproportion",col=colorRamp2(c(1,0,-1),c("red","white","blue")),
top_annotation=col_anno,left_annotation=row_anno_left,right_annotation=row_anno_right,column_title=title,row_title=NULL,show_row_names=FALSE,show_column_names=FALSE,
row_split=row_split,column_split=columns_split))
		}else{
			transcripts <- lapply(unique(group)[which(unique(group)%in%colnames(DTU_table))],function(x){rownames(DTU_table)[which(DTU_table[,x]>0)]})
			names(transcripts) <- unique(group)[which(unique(group)%in%colnames(DTU_table))]
			le <- lapply(transcripts,length)
			row_split <- lapply(unique(group)[which(unique(group)%in%colnames(DTU_table))],function(x){rep(x,length(transcripts[[x]]))})
			columns_split <- group
			is_signif_gene <- unlist(lapply(names(transcripts),function(x){DGE_signif_res[transcripts[[x]],x]}))
			row_anno_left <- rowAnnotation(transcript_n=anno_block(labels=paste(le,sep=" "),
gp=gpar(fill=c("red","darkgreen","lightblue","orange","blue")[1:length(levels(group))])))
			row_anno_right <- rowAnnotation(DGE_results=as.character(is_signif_gene),col=list(DGE_results=c("-1"="blue","0"="grey","1"="red")))		
col_anno <- HeatmapAnnotation(classes=anno_block(labels=sub("only","",levels(group)),labels_gp=gpar(col="white",fontsize=10)
,gp=gpar(fill=c("red","darkgreen","lightblue","orange","blue")[1:length(levels(group))])))
			return(Heatmap(log2proportion[unlist(transcripts),],cluster_rows=FALSE, 
cluster_columns=FALSE,name="L2R\nisoform\nproportion",col=colorRamp2(c(1,0,-1),c("red","white","blue")),
show_column_names=FALSE,row_split=unlist(row_split),column_split=columns_split,top_annotation=col_anno,column_title=title,row_title=NULL,
show_row_names=FALSE,left_annotation=row_anno_left,right_annotation=row_anno_right))

		}
	}

###code
	pdf("out/Heatmaps_DTU_DRIMSeq.pdf")

##MS status
	cells <- colnames(isoforms_cell_bank$isoform_proportions)[which(pData(isoforms_cell_bank$isoform_proportions)$MS_status!="NA")]
	group <- factor(pData(isoforms_cell_bank$isoform_proportions)$MS_status[which(pData(isoforms_cell_bank$isoform_proportions)$MS_status!="NA")])
	DTUHeatmap(cells,group,"MS_status_adjusted","MSI/MSS") -> p
	draw(p)

##CRIS

	group <- factor(pData(isoforms_cell_bank$isoform_proportions)$CRIS_class[grep("only",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)])
	cells <- colnames(isoforms_cell_bank$isoform_proportions)[grep("only",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)]
	DTUHeatmap(cells,group,"NA","CRIS") -> p
	draw(p)

##CMS

	group <- factor(pData(isoforms_cell_bank$isoform_proportions)$CMS_class[which(pData(isoforms_cell_bank$isoform_proportions)$CMS_class!="NA")])
	cells <- colnames(isoforms_cell_bank$isoform_proportions)[which(pData(isoforms_cell_bank$isoform_proportions)$CMS_class!="NA")]
	DTUHeatmap(cells,group,"NA","CMS") -> p
	draw(p)

## cetuximab DGE

	group <- factor(pData(isoforms_cell_bank$isoform_proportions)$cetuximab_sensitivity[which(pData(isoforms_cell_bank$isoform_proportions)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
	cells <- colnames(isoforms_cell_bank$isoform_proportions)[which(pData(isoforms_cell_bank$isoform_proportions)$cetuximab_sensitivity%in%c(0,1)==TRUE)]
	DTUHeatmap(cells,group,"Cetuximab_adjusted","Cetuximab\nsensitivity") -> p
	draw(p)
## mutations 

	pData(isoforms_cell_bank$isoform_proportions)$mutations[which(pData(isoforms_cell_bank$isoform_proportions)$mutations=="NRAS")] <- NA
        group <- factor(pData(isoforms_cell_bank$isoform_proportions)$mutations[which(pData(isoforms_cell_bank$isoform_proportions)$mutations!="NA")])
	cells <- colnames(isoforms_cell_bank$isoform_proportions)[which(pData(isoforms_cell_bank$isoform_proportions)$mutations!="NA")]
	DTUHeatmap(cells,group,"mutations_adjusted","Mutations") -> p
	draw(p)

## CIMP

        sel <- pData(isoforms_cell_bank$isoform_proportions)$CIMP_status[which(pData(isoforms_cell_bank$isoform_proportions)$CIMP_status!="0")]
        group <- factor(ifelse(sel=="CIMP-H","CIMP-H","rest"))
	cells <- colnames(isoforms_cell_bank$isoform_proportions)[which(pData(isoforms_cell_bank$isoform_proportions)$CIMP_status!="0")]
	DTUHeatmap(cells,group,"CIMP_adjusted","CIMP-H/CIMP-L") -> p
	draw(p)

dev.off()
