###library
	library(circlize)
	library(ComplexHeatmap)

###loading
	load("out/annotated_isoforms_cell_bank_proportions.rda")
#	load("out/CRIS_class.rda")
###code
group <- factor(pData(isoforms_cell_bank$isoform_proportions)$CRIS_class[grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)])
sub <- exprs(isoforms_cell_bank$isoform_proportions)[,grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)]	
ord_sub <- lapply(levels(group),function(x){s=colnames(sub)[which(group==x)];s[order(sub["ENST00000256078.10",which(group==x)])]})
co <- c("white","grey","coral4","olivedrab3")
co1 <- c("white","grey","cyan4")
co2 <- c("red","darkgreen","lightblue","orange","blue")
co3 <- c("#FFFF33","#CCFFFF","grey")
names(co3) <- c("MSI","MSS","NA")
names(co2) <- levels(group)
names(co1) <- c(1,"NA",0)
names(co) <- c("WT","NA","KRAS","BRAF")
	pdf("out/KRAS_iso.pdf",width=16,height=12)
draw(Heatmap(sub[c("ENST00000256078.10","ENST00000311936.8","ENST00000685328.1","ENST00000688940.1"),],cluster_rows=FALSE,
column_split=group, column_order=unlist(ord_sub),column_title=NULL,
top_annotation=HeatmapAnnotation(CRIS_CLASS=anno_block(labels=levels(group),gp=gpar(fill=c("red","darkgreen","lightblue","orange","blue")),,labels_gp= 
gpar(col="white")),mutation=anno_simple(pData(isoforms_cell_bank$isoform_proportions)$mutations[grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)],col=co),
cetuximab=anno_simple(pData(isoforms_cell_bank$isoform_proportions)$cetuximab_sensitivity[grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)],col=co1),
MSI_MSS=anno_simple(pData(isoforms_cell_bank$isoform_proportions)$MS_status[grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)],col=co3)),
col=colorRamp2(c(0,1),c("blue","red")),show_row_names=FALSE,split=c(1,2,2,2),row_title=NULL,
left_annotation=rowAnnotation(KRAS=anno_block(labels=c("KRAS4A","KRAS4B"),gp=gpar(fill=c("yellow","red")))),name="% Abundance"))

        dev.off()
ord_sub <- lapply(levels(group),function(x){s=colnames(sub)[which(group==x)];s[order(sub["ENST00000379666.7",which(group==x)])]})
	pdf("out/MRLP33_iso.pdf",width=16,height=12)
draw(Heatmap(sub[c("ENST00000379666.7","ENST00000296102.8"),],cluster_rows=FALSE,
top_annotation=HeatmapAnnotation(CRIS_CLASS=anno_simple(pData(isoforms_cell_bank$isoform_proportions)$CRIS_class[grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)],col=co2),
mutation=anno_simple(pData(isoforms_cell_bank$isoform_proportions)$mutations[grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)],col=co),
cetuximab=anno_simple(pData(isoforms_cell_bank$isoform_proportions)$cetuximab_sensitivity[grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)],col=co1),
MSI_MSS=anno_simple(pData(isoforms_cell_bank$isoform_proportions)$MS_status[grep("onlyCRIS",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)],col=co3)),
col=colorRamp2(c(0,1),c("blue","red")),show_row_names=FALSE,split=c(1,2),row_title=NULL,column_split = 4,
left_annotation=rowAnnotation(KRAS=anno_block(labels=c("MRPL33-tr","MRPL33-fl"),gp=gpar(fill=c("yellow","red")))),name="% Abundance"))

        dev.off()

