###library
	library(circlize)
	library(ComplexHeatmap)

###loading
	load("data/isoforms_cell_bank.rda")
	load("out/CRIS_class.rda")
###code

	tab <- exprs(isoforms_cell_bank$TPM)[which(fData(isoforms_cell_bank$TPM)$Geneid=="KRAS"),]
	su <- apply(tab,2,sum)
	perc <- t(apply(tab,1,function(x){x*100/su}))
	perc[is.na(perc)] <- 0
	col <- rep("grey",length(cla))
	names(col) <- cla
	col["onlyCRISA"] <- "red"
	col["onlyCRISB"] <- "darkgreen"
	col["onlyCRISC"] <- "lightblue"
	col["onlyCRISD"] <- "orange"
	col["onlyCRISE"] <- "blue"
	pdf("out/KRAS_iso.pdf",width=16,height=12)
	
draw(Heatmap(perc[c("ENST00000256078.10","ENST00000311936.8","ENST00000685328.1","ENST00000688940.1"),names(cla)],cluster_rows=FALSE,
top_annotation=HeatmapAnnotation(CRIS_CLASS=anno_simple(cla,col=col)),col=colorRamp2(c(0,100),c("blue","red")),show_row_names=FALSE,split=c(1,2,2,2),row_title=NULL,left_annotation=rowAnnotation(KRAS=anno_block(labels=c("KRAS2A","KRAS2B"),gp=gpar(fill=c("yellow","red")))),name="% 
Abundance"))

	dev.off()
