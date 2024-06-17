###library
library(xlsx)
library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)
###loading
load("data/count_raw_results_cutoff_0.05.rda")
load("../AS_isoform_analysis/out/annotated_genes_counts_cell_bank.rda")
load("../AS_isoform_analysis/out/annotated_isoforms_cell_bank_proportions.rda")
read.xlsx("primers/NOEMI RUSSO_140624_Results.xlsx",sheetName="Mean",header=T,)->PCR
rownames(PCR) <- PCR[,1]
PCR <- PCR[,-1]
###code
iso <- c()
introns <-  c(
  "chr14:24058986:24059134:+:CARMIL3",
  "chr17:64500326:64500548:-:DDX5",
  "chr1:211365468:211369451:+:TRAF5",
  rownames(raw_tab)[4192],
  "chr16:31185179:31186801:+:FUS",
  "chr20:2862934:2863064:+:VPS16"
)

introns_tpm <- raw_tab[introns,c("LS180.STAR","COGA2.STAR","HROC57.STAR","SUN1047.STAR")]
introns_PCR <- apply(PCR,2,function(x){2^x})
introns_PCR_norm <- do.call("rbind",
lapply(rownames(introns_PCR)[1:6],function(x){introns_PCR[x,]*(mean(cpm[which(fData(annotated_genes_counts_cell_bank)$Geneid==x),colnames(introns_tpm)])/
cpm[which(fData(annotated_genes_counts_cell_bank)$Geneid==x),colnames(introns_tpm)])}))
rownames(introns_PCR_norm) <- rownames(introns_PCR)[1:6]
genes_cpm <- do.call("rbind",lapply(rownames(introns_PCR)[1:6],function(x){cpm[which(fData(annotated_genes_counts_cell_bank)$Geneid==x),colnames(introns_tpm)]}))
rownames(genes_cpm) <- rownames(introns_PCR)[1:6]

PCR_proportions <- rbind(c(0.049451501,0.270770882,0.51623765),c(0.950548499,0.729229118,0.48376235))
isoform_TPM <- exprs(isoforms_cell_bank$TPM)[c("ENST00000379666.7","ENST00000296102.8"),c("WIDR.STAR","HROC57.STAR","HDC8.STAR")]
###out
pdf("out/PCR_results.pdf")
par(mfrow=c(2,3),xpd=TRUE)
	for(x in rownames(introns_PCR)[1:6]){
	barplot(introns_PCR[x,c("LS.180","COGA2","HROC57","SNU1047")],ylab="2^(40-normCt)",col=brewer.pal(n=4,"Paired"),main=paste("PCR",x,sep=" "),las=2, ylim=c(0,max(introns_PCR[x,])))
	barplot(introns_tpm[grep(x,rownames(introns_tpm)),],ylab="TPM",col=brewer.pal(n=4,"Paired"),main=paste("Stringtie",x,sep=" "),las=2,
 ylim=c(0,max(introns_tpm[grep(x,rownames(introns_tpm)),])))
	barplot(genes_cpm[x,],ylab="CPM",col=brewer.pal(n=4,"Paired"),main=paste("Gene expression",x,sep=" "),las=2, ylim=c(0, max(genes_cpm[x,])))
	#barplot(introns_PCR_norm[x,],ylab="2^(40-normCt)",col=brewer.pal(n=4,"Paired"),main=paste("PCR norm",x,sep=" "))
}
par(mfrow=c(1,1))
ht_opt(heatmap_border = TRUE)
ht_opt$COLUMN_ANNO_PADDING=unit(0.3,"cm")
draw(Heatmap(rbind(PCR_proportions,exprs(isoforms_cell_bank$isoform_proportions)[c("ENST00000379666.7","ENST00000296102.8"),c("WIDR.STAR","HROC57.STAR","HDC8.STAR")]),name="%",
colorRamp2(c(0,1),c("blue","red")),show_row_names=FALSE,column_split=c(1,2,3),split=c(1,2,3,4),row_title=NULL,cluster_rows=FALSE,cluster_columns=FALSE,show_column_names=FALSE,
column_title=c("WIDR","HROC57","HDC8"),left_annotation=rowAnnotation(iso=anno_block(labels=c("MRPL33-tr","MRPL33-fl","MRPL33-tr","MRPL33-fl"),
gp=gpar(fill=c("yellow","red","yellow","red")),labels_gp=gpar(fontsize=6))),top_annotation=HeatmapAnnotation(barplot=anno_empty(height = unit(4, "cm")),barplot2=anno_empty(height = unit(4, 
"cm")))),
padding = unit(c(10, 20, 20, 20), "mm"),height = unit(5, "cm"))

#decorate_heatmap_body("%",{
#    grid.text("StringTie isoform proportions    PCR proportion",x = unit(-1.25, "cm"),rot = 90, just = c("center","bottom"))})
decorate_annotation("barplot", {
    pushViewport(viewport(xscale = c(0.45, 2.45), yscale = c(0, max(introns_PCR[c(7,8),c("WIDR","HROC57","HDC8")])*1.1)))
    grid.rect(x = 1:2, y = 0, width = 0.7, height = introns_PCR[c(7,8),"WIDR"], just = "bottom",
    gp = gpar(fill = c("yellow","red")), default.units = "native")
    grid.yaxis()
    grid.text("PCR", x = unit(-2.20, "cm"),rot = 90, just = "bottom")
    grid.text("2^(40-normCt)", x = unit(-1.75, "cm"),rot = 90, just = "bottom")
    popViewport()
})
decorate_annotation("barplot", {
    pushViewport(viewport(xscale = c(0.45, 2.45), yscale = c(0, max(introns_PCR[c(7,8),c("WIDR","HROC57","HDC8")])*1.1)))
    grid.rect(x = 3:4, y = 0, width = 0.7, height = introns_PCR[c(7,8),"HROC57"], just = "bottom",
    gp = gpar(fill = c("yellow","red")), default.units = "native")
    grid.yaxis()
    popViewport()
})
decorate_annotation("barplot", {
    pushViewport(viewport(xscale = c(0.45, 2.45), yscale = c(0, max(introns_PCR[c(7,8),c("WIDR","HROC57","HDC8")])*1.1)))
    grid.rect(x = 5:6, y = 0, width = 0.7, height = introns_PCR[c(7,8),"HDC8"], just = "bottom",
    gp = gpar(fill = c("yellow","red")), default.units = "native")
    grid.yaxis()
    popViewport()
})

decorate_annotation("barplot2", {
    pushViewport(viewport(xscale = c(0.45, 2.45), yscale = c(0, max(isoform_TPM)*1.1)))
    grid.rect(x = 1:2, y = 0, width = 0.7, height = isoform_TPM[,"WIDR.STAR"], just = "bottom",
    gp = gpar(fill = c("yellow","red")), default.units = "native")
    grid.yaxis()
    grid.text("StringTie", x = unit(-2, "cm"),rot = 90, just = "bottom")
    grid.text("TPM", x = unit(-1.75, "cm"),rot = 90, just = "bottom")
    popViewport()
})

decorate_annotation("barplot2", {
    pushViewport(viewport(xscale = c(0.45, 2.45), yscale = c(0, max(isoform_TPM)*1.1)))
    grid.rect(x = 3:4, y = 0, width = 0.7, height = isoform_TPM[,"HROC57.STAR"], just = "bottom",
    gp = gpar(fill = c("yellow","red")), default.units = "native")
    grid.yaxis()
    popViewport()
})
decorate_annotation("barplot2", {
    pushViewport(viewport(xscale = c(0.45, 2.45), yscale = c(0, max(isoform_TPM)*1.1)))
    grid.rect(x = 5:6, y = 0, width = 0.7, height = isoform_TPM[,"HDC8.STAR"], just = "bottom",
    gp = gpar(fill = c("yellow","red")), default.units = "native")
    grid.yaxis()
    popViewport()
})
dev.off()
