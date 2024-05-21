#c##library
	library(xlsx)
	library(Biobase)
	library(ggvenn)
	library(eulerr)
	library(cowplot)
	library(ComplexHeatmap)
	library(circlize)
###loading
        load("data/isoforms_cell_bank.rda")
        load("data/CRIS_class.rda")
	anno <- read.xlsx("data/Annotations.xlsx",sheetIndex=1)
	s_r <- read.xlsx("data/Sample_replicates.xlsx",sheetIndex=1,header=F)
###code
	pdf("out/Fig1A_annotations_in_cell_bank.pdf")
	total_anno <- list(Cell_bank=colnames(isoforms_cell_bank$TPM),CRIS=names(cla))
	set.seed(444)
	fit3 <- euler(total_anno,shape="ellipse")
	plot(fit3,fills = list(fill = c("#0073C2FF", "#EFC000FF","#CD534CFF"), alpha = .5),quantities = list(cex = .8),legend=T,cex.lab=1, main=paste("Cell Bank:",length(colnames(isoforms_cell_bank$TPM)),"total cell lines")) -> p
	for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
  	o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
   	if(!is.null(o)){
     	  if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
       		o$children[[paste0("tag.quantity.",i)]]$label <- " "
       		p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
     	  }
	 }
	}
	plot(plot_grid(p,labels="A",label_size=18))
	dev.off()

	mat <- exprs(isoforms_cell_bank$TPM)
##filtering CRIS_samples
	isoforms_cell_bank$TPM <- isoforms_cell_bank$TPM[,names(cla)]
	isoforms_cell_bank$cov <- isoforms_cell_bank$cov[,names(cla)]
	isoforms_cell_bank$FPKM <- isoforms_cell_bank$FPKM[,names(cla)]
	
##sample/replicate correspondance
	anno$LIM1215S22 <- anno$LIM1215
	anno$HDC114S17 <- anno$HDC114
	anno$DLD1S12 <- anno$DLD1
	anno$HROC69S8 <- anno$HROC69
	anno$HROC24S4 <- anno$HROC24
	anno$SW620S14 <- anno$SW620
	anno$SNU1235S21 <- anno$SNU1235
	anno$SNU1040S11 <- anno$SNU1040
	anno$SNU81S15 <- anno$SNU81
	anno$LS411NS23 <- anno$LS411N
	anno$LS180S23 <- anno$LS180
	anno$C10S1 <- anno$C10
	anno$COLO3202 <- anno$COLO320
	anno$COGA122 <- anno$COGA12
	anno$COLO94H2 <- anno$COLO94H
	#anno$HROC24NT <- anno$HROC24
	
##typo correction
	colnames(anno)[which(colnames(anno)=="OUMS23")] <- "OUSM23"
	colnames(anno)[which(colnames(anno)=="IRCC.1")] <- "IRCC1"
	colnames(anno)[which(colnames(anno)=="HDC9")] <- "HCD9"
	colnames(anno)[which(colnames(anno)=="NCIH716")] <- "NCHI716"

##not_annotated_samples
	not_in_anno <- colnames(isoforms_cell_bank$TPM)[which(colnames(isoforms_cell_bank$TPM)%in%paste(colnames(anno),".STAR",sep="")==FALSE)]
	not_in_CRIS <- colnames(isoforms_cell_bank$TPM)[which(colnames(isoforms_cell_bank$TPM)%in%names(cla)==FALSE)]	
	not_in_both <- names(table(c(not_in_anno,not_in_CRIS)))[which(table(c(not_in_anno,not_in_CRIS))==2)]
	anno_sub <- anno[,which(paste(colnames(anno),".STAR",sep="")%in%names(cla)==TRUE)]
	rownames(anno_sub) <- anno[,1]
	colnames(anno_sub) <- paste(colnames(anno_sub),".STAR",sep="")
	anno_sub <- as.matrix(anno_sub)
##annotation

	pdata <- 
data.frame(source=rep("NA",length(colnames(isoforms_cell_bank$TPM))),organism=rep("Homo Sapiens",length(colnames(isoforms_cell_bank$TPM))),
cetuximab_sensitivity=rep("NA",length(colnames(isoforms_cell_bank$TPM))),CIMP_status=rep(0,length(colnames(isoforms_cell_bank$TPM))),
mutations = rep("NA",length(colnames(isoforms_cell_bank$TPM))),MS_status=rep("NA",length(colnames(isoforms_cell_bank$TPM))),
CMS_class=rep("NA",length(colnames(isoforms_cell_bank$TPM))),CRIS_class =rep("NA",length(colnames(isoforms_cell_bank$TPM))))
	rownames(pdata) <- colnames(isoforms_cell_bank$TPM)


	pdata[colnames(anno_sub),1] <- as.vector(anno_sub["source name",])
	pdata[colnames(anno_sub),3] <- as.vector(anno_sub["Cetuximab SENS",])
	pdata[colnames(anno_sub),4] <- as.vector(anno_sub["CIMP status",])
	pdata[colnames(anno_sub),5] <- as.vector(anno_sub["KRAS/BRAF status",])
	pdata[colnames(anno_sub),6] <- as.vector(anno_sub["MSI/MSS",])
	pdata[colnames(anno_sub),7] <- as.vector(anno_sub["CMS",])
	pdata[names(cla),8] <- cla	
	pData(isoforms_cell_bank$TPM) <- pData(isoforms_cell_bank$FPKM) <- pData(isoforms_cell_bank$cov) <- as.data.frame(pdata)

##graphical representation
	pdf("out/Fig1A-F_annotation_within_classes.pdf",width=12,height=18)
	par(mfrow=c(3,2))
	cetuximab <- list(Cell_bank= colnames(isoforms_cell_bank$TPM),Sensitive=colnames(anno_sub)[which(anno_sub["Cetuximab SENS",]==1)],Resistant=colnames(anno_sub)[which(anno_sub["Cetuximab SENS",]==0)])
	fit4 <-euler(cetuximab,shape="ellipse")
	plot(fit4,fills = list(fill = c("seashell2", "turquoise4","tomato3"), alpha = 0.5),quantities = list(cex = 2),legend=T,cex.lab=1, main="Cetuximab \nsensitivity annotations") -> p
	for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
  	o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
   	if(!is.null(o)){
     	  if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
       		o$children[[paste0("tag.quantity.",i)]]$label <- " "
       		p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
     	  }
	 }
	}
	p -> p1
	MS_status <- list(Cell_bank= colnames(isoforms_cell_bank$TPM),MSI=colnames(anno_sub)[which(anno_sub["MSI/MSS",]=="MSI")],MSS=colnames(anno_sub)[which(anno_sub["MSI/MSS",]=="MSS")])
	set.seed(444)
	fit4 <-euler(MS_status,shape="ellipse")
	plot(fit4,fills = list(fill = c("seashell2", "gold","lightskyblue"), alpha = 0.5),quantities = list(cex = 2),legend=T,cex.lab=1, main="Microsatellite\n instability annotations") -> p
	for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
  	o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
   	if(!is.null(o)){
     	  if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
       		o$children[[paste0("tag.quantity.",i)]]$label <- " "
       		p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
     	  }
	 }
	}
	p -> p2	

	CIMP_status <- list(Cell_bank= colnames(isoforms_cell_bank$TPM),CIMP_H=colnames(anno_sub)[which(anno_sub["CIMP status",]=="CIMP-H")],CIMP_L=colnames(anno_sub)[which(anno_sub["CIMP status",]%in%c(0,"CIMP-H")==FALSE)])
	set.seed(444)
	fit4 <-euler(CIMP_status,shape="ellipse")
	plot(fit4,fills = list(fill = c("seashell2", "olivedrab4","royalblue2"), alpha = 0.5),quantities = list(cex = 2),legend=T,cex.lab=1, main="CIMP annotations") -> p
	for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
  	o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
   	if(!is.null(o)){
     	  if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
       		o$children[[paste0("tag.quantity.",i)]]$label <- " "
       		p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
     	  }
	 }
	}
	p -> p3	
	Mutational_status <- list(Cell_bank= colnames(isoforms_cell_bank$TPM),KRAS=colnames(anno_sub)[which(anno_sub["KRAS/BRAF status",]=="KRAS")],BRAF=colnames(anno_sub)[which(anno_sub["KRAS/BRAF status",]=="BRAF")],WT=colnames(anno_sub)[which(anno_sub["KRAS/BRAF status",]=="WT")])
	fit4 <-euler(Mutational_status,shape="ellipse")
	plot(fit4,fills = list(fill = c("seashell2", "indianred4","orange3","grey"), alpha = 0.5),quantities = list(cex = 2),legend=T,cex.lab=1, main="KRAS/BRAF\n mutation annotation") -> p
	for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
  	o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
   	if(!is.null(o)){
     	  if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
       		o$children[[paste0("tag.quantity.",i)]]$label <- " "
       		p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
     	  }
	 }
	}
	p -> p4	
	CMS_class <- list(Cell_bank= colnames(isoforms_cell_bank$TPM),CMS1=colnames(anno_sub)[which(anno_sub["CMS",]=="CMS1")],CMS2=colnames(anno_sub)[which(anno_sub["CMS",]=="CMS2")],CMS3=colnames(anno_sub)[which(anno_sub["CMS",]=="CMS3")],CMS4=colnames(anno_sub)[which(anno_sub["CMS",]=="CMS4")])
	set.seed(666)
	fit4 <-euler(CMS_class,shape="ellipse")
	plot(fit4,fills = list(fill = c("seashell2", "darkorange","steelblue4","orchid4","darkgrey"), alpha = 0.5),quantities = list(cex = 2),legend=T,cex.lab=1, main="CMS class") -> p
	for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
  	o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
   	if(!is.null(o)){
     	  if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
       		o$children[[paste0("tag.quantity.",i)]]$label <- " "
       		p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
     	  }
	 }
	}
	p -> p5	
	CRIS_class <- list(Cell_bank= colnames(isoforms_cell_bank$TPM),CRISA=names(cla)[which(cla=="onlyCRISA")],CRISB=names(cla)[which(cla=="onlyCRISB")],CRISC=names(cla)[which(cla=="onlyCRISC")],CRISD=names(cla)[which(cla=="onlyCRISD")],CRISE=names(cla)[which(cla=="onlyCRISE")],multiCRIS=names(cla)[which(cla%in%c("onlyCRISA","onlyCRISB","onlyCRISC","onlyCRISD","onlyCRISE","NA")==FALSE)])
	set.seed(666)
	fit4 <-euler(CRIS_class,shape="ellipse")
	plot(fit4,fills = list(fill = c("seashell2","red","darkgreen","lightblue","orange","blue","grey"), alpha = 0.5),quantities = list(cex = 2),legend=T,cex.lab=1, main="CRIS class") -> p
	for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
  	o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
   	if(!is.null(o)){
     	  if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
       		o$children[[paste0("tag.quantity.",i)]]$label <- " "
       		p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
     	  }
	 }
	}
	p -> p6	
	plot(plot_grid(p6,p2, p3,p4,p5,p1, labels = c('A','B','C','D','E','F'), label_size = 25,ncol=2,nrow=3))

#	ggvenn(total_anno, fill_color = c("#0073C2FF", "#EFC000FF","#CD534CFF"),stroke_size = 0.5, set_name_size = 4)

	dev.off()
	
	tpm_log <- apply(mat,1,function(x){log2(x+0.1)})
	ltr <- t(apply(tpm_log,2,function(x){x-mean(x)}))
	rownames(ltr) <- rownames(mat)
	v <- apply(mat,1,var)
	ltr_ord <- ltr[order(v,decreasing=T),]
	genes_to_filter_out <- apply(mat,1,function(x){if(quantile(x,0.8)>10){0}else{1}})
        genes_to_filter_out_names <- rownames(mat)[which(genes_to_filter_out==1)]
	ltr_ord_filt <- ltr_ord[setdiff(rownames(ltr_ord),genes_to_filter_out_names),]
        ltr_ord_filt <- ltr_ord_filt[,setdiff(colnames(ltr_ord_filt),c("COGA122.STAR","COLO94H2.STAR","COLO3202.STAR"))]
	col2 <- rep("black",length(colnames(ltr_ord_filt)))
	names(col2) <- colnames(ltr_ord_filt)
	col2[s_r[,1]] <- "orange"
        col2[s_r[,2]] <- "darkred"
        pdf("out/cell_similarity_transcripts.pdf",width=15)
	draw(Heatmap(ltr_ord_filt,cluster_rows=F,show_row_names=F,col=colorRamp2(c(-2,0,2),c("blue","white","red")),column_names_gp=gpar(col=col2,fontsize=5),column_names_side = "top"))
	draw(Heatmap(ltr_ord_filt[,colnames(exprs(isoforms_cell_bank$TPM))],cluster_rows=F,show_row_names=F,col=colorRamp2(c(-2,0,2),c("blue","white","red")),column_names_gp=gpar(col=col2[colnames(exprs(isoforms_cell_bank$TPM))],fontsize=5),column_names_side = "top"))
	dev.off()
	save(isoforms_cell_bank,file="out/annotated_isoforms_cell_bank.rda")
