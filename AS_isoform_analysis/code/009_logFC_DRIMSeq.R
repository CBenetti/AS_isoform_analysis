###library
	library(matrixStats)
	library(Biobase)
	library(DRIMSeq)
	library(ggplot2)
	library(ggvenn)
	library(cowplot)
        library(gridExtra)
        library(grid)
        library(lattice)
	library(AnnotationHub)
	library(ensembldb)
###function

ProportionSD <- function(res){
		cts <- isoform_counts[res$rownam,]
		gid <- txInfo[res$rownam,"gene_id"]
		gene_cts <- rowsum(cts, gid)
		total_cts <- gene_cts[match(gid, rownames(gene_cts)),]
		props <- cts/total_cts
		return(sqrt(rowVars(props)))
	}

###loading
	load("out/DIE_results_DRIMSeq.rda")
	load("out/annotated_isoforms_cell_bank_counts.rda")
###code
        txInfo <- as.data.frame(matrix(data=NA,ncol=2,nrow=length(fData(isoforms_cell_bank$counts_75)$Geneid)))
        colnames(txInfo) = c("isoform_id","gene_id")
        rownames(txInfo) <- txInfo$isoform_id <- rownames(exprs(isoforms_cell_bank$counts_75))
        txInfo$gene_id <- fData(isoforms_cell_bank$counts_75)$Geneid
        txInfo <- subset(txInfo, duplicated(gene_id) | duplicated(gene_id , fromLast=TRUE)) 
	isoform_counts <- exprs(isoforms_cell_bank$counts_75)[rownames(txInfo),] 
	isoform_tpm <- exprs(isoforms_cell_bank$TPM)[rownames(txInfo),] 
	isoform_ratio <- lapply(unique(txInfo$gene_id),function(x){prop.table(isoform_tpm[which(txInfo$gene_id==x),], margin =2)})
	names(isoform_ratio) <- unique(txInfo$gene_id)
	Geneid <- unlist(lapply(unique(txInfo$gene_id),function(x){rep(x,dim(isoform_ratio[[x]])[1])}))	
	isoform_ratio <- do.call(rbind,isoform_ratio)
	isoform_ratio[is.na(isoform_ratio)] <- 0

	
##CRIS
	group <- factor(pData(isoforms_cell_bank$counts_75)$CRIS_class[grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)])
	cells <- grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)
	for(x in names(DIE[["CRIS_adjusted"]])){
		DIE[["CRIS_adjusted"]][[x]]$logFC <- apply(isoform_ratio[DIE[["CRIS_adjusted"]][[x]]$rownam,cells],1,function(a){
			log2(mean(a[which(group==x)]+0.01)/mean(a[which(group%in%setdiff(levels(group),x)==TRUE)]+0.001))
		})
		DIE[["CRIS_adjusted"]][[x]]$proportionSD <- ProportionSD(DIE[["CRIS_adjusted"]][[x]])
	}

	p_dist_FC_cris <- lapply(DIE[["CRIS_adjusted"]],function(dat){
                ggplot(dat, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_cris <- lapply(DIE[["CRIS_adjusted"]],function(dat){
                ggplot(dat, aes(x=FDR)) +
		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_cris <- lapply(DIE[["CRIS_adjusted"]],function(dat){ggplot(dat,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})

##CMS	
	group <- factor(pData(isoforms_cell_bank$counts_75)$CMS_class[which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")])
        cells <- which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")
	for(x in names(DIE[["CMS_adjusted"]])){
		DIE[["CMS_adjusted"]][[x]]$logFC <- apply(isoform_ratio[DIE[["CMS_adjusted"]][[x]]$rownam,cells],1,function(a){
			log2(mean(a[which(group==x)]+0.01)/mean(a[which(group%in%setdiff(levels(group),x)==TRUE)]+0.001))
		})
		DIE[["CMS_adjusted"]][[x]]$proportionSD <- ProportionSD(DIE[["CMS_adjusted"]][[x]])
	}

	p_dist_FC_cms <- lapply(DIE[["CMS_adjusted"]],function(dat){
                ggplot(dat, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_cms <- lapply(DIE[["CMS_adjusted"]],function(dat){
                ggplot(dat, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })                      
        p_volcano_cms <- lapply(DIE[["CMS_adjusted"]],function(dat){ggplot(dat,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})
##cetuximab
	group <- factor(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity[which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
        cells <- which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)
	DIE[["Cetuximab_adjusted"]]$logFC <- apply(isoform_ratio[DIE[["Cetuximab_adjusted"]]$rownam,cells],1,function(a){
		log2(mean(a[which(group==1)]+0.01)/mean(a[which(group%in%setdiff(levels(group),1)==TRUE)]+0.001))
	})
	DIE[["Cetuximab_adjusted"]]$proportionSD <- ProportionSD(DIE[["Cetuximab_adjusted"]])

	ggplot(DIE[["Cetuximab_adjusted"]],aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_cetuximab
        p_FC_cetuximab <- ggplot(DIE[["Cetuximab_adjusted"]], aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        p_FDR_cetuximab <- ggplot(DIE[["Cetuximab_adjusted"]], aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
                
## mutations DIE
	pData(isoforms_cell_bank$counts_75)$mutations[which(pData(isoforms_cell_bank$counts_75)$mutations=="NRAS")] <- NA
        group <- factor(pData(isoforms_cell_bank$counts_75)$mutations[which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")])
	cells <- which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")
	DIE[["mutations_adjusted"]][["BRAF"]]$logFC <- apply(isoform_ratio[DIE[["mutations_adjusted"]][["BRAF"]]$rownam,cells],1,function(a){
		log2(mean(a[which(group=="BRAF")]+0.01)/mean(a[which(group=="WT")]+0.001))
	})
	DIE[["mutations_adjusted"]][["KRAS"]]$logFC <- apply(isoform_ratio[DIE[["mutations_adjusted"]][["KRAS"]]$rownam,cells],1,function(a){
		log2(mean(a[which(group=="KRAS")]+0.01)/mean(a[which(group=="WT")]+0.001))
	})
		DIE[["mutations_adjusted"]][["BRAF"]]$proportionSD <- DIE[["mutations_adjusted"]][["KRAS"]]$proportionSD <- ProportionSD(DIE[["mutations_adjusted"]][["BRAF"]])

	p_dist_FC_mut <- lapply(DIE[["mutations_adjusted"]],function(dat){
                ggplot(dat, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_mut <- lapply(DIE[["mutations_adjusted"]],function(dat){
                ggplot(dat, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_mut <- lapply(DIE[["mutations_adjusted"]],function(dat){ggplot(dat,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})

## MS status
	group <- factor(pData(isoforms_cell_bank$counts_75)$MS_status[which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")])
	cells <- which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")
	DIE[["MS_status_adjusted"]]$logFC <- apply(isoform_ratio[DIE[["MS_status_adjusted"]]$rownam,cells],1,function(a){
		log2(mean(a[which(group=="MSI")]+0.01)/mean(a[which(group=="MSS")]+0.001))
	})
	DIE[["MS_status_adjusted"]]$proportionSD <- ProportionSD(DIE[["MS_status_adjusted"]])

        ggplot(DIE[["MS_status_adjusted"]],aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_MS
        p_FC_MS <- ggplot(DIE[["MS_status_adjusted"]], aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
               geom_density(alpha=.2, fill="#FF6666")
        p_FDR_MS <- ggplot(DIE[["MS_status_adjusted"]], aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")

## CIMP DIE
        sel <- pData(isoforms_cell_bank$counts_75)$CIMP_status[which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")]
        group <- factor(ifelse(sel=="CIMP-H","CIMP-H","rest"))
	cells <- which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")
	DIE[["CIMP_adjusted"]]$logFC <- apply(isoform_ratio[DIE[["CIMP_adjusted"]]$rownam,cells],1,function(a){
		log2(mean(a[which(group=="CIMP-H")]+0.01)/mean(a[which(group=="rest")]+0.001))
	})
	DIE[["CIMP_adjusted"]]$proportionSD <- ProportionSD(DIE[["CIMP_adjusted"]])

        ggplot(DIE[["CIMP_adjusted"]],aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_CIMP
        p_FC_CIMP <- ggplot(DIE[["CIMP_adjusted"]], aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        p_FDR_CIMP <- ggplot(DIE[["CIMP_adjusted"]], aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")

        ah <- AnnotationHub()
        ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 110))
        ahEdb <- ahDb[[1]]
        transcr_data <- transcripts(ahEdb)
	isoforms_cell_bank$isoform_proportions <- ExpressionSet(assayData=as.matrix(isoform_ratio))	
	
	fData(isoforms_cell_bank$isoform_proportions) <- 
data.frame(Geneid=Geneid,biotype=transcr_data[match(rownames(isoform_ratio),transcr_data$tx_id_version,nomatch=0)]$tx_biotype,is_canonical=transcr_data[match(rownames(isoform_ratio),transcr_data$tx_id_version,nomatch=0)]$tx_is_canonical)
	pData(isoforms_cell_bank$isoform_proportions) <- pData(isoforms_cell_bank$counts_75) 
                
###out

	save(DIE,file="out/DIE_results_DRIMSeq_withL2R.rda")
        pdf("out/dispersions_and_volcano_plots_DRIMSeq.pdf")
        grid.arrange(grobs=c(p_dist_FC_cris,p_dist_FDR_cris,p_volcano_cris),layout_matrix=cbind(c(1:5),c(6:10),c(11:15)),top="CRIS A-B")
        grid.arrange(grobs=c(p_dist_FC_cms,p_dist_FDR_cms,p_volcano_cms),layout_matrix=cbind(c(1:4),c(5:8),c(9:12)),top="CMS 1-4")
        grid.arrange(grobs=list(p_FC_cetuximab,p_FDR_cetuximab,volcano_cetuximab),ncol=3,top="Cetuximab sensitivity")
        grid.arrange(grobs=c(p_dist_FC_mut,p_dist_FDR_mut,p_volcano_mut),layout_matrix=cbind(c(1:2),c(3:4),c(5:6)),top="BRAF-WT and KRAS-WT mutations")
        grid.arrange(grobs=list(p_FC_MS,p_FDR_MS,volcano_MS),ncol=3,top="MSS-MSI")
        grid.arrange(grobs=list(p_FC_CIMP,p_FDR_CIMP,volcano_CIMP),ncol=3,top="CIMP-H - CIMP-L")
        dev.off()
	save(isoforms_cell_bank, file="out/annotated_isoforms_cell_bank_proportions.rda")

