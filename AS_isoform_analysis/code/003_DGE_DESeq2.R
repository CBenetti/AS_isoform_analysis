###library	
	library(Biobase)
	library(edgeR)
	library(ggplot2)
	library(cowplot)
	library(gridExtra)
	library(grid)
	library(lattice)
	library(xlsx)
	library(ggvenn)
	library(DESeq2)
###loading
	load("out/annotated_genes_counts_cell_bank.rda")
	CRIS_signature <- read.xlsx("data/CRIS_genes.xlsx",sheetIndex=1)
	CRIS_signature$CRIS.Class <- paste("only",CRIS_signature$CRIS.Class,sep="")	
	CRIS_signature$CRIS.Class <- gsub("-","",CRIS_signature$CRIS.Class)
	MSI_signature <- read.xlsx("data/MSI_signature.xlsx",sheetIndex=1)
###functions	
	custom_DGE <- function(x=counts_matrix,group=group){
		rownames(x)<- fData(annotated_genes_counts_cell_bank)$Geneid
		data <- DESeqDataSetFromMatrix(countData = x,colData=data.frame(condition=group),design= ~ 0+ condition)
		fit <- DESeq(data)
		if(length(levels(group))>3){
			results <- list()
			for(i in 1:length(levels(group))){
				contrast <- rep(-(1/(length(levels(group))-1)),length(levels(group)))
				contrast[i] <- 1
				res <- results(fit, contrast=contrast)
				res <- res[which(res$padj!="NA"),c("log2FoldChange","padj")]
                                colnames(res) <- c("logFC","FDR")
				results[[i]] <- as.data.frame(res)

			}
			names(results) <- levels(group)
		}else{
			if(identical(levels(group),c("BRAF","KRAS","WT"))==TRUE){
				results <- list()
				res <- results(fit,contrast=c(1,0,-1))
                                res <- res[which(res$padj!="NA"),c("log2FoldChange","padj")]
                                colnames(res) <- c("logFC","FDR")
				results[["BRAF"]] <- as.data.frame(res)
				res <- results(fit,contrast=c(0,1,-1))
                                res <- res[which(res$padj!="NA"),c("log2FoldChange","padj")]
                                colnames(res) <- c("logFC","FDR")
				results[["KRAS"]] <- as.data.frame(res)
				
			}else{
				#when in pairwise comparison, contrast is assigned to 1 for first level in aplhabetical order
				res <- results(fit,contrast=c(1,-1))
				res	<- as.data.frame(res[which(res$padj!="NA"),c("log2FoldChange","padj")])
				colnames(res) <- c("logFC","FDR")
				results <- as.data.frame(res)	
			}	
		}
		return(results)
	}
###code
	DGE <- list() 
	pdf("out/dispersions_and_volcano_plots_DESeq2.pdf")
##CRIS DGE
	group <- factor(pData(annotated_genes_counts_cell_bank)$CRIS_class[grep("only",pData(annotated_genes_counts_cell_bank)$CRIS_class)])
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,grep("only",pData(annotated_genes_counts_cell_bank)$CRIS_class)]
	DGE[["CRIS_adjusted"]] <- custom_DGE(counts_matrix,group)
	p_dist_FC_cris <- lapply(DGE[["CRIS_adjusted"]],function(dat){
		ggplot(dat, aes(x=logFC)) + 
    		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
    		geom_density(alpha=.2, fill="#FF6666") 	
	})
	p_dist_FDR_cris <- lapply(DGE[["CRIS_adjusted"]],function(dat){
		ggplot(dat, aes(x=FDR)) + 
    		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
    		geom_density(alpha=.2, fill="#FF6666") 	
	})
	p_volcano_cris <- lapply(DGE[["CRIS_adjusted"]],function(dat){ggplot(dat,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})

	grid.arrange(grobs=c(p_dist_FC_cris,p_dist_FDR_cris,p_volcano_cris),layout_matrix=cbind(c(1:5),c(6:10),c(11:15)),top="CRIS A-B")
	inters_sign <- lapply(levels(group),function(x){dat <- DGE[["CRIS_adjusted"]][[x]];ggvenn(list(DGE_results=rownames(dat)[which(dat$logFC>1.5&dat$FDR<0.01)],
	published_signature=CRIS_signature$Gene.ID[which(CRIS_signature$CRIS.Class==x)]),fill_color = c("#0073C2FF", "#EFC000FF"),stroke_size = 0.5, set_name_size = 
4)+ ggtitle(paste(x,"signature",sep=" "))})

##CMS DGE 
	group <- factor(pData(annotated_genes_counts_cell_bank)$CMS_class[which(pData(annotated_genes_counts_cell_bank)$CMS_class!="NA")])
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$CMS_class!="NA")]
	DGE[["CMS_adjusted"]] <- custom_DGE(counts_matrix,group)
	p_dist_FC_cms <- lapply(DGE[["CMS_adjusted"]],function(dat){
                ggplot(dat, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5, 
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_cms <- lapply(DGE[["CMS_adjusted"]],function(dat){
                ggplot(dat, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_cms <- lapply(DGE[["CMS_adjusted"]],function(dat){ggplot(dat,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})
	grid.arrange(grobs=c(p_dist_FC_cms,p_dist_FDR_cms,p_volcano_cms),layout_matrix=cbind(c(1:4),c(5:8),c(9:12)),top="CMS 1-4")

## cetuximab DGE

	group <- factor(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity[which(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
	counts_matrix <-exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity%in%c(0,1)==TRUE)]
	DGE[["Cetuximab_adjusted"]] <- custom_DGE(counts_matrix,group)
	ggplot(DGE[["Cetuximab_adjusted"]],aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_cetuximab
	p_FC_cetuximab <- ggplot(DGE[["Cetuximab_adjusted"]], aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	p_FDR_cetuximab <- ggplot(DGE[["Cetuximab_adjusted"]], aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_cetuximab,p_FDR_cetuximab,volcano_cetuximab),ncol=3,top="Cetuximab sensitivity")


## mutations DGE
	pData(annotated_genes_counts_cell_bank)$mutations[which(pData(annotated_genes_counts_cell_bank)$mutations=="NRAS")] <- NA
	group <- factor(pData(annotated_genes_counts_cell_bank)$mutations[which(pData(annotated_genes_counts_cell_bank)$mutations!="NA")])
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$mutations!="NA")]
	DGE[["mutations_adjusted"]] <- custom_DGE(counts_matrix,group)
	p_dist_FC_mut <- lapply(DGE[["mutations_adjusted"]],function(dat){
                ggplot(dat, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_mut <- lapply(DGE[["mutations_adjusted"]],function(dat){
                ggplot(dat, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_mut <- lapply(DGE[["mutations_adjusted"]],function(dat){ggplot(dat,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})
	grid.arrange(grobs=c(p_dist_FC_mut,p_dist_FDR_mut,p_volcano_mut),layout_matrix=cbind(c(1:2),c(3:4),c(5:6)),top="BRAF-WT and KRAS-WT mutations")

## microsatellite DGE

	group <- factor(pData(annotated_genes_counts_cell_bank)$MS_status[which(pData(annotated_genes_counts_cell_bank)$MS_status!="NA")])
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$MS_status!="NA")]
	DGE[["MS_status_adjusted"]] <- custom_DGE(counts_matrix,group)
	
	ggplot(DGE[["MS_status_adjusted"]],aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_MS
	p_FC_MS <- ggplot(DGE[["MS_status_adjusted"]], aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        p_FDR_MS <- ggplot(DGE[["MS_status_adjusted"]], aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_MS,p_FDR_MS,volcano_MS),ncol=3,top="MSS-MSI")
	
	inters_sign$MSI <- 
ggvenn(list(DGE_results=rownames(DGE[["MS_status_adjusted"]])[which(DGE[["MS_status_adjusted"]]$logFC>1.5 & DGE[["MS_status_adjusted"]]$FDR<0.01)],
        published_signature=MSI_signature$gene),fill_color = c("#0073C2FF", "#EFC000FF"),stroke_size = 0.5, set_name_size = 4)+ ggtitle("MSI signature")
		

## CIMP DGE

	sel <- pData(annotated_genes_counts_cell_bank)$CIMP_status[which(pData(annotated_genes_counts_cell_bank)$CIMP_status!="0")]
	group <- factor(ifelse(sel=="CIMP-H","CIMP-H","rest"))
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$CIMP_status!="0")]
	DGE[["CIMP_adjusted"]] <- custom_DGE(counts_matrix,group)
	ggplot(DGE[["CIMP_adjusted"]],aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_CIMP
	p_FC_CIMP <- ggplot(DGE[["CIMP_adjusted"]], aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        p_FDR_CIMP <- ggplot(DGE[["CIMP_adjusted"]], aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_CIMP,p_FDR_CIMP,volcano_CIMP),ncol=3,top="CIMP-H - CIMP-L")

###out
	
	dev.off()
	pdf("out/intersections_DGE_results_DESeq2.pdf")
	grid.arrange(grobs=inters_sign,nrow=3,top="Intersections with established signatures")
	dev.off()
	save(DGE, file="out/DGE_results_DESeq2.rda")
