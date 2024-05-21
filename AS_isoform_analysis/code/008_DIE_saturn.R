###library
        library(Biobase)
	library(stageR)
	library(DEXSeq)
        library(satuRn)
        library(ggplot2)
        library(cowplot)
        library(gridExtra)
        library(grid)
        library(lattice)
        library(xlsx)
        library(ggvenn)
	library(edgeR)
	library(SummarizedExperiment)
###function
	isoform_DE <- function(x,group,txinfo){
		keep  <- filterByExpr(x,group=group,min.count=10,min.total.count=10)  
		txinfo <- txinfo[names(keep)[which(keep==TRUE)],]
		txinfo <- subset(txinfo, duplicated(gene_id) | duplicated(gene_id , fromLast=TRUE))
		filt <- x[rownames(txinfo),]
		sumExp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=filt),colData=group,rowData=txinfo)
		metadata(sumExp)$formula <- ~ 0 + group
		sumExp <- fitDTU(object=sumExp,formula=~0+group)
		if(length(levels(group))>2){
                	if(identical(levels(group),c("BRAF","KRAS","WT"))==FALSE){
				contrast <- matrix(-1/(length(levels(group))-1),ncol=length(levels(group)),nrow=length(levels(group)))
				colnames(contrast) <- rownames(contrast) <- levels(group)
				diag(contrast) <- 1
			}else{
				contrast <- cbind(c(1,0,-1),c(0,1,-1))
				colnames(contrast) <- c("BRAF","KRAS")
				rownames(contrast) <- levels(group)
			}
			sumExp <- satuRn::testDTU(object = sumExp,contrasts=contrast,sort = FALSE,diagplot1=FALSE,diagplot2=FALSE)        
			# compute gene level q-values
			geneID <- factor(rowData(sumExp)$gene_id)
			geneSplit <- split(seq(along = geneID), geneID)
			tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
			colnames(tx2gene) <- c("transcript", "gene")
			pGandtheta <- lapply(rowData(sumExp)[grep("fitDTUResult",names(rowData(sumExp)))],function(x){
				pvals <- as.data.frame(x)[,grep("empirical_pval",colnames(as.data.frame(x)))]
				pGene <- sapply(geneSplit, function(i) min(pvals[i]))
				pGene[is.na(pGene)] <- 1
				tetha <- unique(sort(pGene))
				pConfirmation <- matrix(matrix(pvals),ncol = 1,dimnames = list(rownames(tx2gene), "transcript"))
				return(list(pGene,tetha,pConfirmation))
			})
			# gene-level significance testing
			adj_pval <- lapply(pGandtheta, function(x){
				q <- DEXSeq:::perGeneQValueExact(x[[1]], x[[2]], geneSplit) 
				qScreen <- rep(NA_real_, length(x[[1]]))
				qScreen <- q[match(x[[1]], x[[2]])]
				qScreen <- pmin(1, qScreen)
				names(qScreen) <- names(geneSplit)
				stageRObj <- stageR::stageRTx(pScreen = qScreen,pConfirmation = x[[3]], pScreenAdjusted = TRUE,tx2gene = tx2gene)
				stageRObj <- stageR::stageWiseAdjustment(object = stageRObj,method = "dtu",alpha = 0.05,allowNA = TRUE)
				padj <-stageR::getAdjustedPValues(stageRObj,order = TRUE,onlySignificantGenes = FALSE)
				rownames(padj) <- padj$txID
				return(padj)
			})
			nam <- names(adj_pval) <- names(rowData(sumExp))[grep("fitDTUResult",names(rowData(sumExp)))]
			results <- lapply(nam,function(x){
data.frame(logFC=as.data.frame(rowData(sumExp)[[x]])[adj_pval[[x]]$txID,grep("estimates",colnames(as.data.frame(rowData(sumExp)[[x]])))],FDR=adj_pval[[x]]$transcript,Geneid=adj_pval[[x]]$geneID,rownam=adj_pval[[x]]$txID)
			})
			names(results) <- colnames(contrast)                
                }else{
			#when in pairwise comparison, contrast is assigned to 1 for first level in aplhabetical order
			sumExp <- satuRn::testDTU(object = sumExp,contrasts=data.frame(Contrast1 = c(1,-1)),diagplot1=TRUE,diagplot2=TRUE,sort = FALSE)        
			x <- as.data.frame(rowData(sumExp)[[grep("fitDTUResult",names(rowData(sumExp)))]])
			# transcript level p-values from satuRn
			pvals <- x[,grep("empirical_pval",colnames(x))]
			# compute gene level q-values
			geneID <- factor(rowData(sumExp)$gene_id)
			geneSplit <- split(seq(along = geneID), geneID)
			pGene <- sapply(geneSplit, function(i) min(pvals[i]))
			pGene[is.na(pGene)] <- 1
			theta <- unique(sort(pGene))
			# gene-level significance testing
			q <- DEXSeq:::perGeneQValueExact(pGene, theta, geneSplit) 
			qScreen <- rep(NA_real_, length(pGene))
			qScreen <- q[match(pGene, theta)]
			qScreen <- pmin(1, qScreen)
			names(qScreen) <- names(geneSplit)
			# prepare stageR input
			tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
			colnames(tx2gene) <- c("transcript", "gene")
			pConfirmation <- matrix(matrix(pvals),ncol = 1,dimnames = list(rownames(tx2gene), "transcript"))
			# create a stageRTx object
			stageRObj <- stageR::stageRTx(pScreen = qScreen, pConfirmation = pConfirmation, pScreenAdjusted = TRUE,tx2gene = tx2gene)
			# perform the two-stage testing procedure
			stageRObj <- stageR::stageWiseAdjustment(object = stageRObj,method = "dtu",alpha = 0.05,allowNA = TRUE)
			# retrieves the adjusted p-values from the stageRTx object
			padj <- stageR::getAdjustedPValues(stageRObj, order = TRUE, onlySignificantGenes = FALSE)
			rownames(padj) <- padj$txID
			results <- data.frame(logFC=x$estimates,FDR=padj[rownames(x),"transcript"],Geneid=txInfo[rownames(x),"gene_id"],rownam=rownames(x))
			rownames(results) <- results$rownam
		}
		return(results)
	}

###loading
	load("out/annotated_isoforms_cell_bank_counts.rda")
###code

	txInfo <- as.data.frame(matrix(data=NA,ncol=2,nrow=length(fData(isoforms_cell_bank$counts_75)$Geneid)))
	colnames(txInfo) = c("isoform_id","gene_id")
	rownames(txInfo) <- txInfo$isoform_id <- rownames(exprs(isoforms_cell_bank$counts_75))
	txInfo$gene_id <- fData(isoforms_cell_bank$counts_75)$Geneid
	txInfo <- subset(txInfo, duplicated(gene_id) | duplicated(gene_id , fromLast=TRUE))
	isoform_counts <- exprs(isoforms_cell_bank$counts_75)[rownames(txInfo),]
	DIE <- list()
	pdf("out/diaplots_satuRn.pdf")
##CRIS DIE
        group <- factor(pData(isoforms_cell_bank$counts_75)$CRIS_class[grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)])
        counts_matrix <- isoform_counts[,grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)]
        DIE[["CRIS_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
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
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_cris <- lapply(DIE[["CRIS_adjusted"]],function(dat){ggplot(dat,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})


##CMS DIE
        group <- factor(pData(isoforms_cell_bank$counts_75)$CMS_class[which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")])
        counts_matrix <- isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")]
        DIE[["CMS_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
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
        
## cetuximab DIE
        
        group <- factor(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity[which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
        counts_matrix <-isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)]
        DIE[["Cetuximab_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
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
        counts_matrix <- isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")]
        DIE[["mutations_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
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
                
## microsatellite DIE
                
        group <- factor(pData(isoforms_cell_bank$counts_75)$MS_status[which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")])
        counts_matrix <- isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")]
        DIE[["MS_status_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)

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
        counts_matrix <- isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")]
        DIE[["CIMP_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
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
                
	dev.off()
###out
        pdf("out/dispersions_and_volcano_plots_satuRn.pdf")
	grid.arrange(grobs=c(p_dist_FC_cris,p_dist_FDR_cris,p_volcano_cris),layout_matrix=cbind(c(1:5),c(6:10),c(11:15)),top="CRIS A-B")
        grid.arrange(grobs=c(p_dist_FC_cms,p_dist_FDR_cms,p_volcano_cms),layout_matrix=cbind(c(1:4),c(5:8),c(9:12)),top="CMS 1-4")
	grid.arrange(grobs=list(p_FC_cetuximab,p_FDR_cetuximab,volcano_cetuximab),ncol=3,top="Cetuximab sensitivity")
	grid.arrange(grobs=c(p_dist_FC_mut,p_dist_FDR_mut,p_volcano_mut),layout_matrix=cbind(c(1:2),c(3:4),c(5:6)),top="BRAF-WT and KRAS-WT mutations")
        grid.arrange(grobs=list(p_FC_MS,p_FDR_MS,volcano_MS),ncol=3,top="MSS-MSI")
	grid.arrange(grobs=list(p_FC_CIMP,p_FDR_CIMP,volcano_CIMP),ncol=3,top="CIMP-H - CIMP-L")
	dev.off()
        save(DIE, file="out/DIE_results_satuRn.rda")
