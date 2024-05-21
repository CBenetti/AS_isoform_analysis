###library
        library(Biobase)
	library(stageR)
	library(DRIMSeq)
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
		txinfo <- txinfo[,c("gene_id","isoform_id")]
		colnames(txinfo) <- c("gene_id","feature_id")
		data <- dmDSdata(counts = cbind(txinfo,x), samples = data.frame(sample_id=colnames(x),group=group))
		data <- dmFilter(data, min_samps_gene_expr = dim(x)[2], min_samps_feature_expr = min(table(group)),min_gene_expr = 10, min_feature_expr = 10)
		design <- model.matrix(~ 0 + group, data = samples(data))
		set.seed(123)
		data  <- dmPrecision(data, design = design)
		print(plotPrecision(data))
		data <- dmFit(data, design=design, one_way=TRUE, bb_model = TRUE)
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
			res <- lapply(levels(group)[which(levels(group)!="WT")],function(x){
				dx <-dmTest(data, contrast = contrast[,x])
				return(list(results(dx),results(dx,level="feature")))
			})
			adj_pval <-lapply(res,function(x){
				tx2gene <- x[[2]][,c("feature_id", "gene_id")]
				pConfirmation <- matrix(x[[2]]$pvalue, ncol = 1)
				rownames(pConfirmation) <- x[[2]]$feature_id
				pConfirmation[is.na(pConfirmation)] <- 1
				pScreen <- x[[1]]$pvalue
				pScreen[is.na(pScreen)] <- 1
				names(pScreen) <- x[[1]]$gene_id
				stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,pScreenAdjusted = FALSE, tx2gene = tx2gene)
				stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",alpha = 0.05)
				padj <- getAdjustedPValues(stageRObj, order = TRUE,onlySignificantGenes = FALSE)
				rownames(padj) <- padj$txID
				return(padj)					
			})			
			nam <- names(adj_pval) <- names(res) <- levels(group)[which(levels(group)!="WT")]
			results <- lapply(nam,function(x){
data.frame(logFC=res[[x]][[2]]$lr,FDR=adj_pval[[x]][res[[x]][[2]]$feature_id,"transcript"],Geneid=res[[x]][[2]]$gene_id,rownam=res[[x]][[2]]$feature_id)
			})
			names(results) <- levels(group)[which(levels(group)!="WT")]                
                }else{
			#when in pairwise comparison, contrast is assigned to 1 for first level in aplhabetical order
			data <- dmTest(data, contrast = c(1,-1))
			res_g <- results(data)
			res_t <- results(data,level="feature")
			tx2gene <- res_t[,c("feature_id", "gene_id")]
			pScreen <- res_g$pvalue
			pScreen[is.na(pScreen)] <- 1
			names(pScreen) <- res_g$gene_id
			pConfirmation <- matrix(res_t$pvalue, ncol = 1)
			rownames(pConfirmation) <- res_t$feature_id
			pConfirmation[is.na(pConfirmation)] <- 1
			stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,pScreenAdjusted = FALSE, tx2gene = tx2gene)
                        stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",alpha = 0.05)
                        padj <- getAdjustedPValues(stageRObj, order = TRUE,onlySignificantGenes = FALSE)
                        rownames(padj) <- padj$txID 
			results <- data.frame(logFC=res_t$lr,FDR=padj[res_t$feature_id,"transcript"],Geneid=res_t$gene_id,rownam=res_t$feature_id)
		}
		return(results)
	}

###loading
	load("out/annotated_isoforms_cell_bank_counts.rda")
	load("out/DGE_results.rda")
###code

	txInfo <- as.data.frame(matrix(data=NA,ncol=2,nrow=length(fData(isoforms_cell_bank$counts_75)$Geneid)))
	colnames(txInfo) = c("isoform_id","gene_id")
	rownames(txInfo) <- txInfo$isoform_id <- rownames(exprs(isoforms_cell_bank$counts_75))
	txInfo$gene_id <- fData(isoforms_cell_bank$counts_75)$Geneid
	txInfo <- subset(txInfo, duplicated(gene_id) | duplicated(gene_id , fromLast=TRUE))
	isoform_counts <- exprs(isoforms_cell_bank$counts_75)[rownames(txInfo),]
	DIE <- list()
	pdf("out/precision_DRIMSeq.pdf")
##CRIS DIE
        group <- factor(pData(isoforms_cell_bank$counts_75)$CRIS_class[grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)])
        counts_matrix <- isoform_counts[,grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)]
        DIE[["CRIS_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
pdf("out/DTU_genes.pdf")
        db <- data.frame(classes=rep(levels(group),2),status=c(rep("shared",length(levels(group))),rep("not_shared",length(levels(group)))))
        for(i in levels(group)){
                DGE_all <-  rownames(DGE[[1]][[i]])[which(DGE[[1]][[i]]$table$FDR < 0.05 & abs(DGE[[1]][[i]]$table$logFC)>2)]
                DTU_transcripts <- DIE[[1]][[i]]$rownam[which(DIE[[1]][[i]]$FDR < 0.01)]
                DTU_genes <- unique(fData(isoforms_cell_bank$cov)$Geneid[match(DTU_transcripts,rownames(isoforms_cell_bank$cov))])
		plot(ggvenn(list(DGE_genes = DGE_all,DTU_genes=DTU_genes),fill_color = c("#00BFC4", "#F8766D"),stroke_size = 0.5, set_name_size =4)+ 
ggtitle(paste(i,"results",sep=" ")))
                db[which(db$status=="shared"& db$class==i),"gene_number"] <- length(intersect(DTU_genes,DGE_all))
                db[which(db$status=="not_shared"& db$class==i),"gene_number"] <- length(setdiff(DTU_genes,DGE_all))
        } 
        gsub("not_shared","Not common to significative DGE related genes",db$status) ->db$status
        gsub("shared","Common to significative DGE related genes",db$status) ->db$status
        db$classes <- sapply(strsplit(db$classes,"only"),"[[",2)
        plot(ggplot(data=db, aes(x=classes, y=gene_number, fill=status)) +
  geom_bar(stat="identity")+ggtitle("DTU significative genes"))
        dev.off()
##CMS DIE
        group <- factor(pData(isoforms_cell_bank$counts_75)$CMS_class[which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")])
        counts_matrix <- isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")]
        DIE[["CMS_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
        
## cetuximab DIE
        
        group <- factor(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity[which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
        counts_matrix <-isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)]
        DIE[["Cetuximab_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
        
                
## mutations DIE
        pData(isoforms_cell_bank$counts_75)$mutations[which(pData(isoforms_cell_bank$counts_75)$mutations=="NRAS")] <- NA
        group <- factor(pData(isoforms_cell_bank$counts_75)$mutations[which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")])
        counts_matrix <- isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")]
        DIE[["mutations_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
                
## microsatellite DIE
                
        group <- factor(pData(isoforms_cell_bank$counts_75)$MS_status[which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")])
        counts_matrix <- isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")]
        DIE[["MS_status_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)

          
## CIMP DIE
        
        sel <- pData(isoforms_cell_bank$counts_75)$CIMP_status[which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")]
        group <- factor(ifelse(sel=="CIMP-H","CIMP-H","rest"))
        counts_matrix <- isoform_counts[,which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")]
        DIE[["CIMP_adjusted"]] <- isoform_DE(counts_matrix,group,txInfo)
	dev.off()
###out
        save(DIE, file="out/DIE_results_DRIMSeq.rda")
