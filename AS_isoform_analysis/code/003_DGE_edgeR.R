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
###loading
	load("out/annotated_genes_counts_cell_bank.rda")
	CRIS_signature <- read.xlsx("data/CRIS_genes.xlsx",sheetIndex=1)
	CRIS_signature$CRIS.Class <- paste("only",CRIS_signature$CRIS.Class,sep="")	
	CRIS_signature$CRIS.Class <- gsub("-","",CRIS_signature$CRIS.Class)
	MSI_signature <- read.xlsx("data/MSI_signature.xlsx",sheetIndex=1)
###functions	
	custom_DGE <- function(x=counts_matrix,group=group){
		rownames(x)<- fData(annotated_genes_counts_cell_bank)$Geneid
		data <- DGEList(counts=x,group=group)
		keep <- filterByExpr(data)
		filt_data <- data[keep,,keep.lib.sizes=FALSE]
		filt_data <- normLibSizes(filt_data)
		design <- model.matrix(~0+group)
		filt_data <- estimateDisp(filt_data,design)
		if(length(levels(group))>3){
			results <- list()
			for(i in 1:length(levels(group))){
				contrast <- rep(-(1/(length(levels(group))-1)),length(levels(group)))
				contrast[i] <- 1
				fit <- glmQLFit(filt_data, design)
				results[[i]] <- glmQLFTest(fit, contrast=contrast)
			}
			names(results) <- levels(group)
		}else{
			if(identical(levels(group),c("BRAF","KRAS","WT"))==TRUE){
				fit <- glmQLFit(filt_data, design)
				results <- list()
				results[["BRAF"]] <- glmQLFTest(fit, contrast=c(1,0,-1))
				results[["KRAS"]] <- glmQLFTest(fit, contrast=c(0,1,-1))
				
			}else{
				#remember, when in pairwise comparison, contrast is assigned to -1 for first level in aplhabetical order
				fit <- glmQLFit(filt_data, design)
				results <- glmQLFTest(fit,contrast=c(1,-1))		
			}	
		}
		return(results)
	}
###code
	DGE <- list() 
	pdf("out/dispersions_and_volcano_plots_EdgeR.pdf")
##CRIS DGE
	group <- factor(pData(annotated_genes_counts_cell_bank)$CRIS_class[grep("only",pData(annotated_genes_counts_cell_bank)$CRIS_class)])
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,grep("only",pData(annotated_genes_counts_cell_bank)$CRIS_class)]
	CRIS <- custom_DGE(counts_matrix,group)
	DGE[["CRIS_adjusted"]] <- lapply(CRIS, function(x){topTags(x,n=dim(x$table)[1])})
	p_dist_FC_cris <- lapply(DGE[["CRIS_adjusted"]],function(dat){
		ggplot(dat$table, aes(x=logFC)) + 
    		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
    		geom_density(alpha=.2, fill="#FF6666") 	
	})
	p_dist_FDR_cris <- lapply(DGE[["CRIS_adjusted"]],function(dat){
		ggplot(dat$table, aes(x=FDR)) + 
    		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
    		geom_density(alpha=.2, fill="#FF6666") 	
	})
	p_volcano_cris <- lapply(DGE[["CRIS_adjusted"]],function(dat){ggplot(dat$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})

	grid.arrange(grobs=c(p_dist_FC_cris,p_dist_FDR_cris,p_volcano_cris),layout_matrix=cbind(c(1:5),c(6:10),c(11:15)),top="CRIS A-E")
	inters_sign <- lapply(levels(group),function(x){dat <- DGE[["CRIS_adjusted"]][[x]]$table;ggvenn(list(DGE_results=rownames(dat)[which(dat$logFC>2 & dat$FDR<0.05)],
	published_signature=CRIS_signature$Gene.ID[which(CRIS_signature$CRIS.Class==x)]),fill_color = c("#0073C2FF", "#EFC000FF"),stroke_size = 0.5, set_name_size = 
4)+ ggtitle(paste(x,"signature",sep=" "))})

##CMS DGE 
	group <- factor(pData(annotated_genes_counts_cell_bank)$CMS_class[which(pData(annotated_genes_counts_cell_bank)$CMS_class!="NA")])
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$CMS_class!="NA")]
	CMS <- custom_DGE(counts_matrix,group)
	DGE[["CMS_adjusted"]] <- lapply(CMS, function(x){topTags(x,n=dim(x$table)[1])})
	p_dist_FC_cms <- lapply(DGE[["CMS_adjusted"]],function(dat){
                ggplot(dat$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5, 
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_cms <- lapply(DGE[["CMS_adjusted"]],function(dat){
                ggplot(dat$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_cms <- lapply(DGE[["CMS_adjusted"]],function(dat){ggplot(dat$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})
	grid.arrange(grobs=c(p_dist_FC_cms,p_dist_FDR_cms,p_volcano_cms),layout_matrix=cbind(c(1:4),c(5:8),c(9:12)),top="CMS 1-4")

## cetuximab DGE

	group <- factor(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity[which(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
	counts_matrix <-exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$cetuximab_sensitivity%in%c(0,1)==TRUE)]
	Cetuximab <- custom_DGE(counts_matrix,group)
	DGE[["Cetuximab_adjusted"]] <- topTags(Cetuximab,n=dim(Cetuximab$table)[1])
	ggplot(DGE[["Cetuximab_adjusted"]]$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_cetuximab
	p_FC_cetuximab <- ggplot(DGE[["Cetuximab_adjusted"]]$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	p_FDR_cetuximab <- ggplot(DGE[["Cetuximab_adjusted"]]$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_cetuximab,p_FDR_cetuximab,volcano_cetuximab),ncol=3,top="Cetuximab sensitivity")


## mutations DGE
	pData(annotated_genes_counts_cell_bank)$mutations[which(pData(annotated_genes_counts_cell_bank)$mutations=="NRAS")] <- NA
	group <- factor(pData(annotated_genes_counts_cell_bank)$mutations[which(pData(annotated_genes_counts_cell_bank)$mutations!="NA")])
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$mutations!="NA")]
	mutations <- custom_DGE(counts_matrix,group)
	DGE[["mutations_adjusted"]] <- lapply(mutations, function(x){topTags(x,n=dim(x$table)[1])})
	p_dist_FC_mut <- lapply(DGE[["mutations_adjusted"]],function(dat){
                ggplot(dat$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_mut <- lapply(DGE[["mutations_adjusted"]],function(dat){
                ggplot(dat$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_mut <- lapply(DGE[["mutations_adjusted"]],function(dat){ggplot(dat$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})
	grid.arrange(grobs=c(p_dist_FC_mut,p_dist_FDR_mut,p_volcano_mut),layout_matrix=cbind(c(1:2),c(3:4),c(5:6)),top="BRAF-WT and KRAS-WT mutations")

## microsatellite DGE

	group <- factor(pData(annotated_genes_counts_cell_bank)$MS_status[which(pData(annotated_genes_counts_cell_bank)$MS_status!="NA")])
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$MS_status!="NA")]
	MS_status <- custom_DGE(counts_matrix,group)
	DGE[["MS_status_adjusted"]] <- topTags(MS_status,n=dim(MS_status$table)[1])
	
	ggplot(DGE[["MS_status_adjusted"]]$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_MS
	p_FC_MS <- ggplot(DGE[["MS_status_adjusted"]]$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        p_FDR_MS <- ggplot(DGE[["MS_status_adjusted"]]$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_MS,p_FDR_MS,volcano_MS),ncol=3,top="MSS-MSI")
	
	inters_sign$MSI <- 
ggvenn(list(DGE_results=rownames(DGE[["MS_status_adjusted"]]$table)[which(DGE[["MS_status_adjusted"]]$table$logFC>1.5 & DGE[["MS_status_adjusted"]]$table$FDR<0.01)],
        published_signature=MSI_signature$gene),fill_color = c("#0073C2FF", "#EFC000FF"),stroke_size = 0.5, set_name_size = 4)+ ggtitle("MSI signature")
		

## CIMP DGE

	sel <- pData(annotated_genes_counts_cell_bank)$CIMP_status[which(pData(annotated_genes_counts_cell_bank)$CIMP_status!="0")]
	group <- factor(ifelse(sel=="CIMP-H","CIMP-H","rest"))
	counts_matrix <- exprs(annotated_genes_counts_cell_bank)[,which(pData(annotated_genes_counts_cell_bank)$CIMP_status!="0")]
	CIMP <- custom_DGE(counts_matrix,group)
	DGE[["CIMP_adjusted"]] <- topTags(CIMP,n=dim(CIMP$table)[1])
	ggplot(DGE[["CIMP_adjusted"]]$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_CIMP
	p_FC_CIMP <- ggplot(DGE[["CIMP_adjusted"]]$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        p_FDR_CIMP <- ggplot(DGE[["CIMP_adjusted"]]$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_CIMP,p_FDR_CIMP,volcano_CIMP),ncol=3,top="CIMP-H - CIMP-L")

###out
	
	dev.off()
	pdf("out/intersections_DGE_results.pdf")
	grid.arrange(grobs=inters_sign,nrow=3,top="Intersections with established signatures")
	dev.off()
	save(DGE, file="out/DGE_results.rda")
