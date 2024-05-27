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
	load("out/annotated_isoforms_cell_bank_counts.rda")
	load("out/DGE_results.rda")

###functions	
	custom_DTE <- function(x=counts_matrix,group=group){
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
	DGE_signif <- lapply(DGE,function(x){
		if("table"%in%names(x)==FALSE){
			z=lapply(x,function(z){	
					unlist(lapply(fData(isoforms_cell_bank$counts_75)$Geneid,function(nam){
						if(nam%in%rownames(z$table)==TRUE){
							if(z$table[nam,"FDR"]<0.01 & abs(z$table[nam,"logFC"])>2){
								if(z$table[nam,"logFC"]>2){1}else{-1}
							}else{0}
						}else{0}
					}))
			})
			do.call(cbind,z)
		}else{
			unlist(lapply(fData(isoforms_cell_bank$counts_75)$Geneid,function(nam){
				if(nam%in%rownames(x$table)==TRUE){
					if(x$table[nam,"FDR"]<0.01&abs(x$table[nam,"logFC"])>2){
						if(x$table[nam,"logFC"]>2){1}else{-1}
					}else{0}
				}else{0}
			}))
		}
	})
	DGE_signif_res <- do.call(cbind,DGE_signif)
	colnames(DGE_signif_res) <- unlist(lapply(names(DGE),function(x){if("table"%in%names(DGE[[x]])==TRUE){x}else{names(DGE[[x]])}}))
	rownames(DGE_signif_res) <- rownames(exprs(isoforms_cell_bank$counts_75))
	DTE <- list() 
	pdf("out/dispersions_and_volcano_plots_EdgeR.pdf")
##CRIS DTE
	group <- factor(pData(isoforms_cell_bank$counts_75)$CRIS_class[grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)])
	counts_matrix <- exprs(isoforms_cell_bank$counts_75)[,grep("only",pData(isoforms_cell_bank$counts_75)$CRIS_class)]
	CRIS <- custom_DTE(counts_matrix,group)
	DTE[["CRIS_adjusted"]] <- lapply(CRIS, function(x){topTags(x,n=dim(x$table)[1])})
	p_dist_FC_cris <- lapply(DTE[["CRIS_adjusted"]],function(dat){
		ggplot(dat$table, aes(x=logFC)) + 
    		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
    		geom_density(alpha=.2, fill="#FF6666") 	
	})
	p_dist_FDR_cris <- lapply(DTE[["CRIS_adjusted"]],function(dat){
		ggplot(dat$table, aes(x=FDR)) + 
    		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
    		geom_density(alpha=.2, fill="#FF6666") 	
	})
	p_volcano_cris <- lapply(DTE[["CRIS_adjusted"]],function(dat){ggplot(dat$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})
	for(x in names(DTE[["CRIS_adjusted"]])){DTE[["CRIS_adjusted"]][[x]]$table$DGE_res <- DGE_signif_res[rownames(DTE[["CRIS_adjusted"]][[x]]$table),x]}	
	p_volcano_int_cris <- lapply(names(DTE[["CRIS_adjusted"]]),function(x){ggplot(DTE[["CRIS_adjusted"]][[x]]$table,aes(x=logFC,y=-log(FDR),col=as.character(DGE_res)))+
		geom_point(size=0.5,)+scale_color_manual(values=c("-1"="lightblue","0"="black","1"="orange"),name="DGE\nresults")+ggtitle(x)})
	grid.arrange(grobs=c(p_dist_FC_cris,p_dist_FDR_cris,p_volcano_cris),layout_matrix=cbind(c(1:5),c(6:10),c(11:15)),top="CRIS A-B")
	pdf("out/DGE_genes.pdf")
	db <- data.frame(classes=rep(levels(group),2),status=c(rep("shared",length(levels(group))),rep("not_shared",length(levels(group)))))
	db_t <- db
	for(i in levels(group)){
		DGE_all <-  rownames(DGE[[1]][[i]])[which(DGE[[1]][[i]]$table$FDR < 0.05 & abs(DGE[[1]][[i]]$table$logFC)>2)]
		DTE_transcripts <- rownames(DTE[[1]][[i]])[which(DTE[[1]][[i]]$table$FDR < 0.05 & abs(DTE[[1]][[i]]$table$logFC)>2)]
		DTE_genes <- unique(fData(isoforms_cell_bank$cov)$Geneid[match(DTE_transcripts,rownames(isoforms_cell_bank$cov))])
		db[which(db$status=="shared"& db$class==i),"gene_number"] <- length(intersect(DGE_all,DTE_genes))
		db[which(db$status=="not_shared"& db$class==i),"gene_number"] <- length(setdiff(DGE_all,DTE_genes))		
		db_t[which(db_t$status=="shared"& db$class==i),"gene_number"] <- length(intersect(DTE_genes,DGE_all))
		db_t[which(db_t$status=="not_shared"& db$class==i),"gene_number"] <- length(setdiff(DTE_genes,DGE_all))		


	}
	gsub("not_shared","Not common to significative DTE related genes",db$status) ->db$status 
	gsub("not_shared","Not common to significative DGE genes",db_t$status) ->db_t$status 
	gsub("shared","Common to significative DTE related genes",db$status) ->db$status
	gsub("shared","Common to significative DGE genes",db_t$status) ->db_t$status
	db_t$classes <- db$classes <- sapply(strsplit(db$classes,"only"),"[[",2)
	plot(ggplot(data=db, aes(x=classes, y=gene_number, fill=status)) +
  geom_bar(stat="identity")+ggtitle("DGE significative genes"))
	        plot(ggplot(data=db_t, aes(x=classes, y=gene_number, fill=status)) +
  geom_bar(stat="identity")+ggtitle("DTE significative genes"))
	dev.off()
##CMS DTE 
	group <- factor(pData(isoforms_cell_bank$counts_75)$CMS_class[which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")])
	counts_matrix <- exprs(isoforms_cell_bank$counts_75)[,which(pData(isoforms_cell_bank$counts_75)$CMS_class!="NA")]
	CMS <- custom_DTE(counts_matrix,group)
	DTE[["CMS_adjusted"]] <- lapply(CMS, function(x){topTags(x,n=dim(x$table)[1])})
	p_dist_FC_cms <- lapply(DTE[["CMS_adjusted"]],function(dat){
                ggplot(dat$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5, 
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_cms <- lapply(DTE[["CMS_adjusted"]],function(dat){
                ggplot(dat$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_cms <- lapply(DTE[["CMS_adjusted"]],function(dat){ggplot(dat$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})
	for(x in names(DTE[["CMS_adjusted"]])){DTE[["CMS_adjusted"]][[x]]$table$DGE_res <- DGE_signif_res[rownames(DTE[["CMS_adjusted"]][[x]]$table),x]}
	p_volcano_int_cms <- lapply(names(DTE[["CMS_adjusted"]]),function(x){ggplot(DTE[["CMS_adjusted"]][[x]]$table,aes(x=logFC,y=-log(FDR),col=as.character(DGE_res)))+
		geom_point(size=0.5)+scale_color_manual(values=c("-1"="lightblue","0"="black","1"="red"),name="DGE\nresults")+ggtitle(x)})

	grid.arrange(grobs=c(p_dist_FC_cms,p_dist_FDR_cms,p_volcano_cms),layout_matrix=cbind(c(1:4),c(5:8),c(9:12)),top="CMS 1-4")

## cetuximab DTE

	group <- factor(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity[which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)])
	counts_matrix <-exprs(isoforms_cell_bank$counts_75)[,which(pData(isoforms_cell_bank$counts_75)$cetuximab_sensitivity%in%c(0,1)==TRUE)]
	Cetuximab <- custom_DTE(counts_matrix,group)
	DTE[["Cetuximab_adjusted"]] <- topTags(Cetuximab,n=dim(Cetuximab$table)[1])
	ggplot(DTE[["Cetuximab_adjusted"]]$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_cetuximab
	p_FC_cetuximab <- ggplot(DTE[["Cetuximab_adjusted"]]$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	p_FDR_cetuximab <- ggplot(DTE[["Cetuximab_adjusted"]]$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_cetuximab,p_FDR_cetuximab,volcano_cetuximab),ncol=3,top="Cetuximab sensitivity")
	DTE[["Cetuximab_adjusted"]]$table$DGE_res <- DGE_signif_res[rownames(DTE[["Cetuximab_adjusted"]]$table),"Cetuximab_adjusted"]
ggplot(DTE[["Cetuximab_adjusted"]]$table,aes(x=logFC,y=-log(FDR),col=as.character(DGE_res)))+geom_point(size=0.5)+scale_color_manual(values=c("-1"="lightblue","0"="black","1"="red"),name="DGE\nresults")+ggtitle("Cetuximab\nresponse") -> volcano_int_cetuximab


## mutations DTE
	pData(isoforms_cell_bank$counts_75)$mutations[which(pData(isoforms_cell_bank$counts_75)$mutations=="NRAS")] <- NA
	group <- factor(pData(isoforms_cell_bank$counts_75)$mutations[which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")])
	counts_matrix <- exprs(isoforms_cell_bank$counts_75)[,which(pData(isoforms_cell_bank$counts_75)$mutations!="NA")]
	mutations <- custom_DTE(counts_matrix,group)
	DTE[["mutations_adjusted"]] <- lapply(mutations, function(x){topTags(x,n=dim(x$table)[1])})
	p_dist_FC_mut <- lapply(DTE[["mutations_adjusted"]],function(dat){
                ggplot(dat$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_dist_FDR_mut <- lapply(DTE[["mutations_adjusted"]],function(dat){
                ggplot(dat$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        })
        p_volcano_mut <- lapply(DTE[["mutations_adjusted"]],function(dat){ggplot(dat$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5)})
	for(x in names(DTE[["mutations_adjusted"]])){DTE[["mutations_adjusted"]][[x]]$table$DGE_res <- DGE_signif_res[rownames(DTE[["mutations_adjusted"]][[x]]$table),x]}
	p_volcano_int_mut <- lapply(names(DTE[["mutations_adjusted"]]),function(x){ggplot(DTE[["mutations_adjusted"]][[x]]$table,aes(x=logFC,y=-log(FDR),col=as.character(DGE_res)))+
	geom_point(size=0.5)+scale_color_manual(values=c("-1"="lightblue","0"="black","1"="red"),name="DGE\nresults")+ggtitle(x)})

	grid.arrange(grobs=c(p_dist_FC_mut,p_dist_FDR_mut,p_volcano_mut),layout_matrix=cbind(c(1:2),c(3:4),c(5:6)),top="BRAF-WT and KRAS-WT mutations")

## microsatellite DTE

	group <- factor(pData(isoforms_cell_bank$counts_75)$MS_status[which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")])
	counts_matrix <- exprs(isoforms_cell_bank$counts_75)[,which(pData(isoforms_cell_bank$counts_75)$MS_status!="NA")]
	MS_status <- custom_DTE(counts_matrix,group)
	DTE[["MS_status_adjusted"]] <- topTags(MS_status,n=dim(MS_status$table)[1])
	
	ggplot(DTE[["MS_status_adjusted"]]$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_MS
	p_FC_MS <- ggplot(DTE[["MS_status_adjusted"]]$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        p_FDR_MS <- ggplot(DTE[["MS_status_adjusted"]]$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_MS,p_FDR_MS,volcano_MS),ncol=3,top="MSS-MSI")
	DTE[["MS_status_adjusted"]]$table$DGE_res <- DGE_signif_res[rownames(DTE[["MS_status_adjusted"]]$table),"MS_status_adjusted"]
ggplot(DTE[["MS_status_adjusted"]]$table,aes(x=logFC,y=-log(FDR),col=as.character(DGE_res)))+geom_point(size=0.5)+scale_color_manual(values=c("-1"="lightblue","0"="black","1"="red"),name="DGE\nresults")+ggtitle("MSI / MSS") -> volcano_int_MS
	
	
		

## CIMP DTE

	sel <- pData(isoforms_cell_bank$counts_75)$CIMP_status[which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")]
	group <- factor(ifelse(sel=="CIMP-H","CIMP-H","rest"))
	counts_matrix <- exprs(isoforms_cell_bank$counts_75)[,which(pData(isoforms_cell_bank$counts_75)$CIMP_status!="0")]
	CIMP <- custom_DTE(counts_matrix,group)
	DTE[["CIMP_adjusted"]] <- topTags(CIMP,n=dim(CIMP$table)[1])
	ggplot(DTE[["CIMP_adjusted"]]$table,aes(x=logFC,y=-log(FDR)))+geom_point(size=0.5) -> volcano_CIMP
	p_FC_CIMP <- ggplot(DTE[["CIMP_adjusted"]]$table, aes(x=logFC)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.5,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
        p_FDR_CIMP <- ggplot(DTE[["CIMP_adjusted"]]$table, aes(x=FDR)) +
                geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                binwidth=.01,
                colour="black", fill="white") +
                geom_density(alpha=.2, fill="#FF6666")
	grid.arrange(grobs=list(p_FC_CIMP,p_FDR_CIMP,volcano_CIMP),ncol=3,top="CIMP-H - CIMP-L")
	DTE[["CIMP_adjusted"]]$table$DGE_res <- DGE_signif_res[rownames(DTE[["CIMP_adjusted"]]$table),"CIMP_adjusted"]
ggplot(DTE[["CIMP_adjusted"]]$table,aes(x=logFC,y=-log(FDR),col=as.character(DGE_res)))+geom_point(size=0.5)+scale_color_manual(values=c("-1"="lightblue","0"="black","1"="red"),name="DGE\nresults")+ggtitle("CIMP-H / CIMP-L") -> volcano_int_CIMP
	


###out
	
	dev.off()
	pdf("out/DTE_DGE_comparison.pdf",width=16,height=12)
	grid.arrange(grobs=p_volcano_int_cris, ncol = 3, nrow=3)
	grid.arrange(grobs = p_volcano_int_cms, ncol = 3, nrow=3)
	grid.arrange(grobs = p_volcano_int_mut, ncol = 3, nrow=3)
	grid.arrange(grobs=list(volcano_int_cetuximab,volcano_int_MS,volcano_int_CIMP),ncol=3,nrow=3)	
	dev.off()
	save(DTE, file="out/DTE_results.rda")
	save(DGE_signif_res,file="out/DGE_signif_res.rda")
