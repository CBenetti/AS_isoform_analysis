###library
        library(cowplot)
        library(gridExtra)
        library(grid)
        library(lattice)
	library(venn)
	library(ComplexHeatmap)
	library(circlize)
	library(Biobase)
	library(enrichR)
	library("RColorBrewer")
	library(ggplot2)
	library(reshape2)
###functions
color_labels <- function(x,tr){
		if(length(unlist(strsplit(x,":")))>1){
			classes <- unlist(strsplit(x,":"))
			tab <- matrix(0,ncol=length(classes),nrow=length(tr)) 
			colnames(tab) <- classes
			rownames(tab) <- tr
			tab[tr[which(tr%in%rownames(DTU_table)==TRUE)],classes] <- as.matrix(DTU_table[tr[which(tr%in%rownames(DTU_table)==TRUE)],classes])
			comb <- expand.grid(rep(list((0:1)), length(classes)))
			comb <- comb[which(rowSums(comb)<3),]
			comb <- comb[which(rowSums(comb)!=0),]
			li <-lapply(c(1:dim(comb)[1]),function(z){comb[z,]})
			cols <- brewer.pal(n = length(li), name = "Paired")
			names(cols) <- lapply(li,function(z){paste(c(classes[which(z!=0)]),collapse=":")}) 
			colfun <- apply(tab,1,function(x){if(length(which(x!=0))>0){cols[paste(c(names(x)[which(x!=0)]),collapse=":")]}else{"black"}})
			list(colfun,cols[which(unlist(cols)%in%colfun==TRUE)])
		}else{colfun <-ifelse(tr%in%rownames(DTU_table)==TRUE,"red","black")
		names(colfun) <- tr
		colfun
		}
	}

###load
	load("out/DIE_results_DRIMSeq_withL2R.rda")
	load("out/annotated_isoforms_cell_bank_proportions.rda")
	load("out/annotated_genes_counts_cell_bank.rda")

###code
	cpm <- apply(exprs(annotated_genes_counts_cell_bank),2,function(x){(x*1000000)/sum(x)})
	rownames(cpm) <- fData(annotated_genes_counts_cell_bank)$Geneid	
	DTU_significative <- lapply(DIE,function(x){
		if(is.data.frame(x)==FALSE){lapply(x,function(y){
				y[which(y$FDR<0.01 & y$proportionSD > 0.1),]
			})
		}else{
			x[which(x$FDR<0.01 & x$proportionSD > 0.1),]
		}
	})
	DTU_significative <- c(DTU_significative[[1]],DTU_significative[[2]],DTU_significative[3],DTU_significative[[4]],DTU_significative[5],DTU_significative[6])
	names(DTU_significative) <- unlist(lapply(names(DIE),function(x){if(is.data.frame(DIE[[x]])==FALSE){names(DIE[[x]])}else{x}}))
	significative_transcripts <- unique(unlist(sapply(DTU_significative,"[[",4)))
	DTU_table <- matrix(0,ncol=length(DTU_significative),nrow=length(significative_transcripts))
	rownames(DTU_table) <- significative_transcripts
	colnames(DTU_table) <- names(DTU_significative)
	for(i in names(DTU_significative)){DTU_table[DTU_significative[[i]]$rownam,i] <- DTU_significative[[i]]$logFC}
	txInfo <- as.data.frame(matrix(data=NA,ncol=2,nrow=length(fData(isoforms_cell_bank$counts_75)$Geneid)))
        colnames(txInfo) <- c("isoform_id","gene_id")
        rownames(txInfo) <- txInfo$isoform_id <- rownames(exprs(isoforms_cell_bank$counts_75))
	txInfo$gene_id <- fData(isoforms_cell_bank$counts_75)$Geneid
	DTU_table <- data.frame(DTU_table,gene_id=txInfo[rownames(DTU_table),"gene_id"])
	signif_genes <- sapply(DTU_significative,"[[",3)
	signif_genes <- lapply(signif_genes,unique)
	names(signif_genes) <- names(DTU_significative)
##intersections
	pdf("out/DTU_genes_intersections.pdf")
	venn(signif_genes[grep("onlyCRIS",names(signif_genes))],ilabels="counts",ellipse=TRUE,zcolor=c("red","lightgreen","lightblue","yellow","blue"))
	venn(signif_genes[grep("CMS",names(signif_genes))],ilabels="counts",ellipse=TRUE,zcolor=c("yellow","steelblue","purple","grey"))
	venn(signif_genes[grep("RA",names(signif_genes))],ilabels="counts",ellipse=TRUE,zcolor=c("orange","brown"))
	dev.off()
#ontologies of unique genes per class (B-D-E)
	info <- extractInfo(signif_genes[grep("CRIS",names(signif_genes))],what="intersections",use.names=TRUE)
	le <- lapply(info,length)
	info_onthology <- c(info[which(le>5)])
	dbs <- c("MSigDB_Hallmark_2020","GO_Molecular_Function_2015","GO_Cellular_Component_2015","GO_Biological_Process_2015","KEGG_2021_Human",
"BioCarta_2015","Panther_2015","WikiPathways_2015","Reactome_2015","ENCODE_TF_ChIP-seq_2015")
	onthology_res_full <- lapply(signif_genes[grep("CRIS",names(signif_genes))],function(x){enrichr(x,dbs)})
	names(onthology_res_full) <- paste(names(onthology_res_full),"full",sep="_")
	onthology_res <- lapply(info_onthology,function(x){enrichr(x,dbs)})
	c(onthology_res,onthology_res_full) -> onthology_res
	onthology_res_signif <- lapply(onthology_res,function(x){
		lapply(x,function(y){if(dim(y)[1]>0){
				int <- as.numeric(unlist(sapply(strsplit(y$Overlap,"/"),"[[",1)))
				y[which(y$Adjusted.P.value<0.05 & y$Odds.Ratio > 1 & int>2),]
			}else{y}
		})
	})	
	pdf("out/onthology.pdf",width=18, height=12)
	lapply(names(onthology_res_signif),function(x){
		p <- lapply(names(onthology_res_signif[[x]]),function(y){
			if(dim(onthology_res_signif[[x]][[y]])[1]>0){
				plotEnrich(onthology_res_signif[[x]][[y]], y = "Count", orderBy = "P.value")+ggtitle(y)	
			}else{"NA"}
		})
		grid.arrange(grobs=p[which(p!="NA")],ncol=3,nrow=3,top=x)
	})
	dev.off()
#gene_level_plots_CRIS
	group <- pData(isoforms_cell_bank$isoform_proportions)$CRIS_class[grep("only",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)]
	cells <- grep("only",pData(isoforms_cell_bank$isoform_proportions)$CRIS_class)
	group_colors <- c("red","darkgreen","lightblue","orange","blue")
	names(group_colors) <- unique(group)
	lapply(names(info),function(x){
		pdf(paste("out/gene_plots/",x,".pdf",sep=""),width=15,height=8)
		lapply(info[[x]],function(y){
			proportions <- exprs(isoforms_cell_bank$isoform_proportions)[which(fData(isoforms_cell_bank$isoform_proportions)$Geneid==y),cells]
			ma <- apply(proportions,1,max)
			proportions <- proportions[which(ma>0.10),]
			prop_samp <- data.frame(feature_id= rownames(proportions),proportions)
			oo <- order(apply(aggregate(t(prop_samp[, -1]),by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)
    			feature_levels <- rownames(prop_samp)[oo]  
			sample_levels <- colnames(proportions)[order(group)]
			prop_samp <- melt(prop_samp, id.vars = "feature_id", variable.name = "sample_id", value.name = "proportions")
			prop_samp$feature_id <- factor(prop_samp$feature_id, levels = feature_levels)
			prop_samp$group <- rep(group, each = nrow(proportions))
			prop_samp$sample_id <- factor(prop_samp$sample_id, levels = sample_levels)
			prop_samp$feature_id <- factor(prop_samp$feature_id, levels = feature_levels)
			prop_samp$group <- rep(group, each = nrow(proportions))
			prop_samp$sample_id <- factor(prop_samp$sample_id, levels = sample_levels)
			colfun  <- color_labels(x,rownames(proportions))
			if(is.list(colfun)==TRUE){
				feature_colors <- colfun[[1]][feature_levels]
			}else{
				feature_colors <- colfun[feature_levels]
			}
			g <- ggplot() +
			geom_jitter(data = prop_samp, aes_string(x = "feature_id", 
      			y = "proportions", fill = "group", colour = "group"), 
      			position = position_jitterdodge(), 
      			alpha = 0.9, size = 2, show.legend = FALSE, na.rm = TRUE) +
    			geom_boxplot(data = prop_samp, aes_string(x = "feature_id", 
      			y = "proportions", colour = "group", fill = "group"), 
      			outlier.size = NA, alpha = 0.4, lwd = 0.5) +
    			theme_bw() + 
    			theme(axis.text.x = element_text(angle = 90, vjust = 0.5,colour=feature_colors), 
      			axis.text=element_text(size=16), 
      			axis.title=element_text(size=14, face="bold"), 
      			plot.title = element_text(size=16), 
      			legend.position = "right", 
      			legend.title = element_text(size = 14), 
      			legend.text = element_text(size = 14)) +
    			ggtitle(paste(y,"; mean cpm expression: ",mean(cpm[y,]),sep="")) +     
    			scale_fill_manual(name = "Groups", values = group_colors, 
      			breaks = names(group_colors)) +
    			scale_colour_manual(name = "Groups", values = group_colors, 
      			breaks = names(group_colors)) +
    			xlab("Features") +
    			ylab("Proportions")
			print(g)
			if(is.list(colfun)==TRUE){		
				plot.new()
				legend("topleft",legend=names(colfun[[2]]),pch=18,col=colfun[[2]],title=paste(x,y,sep=":"))
			}
		})
		dev.off()
	})

##heatmaps
	shared_genes <- names(info)[grep(":",names(info))]
        pdf("out/transcripts_in_shared_genes.pdf")
        lapply(shared_genes,function(x){
                classes <- unlist(strsplit(x,":"))
                DTU_tab <- DTU_table[which(DTU_table$gene_id%in%info[[x]]==TRUE),classes]
                genes <- DTU_table$gene_id[which(DTU_table$gene_id%in%info[[x]]==TRUE)]  
                DTU_tab <- DTU_tab[order(genes),]
                genes <- genes[order(genes)]
                draw(Heatmap(as.matrix(DTU_tab),cluster_rows=FALSE,
                        cluster_columns=FALSE,
                        name="L2R\nisoform\nproportions",
                        col=colorRamp2(c(0.5,0,-0.5),c("red","white","blue")),
top_annotation=HeatmapAnnotation(classes=anno_block(labels=classes,gp=gpar(fill=group_colors[classes]))),
                        left_annotation=rowAnnotation(genes=anno_block(labels=unique(genes),labels_gp=gpar(fontsize=if(length(genes)<10){10}else{5}))),
                        row_split=genes,show_column_names=FALSE,column_split=classes,  
                        row_names_side = "right",row_title=NULL,
                column_title=paste(c(classes,"shared genes"),collapse=" ")))
        })
        dev.off()

###out
	save(info, file="out/DTU_intersections_CRIS.rda")
	save(DTU_significative,file="out/DTU_significative.rda")
	write.table(DTU_table,"out/DTU_table.txt",sep="\t",col.names=T,row.names=T)
	
