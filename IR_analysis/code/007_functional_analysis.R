###library
	library(rWikiPathways)
	library(clusterProfiler)
	library(ggplot2)
###load
	load("out/differential_transcripts.rda")
	load("data/total_isoforms.rda")
###code
	paths <- readPathwayGMT("data/wikipathways-20240311-gmt-Homo_sapiens.gmt")
	U <- bitr(rownames(total_isoforms), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
	paths <- paths[which(paths$gene%in%U==TRUE),]
	signif <- lapply(differential_transcripts,function(x){rownames(x)[which(x$ES>1 & x$p_val<0.05)]})
	entrez <- lapply(signif,function(x){bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID})
	path_names <- unique(paths$name)
	for(i in names(signif)){
		contingency_tables <- lapply(path_names,function(x){cbind(c(length(which(entrez[[i]]%in%paths[which(paths$name==x),"gene"]==TRUE)),length(which(entrez[[i]]%in%paths[which(paths$name==x),"gene"]==FALSE))),c(length(which(setdiff(U,entrez[[i]])%in%paths[which(paths$name==x),"gene"]==TRUE)),length(which(setdiff(U,entrez[[i]])%in%paths[which(paths$name==x),"gene"]==FALSE))))})
		names(contingency_tables) <- path_names
		ES <- lapply(path_names, function(x){log2((contingency_tables[[x]][1,1]/sum(contingency_tables[[x]][,1]))/(length(paths[which(paths$name==x),"gene"])/length(U)))})
		p_value <- lapply(contingency_tables,function(x){fisher.test(x)$p.value})
		le_count <- as.numeric(unlist(lapply(path_names,function(x){length(which(entrez[[i]]%in%paths[which(paths$name==x),"gene"]==TRUE))})))
		le <- lapply(path_names,function(x){paste(c(entrez[[i]][which(entrez[[i]]%in%paths[which(paths$name==x),"gene"]==TRUE)]),collapse=";")})
		le[which(le_count < 1)] <- "NA"
		functional_results <- data.frame(ES=unlist(ES),p_val=unlist(unname(p_value)),hits=unlist(unname(le)),hits_n=unname(le_count),path=path_names)
		rownames(functional_results) <- path_names
		write.table(functional_results,paste("out/wikiPathways_",i,"_enrichment_results.txt",sep=""),sep="\t")
		pdf(paste("out/WikiPathways_",i,".pdf",sep=""))
		plot(ggplot(functional_results[which(ES>1 & p_value<0.05),], aes(x=ES, y=path, colour = p_val)) + scale_colour_gradient(low = "red", high = "blue")+
 		geom_point(aes(size=hits_n))+ggtitle("WikiPathways")+labs(y = ""))		
		dev.off()
	}
