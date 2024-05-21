###library

###loading
	coverage <- read.table("data/coverage.txt", sep="\t")
	rownames(coverage) <- coverage[,1]
	load("out/CRIS_class.rda")
	load("data/IR_ratio.rda")
# Define the pairs of samples and replicates
	samples <- c("C10.STAR","HDC114.STAR", "LS180.STAR", "LS411N.STAR","SNU1040.STAR", "SNU1235.STAR", "SNU81.STAR", "SW620.STAR")
	replicates <- c("C10S1.STAR", "HDC114S17.STAR","LS180S23.STAR", "LS411NS23.STAR", "SNU1040S11.STAR", "SNU1235S21.STAR", "SNU81S15.STAR", "SW620S14.STAR")

###code
	cutoffs<- c(0.05,0.1,0.2,0.3)
	for(c in cutoffs){
		load(paste("data/count_raw_results_cutoff_",c,".rda",sep=""))
		intron_tot <- apply(IR_ratio,2,function(x){length(intersect(which(x>c),grep("clean",rownames(IR_ratio))))})
          	intron_tot_n <- length(unique(unlist(apply(IR_ratio,2,function(x){intersect(which(x>c),grep("clean",rownames(IR_ratio)))}))))
		n_i <- apply(IR_ratio,1,function(x){length(which(x>0))})
                unique_i <- apply(IR_ratio[intersect(which(n_i==1),grep("clean",rownames(IR_ratio))),],2,function(x){length(which(x>c))})
		ni_tot <- dim(raw_tab)[1]
		ni_tot_cris <- apply(raw_tab[,names(cla)],1,sum)
		ni_tot_cris <- ni_tot_cris[which(ni_tot_cris>0)]
		n_positive_samples <- apply(raw_tab,1,function(x){length(which(x>0))})
		unique_tab <- raw_tab[which(n_positive_samples==1),]
		ni_unique <- dim(unique_tab)[1]
		ni_unique_cris <- apply(unique_tab[,names(cla)],1,sum)
                ni_unique_cris <- ni_unique_cris[which(ni_unique_cris>0)]
		n_i_samples <- apply(raw_tab,2,function(x){length(which(x>0))})
		n_unique_i_samples <- apply(unique_tab,2,function(x){length(which(x>0))})
		n_i_samples_cris <- n_i_samples[names(cla)]
		n_unique_i_samples_cris <- n_unique_i_samples[names(cla)]
	##plot
		pdf(paste("out/report_of_pipeline_cutoff_",c,".pdf",sep=""),width=15)
		df <- data.frame(type=factor(c(rep("IR Finder",length(n_i_samples)),rep("Stringtie - IRFinder",length(n_i_samples))),levels=c("IR Finder","Stringtie - IRFinder")),
	num= c(intron_tot[names(n_i_samples)]-(n_i_samples),n_i_samples),
	cells=factor(rep(names(n_i_samples),2),levels=names(sort(n_i_samples))))
		plot(ggplot(df,aes(y=num,x=cells,fill=type))+ 
geom_bar(stat="identity",alpha=ifelse(df$cells%in%names(cla)==TRUE,1,0.5))+theme(axis.text.x=element_blank())+ggtitle(paste("Intron retention events with at least",c,"coverage",sep=" ")))

##In sample and replicates
		plot(n_i_samples[samples],xlab="samples",ylab="number of introns > 0 detected",col="black",pch=18,xaxt="n", main="Intersected introns",ylim=c(0,max(n_i_samples)))
		points(n_i_samples[replicates],col="darkgreen", pch=18)
		legend("topleft",pch=c(18,18),col=c("black","darkgreen"),legend=c("Sample","Replicate"))
		
        	plot(sort(n_unique_i_samples),xlab="samples",ylab="number of unique introns detected",pch=18,xaxt="n",main=paste("Total:",ni_unique,"introns and",length(ni_unique_cris),"for CRIS classified",sep=" "),col=ifelse(names(sort(n_unique_i_samples))%in%names(n_i_samples_cris)==TRUE,"red","black"))
        	axis(2,at=sort(c(0,50,100,150,200,250,300,350,max(n_unique_i_samples),min(n_unique_i_samples))))
        	abline(h=max(n_unique_i_samples),col="darkgreen")
        	abline(h=min(n_unique_i_samples),col="darkred")
	
		par(mfrow=c(1,2))
		plot(n_i_samples,n_unique_i_samples,xlab="number of introns >0",ylab="number of unique introns",col=ifelse(names(n_i_samples)%in%names(n_i_samples_cris)==TRUE,"red","black"),main=paste("Pearson correlation between number of total introns \n  and unique introns per sample:",cor(n_i_samples,n_unique_i_samples,method="pearson"),"and\n",cor(n_i_samples_cris,n_unique_i_samples_cris,method="pearson"),"for those classified for CRIS",sep=" "))

		plot(n_i_samples,n_unique_i_samples,xlab="number of introns >0",ylab="number of unique introns",col=ifelse(names(n_i_samples)%in%names(n_i_samples_cris)==TRUE,"red","black"),main=paste("Spearman correlation between number of total introns \n and unique introns per sample:",cor(n_i_samples,n_unique_i_samples,method="spearman"),"and\n",cor(n_i_samples_cris,n_unique_i_samples_cris,method="spearman"),"for those classified for CRIS",sep=" "))

		par(mfrow=c(1,2))
		plot(n_i_samples,coverage[names(n_i_samples),2],xlab="number of introns >0",col=ifelse(names(n_i_samples)%in%names(n_i_samples_cris)==TRUE,"#00AFBB","black"),pch=20,ylab="number of assigned reads per sample",main=paste("Pearson correlation between number \n of total introns and assigned reads:",cor(n_i_samples,coverage[names(n_i_samples),2],method="pearson"),"and\n",cor(n_i_samples_cris,coverage[names(n_i_samples_cris),2],method="pearson"),"for those classified for CRIS",sep=" "))
	
		plot(n_i_samples,coverage[names(n_i_samples),2],xlab="number of introns >0",col=ifelse(names(n_i_samples)%in%names(n_i_samples_cris)==TRUE,"#00AFBB","black"),pch=20,ylab="number of assigned reads per sample",main=paste("Spearman correlation between number \n of total introns and assigned reads:",cor(n_i_samples,coverage[names(n_i_samples),2],method="spearman"),"and\n",cor(n_i_samples_cris,coverage[names(n_i_samples_cris),2],method="spearman"),"for those classified for CRIS",sep=" "))

		par(mfrow=c(1,2))
		plot(n_unique_i_samples,coverage[names(n_unique_i_samples),2],xlab="number of unique introns",col=ifelse(names(n_i_samples)%in%names(n_i_samples_cris)==TRUE,"red","black"),ylab="number of assigned reads per sample",main=paste("Pearson correlation between number of total \n introns and assigned reads:",cor(n_unique_i_samples,coverage[names(n_unique_i_samples),2],method="pearson"),"and\n",cor(n_unique_i_samples_cris,coverage[names(n_unique_i_samples_cris),2],method="pearson"),"for those classified for CRIS",sep=" "))

		plot(n_unique_i_samples,coverage[names(n_unique_i_samples),2],xlab="number of unique introns",ylab="number of assigned reads per sample",main=paste("Spearman correlation between number of total \n introns and assigned reads:",cor(n_unique_i_samples,coverage[names(n_unique_i_samples),2],method="spearman"),"and\n",cor(n_unique_i_samples_cris,coverage[names(n_unique_i_samples_cris),2],method="spearman"),"for those classified for CRIS",sep=""),col=ifelse(names(n_i_samples)%in%names(n_i_samples_cris)==TRUE,"red","black"))

		dev.off()
	}
