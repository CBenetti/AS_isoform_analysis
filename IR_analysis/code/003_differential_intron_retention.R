###library
	library(xlsx)
	library(ggplot2)
###loading
	load("out/CRIS_D_classification.rda")
	load("out/CRIS_class.rda")
	load("data/IR_ratio.rda")
###code
	cutoffs <- c(0.05,0.1,0.2,0.3)
	cclass <- c("CRISA","CRISB","CRISC","CRISD","CRISE")
	for(c in cutoffs){
		load(paste("data/count_raw_results_cutoff_",c,".rda",sep=""))
		intron_tot <- apply(IR_ratio,2,function(x){length(intersect(which(x>c),grep("clean",rownames(IR_ratio))))})
                n_i <- apply(IR_ratio,1,function(x){length(which(x>0))})
                unique_i <- apply(IR_ratio[intersect(which(n_i==1),grep("clean",rownames(IR_ratio))),],2,function(x){length(which(x>c))})
		n_positive_samples <- apply(raw_tab,1,function(x){length(which(x>0))})
		unique_tab <- raw_tab[which(n_positive_samples==1),]
	        n_i_samples <- apply(raw_tab,2,function(x){length(which(x>0))})
	        n_unique_i_samples <- apply(unique_tab,2,function(x){length(which(x>0))})
		count_s <- apply(raw_tab,2,sum)
		for(i in cclass){
			 not <- cla[-grep(i,cla)]
                        class_spec <- apply(raw_tab[,names(cla)],1,function(x){if(sum(x[names(not)])==0){length(which(x[names(cla)[which(cla==paste("only",i,sep=""))]]>0))}else{0}})
                        sel <- c(cla[which(cla==paste("only",i,sep=""))],cla[grep(paste("primary",i,sep=""),cla)])
                        primary_spec <- apply(raw_tab[,names(cla)],1,function(x){if(sum(x[names(not)])==0){length(which(x[names(sel)]>0))}else{0}})
                        df_only <-  
data.frame(intron=class_spec[which(class_spec>0)],samples=unlist(lapply(names(class_spec[which(class_spec>0)]),
function(x){paste(c(names(which(raw_tab[x,names(cla)[which(cla==paste("only",i,sep=""))]]>0))),collapse=";")})))
                        df_primary <- data.frame(intron=primary_spec[which(primary_spec>0)],samples=unlist(lapply(names(primary_spec[which(primary_spec>0)]),
function(x){paste(c(names(which(raw_tab[x,names(sel)]>0))),collapse=";")})))
                       # write.xlsx(df_only,file=paste("out/intron_list_only_cris_cutoff_",c,".xlsx",sep=""),sheetName=i,append=TRUE)
                       # write.xlsx(df_primary,file=paste("out/intron_list_primary_cris_cutoff_",c,".xlsx",sep=""),sheetName=i,append=TRUE)

		}
		Single_Cris <- cla[grep("only",cla)]
	##plot
		df <- data.frame(class=Single_Cris,intron_number=intron_tot[names(Single_Cris)])
		pdf(paste("out/Discrete_number_class_distribution_cutoff_",c,".pdf",sep=""))
		plot(ggplot(df,aes(x=class,y=intron_number,fill=class))+	
		geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
		scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
		scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
		ggtitle("Distribution of total intron for sample in each class"))
		df <- data.frame(class=Single_Cris,intron_number=n_i_samples[names(Single_Cris)])
		plot(ggplot(df,aes(x=class,y=intron_number,fill=class))+ 
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
		scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Distribution of selected intron for sample in each class"))
		#boxplot(intron_tot[names(Single_Cris)] ~ Single_Cris, xlab="CRIS class",ylab="total IRFinder clean intron number", main="Distribution of total intron \n for sample in each class", col=c("red","darkgreen","lightblue","orange","blue"))
		#boxplot(n_i_samples[names(Single_Cris)] ~ Single_Cris, xlab="CRIS class",ylab="total intron number", main="Distribution of total intron \n for sample in each class", col=c("red","darkgreen","lightblue","orange","blue"))
		boxplot(n_unique_i_samples[names(Single_Cris)] ~ Single_Cris, xlab="CRIS class",ylab="unique intron number", main="Distribution of unique intron \n for sample in each class", col=c("red","darkgreen","lightblue","orange","blue"))
		boxplot(count_s[names(Single_Cris)] ~ Single_Cris, xlab="CRIS class",ylab="total intron count", main="Distribution of total \n intron count in each class", col=c("red","darkgreen","lightblue","orange","blue"))
		dev.off()
	}
