###loading
        load("data/total_isoforms.rda")
        load("out/CRIS_class.rda")
###code
        total_isoforms <- total_isoforms[,names(cla)]
	perc <- apply(total_isoforms,1,function(x){quantile(x)[3]})
	total_isoforms <- total_isoforms[which(perc>0),]
	cclass <- c("CRISA","CRISB","CRISC","CRISD","CRISE")
        differential_transcripts <- as.list(cclass)
	names(differential_transcripts) <- cclass 
	for (i in cclass){
		#sel <- c(cla[which(cla==paste("only",i,sep=""))],cla[grep(paste("primary",i,sep=""),cla)])
		sel <- cla[which(cla==paste("only",i,sep=""))]
		not <- cla[-grep(i,cla)]
		class <- c(rep(1,length(sel)),rep(0,length(not)))
		p_val <- apply(total_isoforms[,names(c(sel,not))],1,function(x){if(length(unique(x))>1){summary(glm(x ~ class,family="poisson"))$coefficients[8]}else{"NA"}})
		ES <- apply(total_isoforms,1,function(x){log2((mean(x[names(sel)])/mean(x[names(not)])))})
		differential_transcripts[[i]] <- data.frame(ES=ES,p_val=p_val)
		rownames(differential_transcripts[[i]]) <- names(ES)
	}
	save(differential_transcripts,file="out/differential_transcripts.rda")


