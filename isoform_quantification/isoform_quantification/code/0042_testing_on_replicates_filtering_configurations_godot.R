###loading
setwd("out")
Sample_Replicates <- read.table("replicates.txt",sep="\t",header=T)

###code
cutoffs <- c(seq(0.50,0.95,0.05),0.99)
index <- c("IntronDepth","ExonToIntronReadsLeft","ExonToIntronReadsRight")
for(a in 1:length(index)){
	
pdf(paste("IRFinder/testing_on_replicates/selecting",index[a],".pdf",sep=""))
	for(c in 1:dim(Sample_Replicates)[1]){
		tab_S <- read.table(paste("IRFinder/Results/",Sample_Replicates[c,1],"/IRFinder-IR-nondir-val.txt",sep=""),sep="\t",header=T)
		tab_R <- read.table(paste("IRFinder/Results/",Sample_Replicates[c,2],"/IRFinder-IR-nondir-val.txt",sep=""),sep="\t",header=T)
		TP_FP <- vector()
		for(i in 1:length(cutoffs)){
			S <- which(which(tab_S$CNN_IRscore > 0.9 & tab_S$IRratio > 0.3 & tab_S[,index[a]] > quantile(tab_S[,index[a]],cutoffs[i]) & tab_S[,setdiff(index,index[a])[1]] > quantile(tab_S[,setdiff(index,index[a])[1]],cutoffs[1]) & tab_S[,setdiff(index,index[a])[2]] > quantile(tab_S[,setdiff(index,index[a])[2]],cutoffs[1]))%in%grep("clean",tab_S[,4])==TRUE)
			R <- which(which(tab_R$CNN_IRscore > 0.9 & tab_R$IRratio > 0.3 & tab_R[,index[a]] > quantile(tab_R[,index[a]],cutoffs[i]) & tab_R[,setdiff(index,index[a])[1]] > quantile(tab_R[,setdiff(index,index[a])[1]],cutoffs[1]) & tab_R[,setdiff(index,index[a])[2]] > quantile(tab_R[,setdiff(index,index[a])[2]],cutoffs[1]))%in%grep("clean",tab_S[,4])==TRUE)
			nam_S <- paste(tab_S[S,1],tab_S[S,2],tab_S[S,3],tab_S[S,4],tab_S[S,6],sep=":")
			nam_R <- paste(tab_R[R,1],tab_R[R,2],tab_R[R,3],tab_R[R,4],tab_R[R,6],sep=":")
			TP <- length(which(nam_S%in%nam_R==TRUE))
			TP_FP[i] <- TP/(min(length(S),length(R)))
		}
		plot(cutoffs,TP_FP,main=Sample_Replicates[c,1],xlab="Percentile cutoff",ylab="TP/FP ratio")
	}
	dev.off()
} 
setwd(../)
