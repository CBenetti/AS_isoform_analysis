###loading
setwd("out")
load("IRFinder/IR_score.rda")
load("IRFinder/IR_ratio.rda")
Sample_Replicates <- read.table("replicates.txt",sep="\t",header=T)

###code
IR_score_clean <- IR_score[grep("clean",rownames(IR_score)),]
IR_ratio_clean <- IR_ratio[grep("clean",rownames(IR_ratio)),]
load("IRFinder/v1/IR_score.rda")
load("IRFinder/v1/IR_ratio.rda")
IR_score_clean_old_alignment <- IR_score[grep("clean",rownames(IR_score)),colnames(IR_score_clean)]
IR_ratio_clean_old_alignment <- IR_ratio[grep("clean",rownames(IR_ratio)),colnames(IR_ratio_clean)]

###plot
##comparison within same alignment
pdf("IRFinder/testing_on_replicates/reproducibility_plots_44.pdf")
for(i in 1:(dim(IR_score_clean)[2]/2)){
	plot(IR_score_clean[,Sample_Replicates[i,1]],IR_score_clean[,Sample_Replicates[i,2]],col=ifelse(IR_ratio_clean[,Sample_Replicates[i,1]]>0.3 & IR_ratio_clean[,Sample_Replicates[i,2]]>0.3,"red","black"),pch=ifelse(IR_ratio_clean[,Sample_Replicates[i,1]]>0.3 | IR_ratio_clean[,Sample_Replicates[i,2]]>0.3,"O","."),xlab=paste(Sample_Replicates[i,1],"score",sep=" "),ylab=paste(Sample_Replicates[i,2],"score",sep=" "))
}
dev.off()

##comparison with different alignments
pdf("IRFinder/testing_on_replicates/reproducibility_plots_44_vs_27.pdf")
index <- rownames(IR_score_clean)[which(rownames(IR_score_clean)%in%rownames(IR_score_clean_old_alignment)==TRUE)]
for(i in 1:dim(IR_score_clean)[2]){
	plot(IR_score_clean[index,i],IR_score_clean_old_alignment[index,i],main=paste(colnames(IR_score_clean)[i],"score",sep=" "),xlab="Genecode 44 - ribo excluded",ylab="Genecode 27")
}
dev.off()

##comparison of length within the selected
pdf("IRFinder/testing_on_replicates/length_replicates.pdf")
for(i in 1:(dim(IR_score_clean)[2]/2)){
	par(mfrow=c(1,2))
	sel <- rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0.9 & IR_score_clean[,Sample_Replicates[i,2]] > 0.9 & IR_ratio_clean[,Sample_Replicates[i,1]]>0.3 & IR_ratio_clean[,Sample_Replicates[i,2]] > 0.3)]
	false_positive_S <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0.9 & IR_ratio_clean[,Sample_Replicates[i,1]]>0.3)],sel)
	all_S <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0)],sel)
	false_positive_R <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,2]]>0.9 & IR_ratio_clean[,Sample_Replicates[i,2]]>0.3)],sel)
	all_R <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,2]]>0)],sel)
	tab <- read.table(paste("IRFinder/results/",Sample_Replicates[i,1],"/IRFinder-IR-nondir-val.txt",sep=""),header=T,sep="\t")
	rownames(tab) <- gsub("/",":",paste(tab[,1],tab[,2],tab[,3],tab[,4],tab[,6],sep=":"))
	tab_R <- read.table(paste("IRFinder/results/",Sample_Replicates[i,2],"/IRFinder-IR-nondir-val.txt",sep=""),header=T,sep="\t")
	rownames(tab_R) <- gsub("/",":",paste(tab_R[,1],tab_R[,2],tab_R[,3],tab_R[,4],tab_R[,6],sep=":"))
	boxplot(tab[sel,3]-tab[sel,2],tab[false_positive_S,3]-tab[false_positive_S,2],tab_R[sel,3]-tab_R[sel,2],tab_R[false_positive_R,3]-tab_R[false_positive_R,2],at=c(1,2,4,5),col=c("goldenrod1","goldenrod1","yellowgreen","yellowgreen"),names=c("TP-Sample","FP-Sample","TP-Rep","FP-Rep"),main=Sample_Replicates[i,1])
	boxplot(tab[sel,3]-tab[sel,2],tab[all_S,3]-tab[all_S,2],tab[sel,3]-tab[sel,2],tab[all_R,3]-tab[all_R,2],at=c(1,2,4,5),col=c("goldenrod1","goldenrod1","yellowgreen","yellowgreen"),names=c("TP-Sample","all-Sample","TP-Rep","all-Rep"),main=Sample_Replicates[i,1])
}
dev.off()

##comparison of depth within the selected
pdf("IRFinder/testing_on_replicates/depth_replicates.pdf")
for(i in 1:(dim(IR_score_clean)[2]/2)){
        par(mfrow=c(1,2))
	sel <- rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0.9 & IR_score_clean[,Sample_Replicates[i,2]] > 0.9 & IR_ratio_clean[,Sample_Replicates[i,1]]>0.3 & IR_ratio_clean[,Sample_Replicates[i,2]] > 0.3)]
	false_positive_S <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0.9 & IR_ratio_clean[,Sample_Replicates[i,1]]>0.3)],sel)
        all_S <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0)],sel)
        false_positive_R <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,2]]>0.9 & IR_ratio_clean[,Sample_Replicates[i,2]]>0.3)],sel)
        all_R <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,2]]>0)],sel)
        tab <- read.table(paste("IRFinder/results/",Sample_Replicates[i,1],"/IRFinder-IR-nondir-val.txt",sep=""),header=T,sep="\t")
        rownames(tab) <- gsub("/",":",paste(tab[,1],tab[,2],tab[,3],tab[,4],tab[,6],sep=":"))
        tab_R <- read.table(paste("IRFinder/results/",Sample_Replicates[i,2],"/IRFinder-IR-nondir-val.txt",sep=""),header=T,sep="\t")
       	rownames(tab_R) <- gsub("/",":",paste(tab_R[,1],tab_R[,2],tab_R[,3],tab_R[,4],tab_R[,6],sep=":"))
	boxplot(tab[sel,9],tab[false_positive_S,9],tab_R[sel,9],tab_R[false_positive_R,9],at=c(1,2,4,5),col=c("goldenrod1","goldenrod1","yellowgreen","yellowgreen"),names=c("TP-Sample","FP-Sample","TP-Rep","FP-Rep"),main=Sample_Replicates[i,1])
	boxplot(tab[sel,9],tab[all_S,9],tab_R[sel,9],tab_R[all_R,9],at=c(1,2,4,5),col=c("goldenrod1","goldenrod1","yellowgreen","yellowgreen"),names=c("TP-Sample","all-Sample","TP-Rep","all-Rep"),main=Sample_Replicates[i,1])
}
dev.off()


##comparison of ExonToIntronReadsRight within the selected
pdf("IRFinder/testing_on_replicates/ExonToIntronReadsRight.pdf")
for(i in 1:(dim(IR_score_clean)[2]/2)){
        par(mfrow=c(1,2))
	sel <- rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0.9 & IR_score_clean[,Sample_Replicates[i,2]] > 0.9 & IR_ratio_clean[,Sample_Replicates[i,1]]>0.3 & IR_ratio_clean[,Sample_Replicates[i,2]] > 0.3)]
	false_positive_S <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0.9 & IR_ratio_clean[,Sample_Replicates[i,1]]>0.3)],sel)
        all_S <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0)],sel)
        false_positive_R <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,2]]>0.9 & IR_ratio_clean[,Sample_Replicates[i,2]]>0.3)],sel)
        all_R <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,2]]>0)],sel)
        tab <- read.table(paste("IRFinder/results/",Sample_Replicates[i,1],"/IRFinder-IR-nondir-val.txt",sep=""),header=T,sep="\t")
        rownames(tab) <- gsub("/",":",paste(tab[,1],tab[,2],tab[,3],tab[,4],tab[,6],sep=":"))
        tab_R <- read.table(paste("IRFinder/results/",Sample_Replicates[i,2],"/IRFinder-IR-nondir-val.txt",sep=""),header=T,sep="\t")
        rownames(tab_R) <- gsub("/",":",paste(tab_R[,1],tab_R[,2],tab_R[,3],tab_R[,4],tab_R[,6],sep=":"))
	boxplot(tab[sel,14],tab[false_positive_S,14],tab_R[sel,14],tab_R[false_positive_R,14],at=c(1,2,4,5),col=c("goldenrod1","goldenrod1","yellowgreen","yellowgreen"),names=c("TP-Sample","FP-Sample","TP-Rep","FP-Rep"),main=Sample_Replicates[i,1])
	boxplot(tab[sel,14],tab[all_S,14],tab_R[sel,14],tab_R[all_R,14],at=c(1,2,4,5),col=c("goldenrod1","goldenrod1","yellowgreen","yellowgreen"),names=c("TP-Sample","all-Sample","TP-Rep","all-Rep"),main=Sample_Replicates[i,1])
}
dev.off()
   
##comparison of ExonToIntronReadsLeft within the selected
pdf("IRFinder/testing_on_replicates/ExonToIntronReadsLeft.pdf")
for(i in 1:(dim(IR_score_clean)[2]/2)){
        par(mfrow=c(1,2))
	sel <- rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0.9 & IR_score_clean[,Sample_Replicates[i,2]] > 0.9 & IR_ratio_clean[,Sample_Replicates[i,1]]>0.3 & IR_ratio_clean[,Sample_Replicates[i,2]] > 0.3)]
	false_positive_S <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0.9 & IR_ratio_clean[,Sample_Replicates[i,1]]>0.3)],sel)
        all_S <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,1]]>0)],sel)
        false_positive_R <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,2]]>0.9 & IR_ratio_clean[,Sample_Replicates[i,2]]>0.3)],sel)
        all_R <- setdiff(rownames(IR_score_clean)[which(IR_score_clean[,Sample_Replicates[i,2]]>0)],sel)
        tab <- read.table(paste("IRFinder/results/",Sample_Replicates[i,1],"/IRFinder-IR-nondir-val.txt",sep=""),header=T,sep="\t")
        rownames(tab) <- gsub("/",":",paste(tab[,1],tab[,2],tab[,3],tab[,4],tab[,6],sep=":"))
        tab_R <- read.table(paste("IRFinder/results/",Sample_Replicates[i,2],"/IRFinder-IR-nondir-val.txt",sep=""),header=T,sep="\t")
        rownames(tab_R) <- gsub("/",":",paste(tab_R[,1],tab_R[,2],tab_R[,3],tab_R[,4],tab_R[,6],sep=":"))
	boxplot(tab[sel,13],tab[false_positive_S,13],tab_R[sel,13],tab_R[false_positive_R,13],at=c(1,4,7,10),col=c("goldenrod1","goldenrod1","yellowgreen","yellowgreen"),names=c("TP-Sample","FP-Sample","TP-Rep","FP-Rep"),main=Sample_Replicates[i,1])
	boxplot(tab[sel,13],tab[all_S,13],tab_R[sel,13],tab_R[all_R,13],at=c(1,4,7,10),col=c("goldenrod1","goldenrod1","yellowgreen","yellowgreen"),names=c("TP-Sample","all-Sample","TP-Rep","all-Rep"),main=Sample_Replicates[i,1])
}
dev.off()

setwd("../")
