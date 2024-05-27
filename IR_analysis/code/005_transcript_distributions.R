###library
	library(MASS)
	library(vioplot)
	library(tidyr)
	library(ggplot2)	
###loading
	load("data/total_isoforms.rda")
	load("out/CRIS_class.rda")
	summary <- read.table("data/summary.txt",sep="\t",header=T)
	long_df <- pivot_longer(summary,cols=-sample,names_to="transcript_type",values_to="number")
	#rownames(summary)<- summary[,1]
	#summary <- summary[,-1]
	coverage<-read.table("data/coverage.txt")
	rownames(coverage) <- coverage[,1]
	colnames(coverage) <- c("sample","coverage")
###code
	pdf("out/transcript_quantitative_distribution.pdf")
	df <- long_df[which(long_df$transcript_type%in%c("total_transcripts","assigned_annotated")),]
	df[which(df$transcript_type=="total_transcripts"),"number"] <- 
df[which(df$transcript_type=="total_transcripts"),"number"]-df[which(df$transcript_type=="assigned_annotated"),"number"]
plot(ggplot(df,aes(x=factor(sample,levels=summary$sample[order(summary$total_transcripts)]),y=number,fill=factor(transcript_type,levels=c("total_transcripts","assigned_annotated"))))+
geom_bar(stat="identity",alpha=ifelse(df$sample%in%names(cla)==TRUE,1,0.5))+guides(fill=guide_legend(title=""))+
xlab("sample")+theme(axis.text.x=element_blank())+ggtitle("Transcript quantitative evaluation"))

df <- long_df[which(long_df$transcript_type%in%c("total_transcripts","assigned_annotated","assigned_unannotated")),]
        df[which(df$transcript_type=="total_transcripts"),"number"] <-
df[which(df$transcript_type=="total_transcripts"),"number"]-df[which(df$transcript_type=="assigned_annotated"),"number"]
        df[which(df$transcript_type=="assigned_annotated"),"number"] <-
df[which(df$transcript_type=="assigned_annotated"),"number"]-df[which(df$transcript_type=="assigned_unannotated"),"number"]
plot(
ggplot(df,aes(x=factor(sample,levels=summary$sample[order(summary$total_transcripts)]),y=number,fill=factor(transcript_type,levels=c("total_transcripts","assigned_annotated","assigned_unannotated"))))+
geom_bar(stat="identity",alpha=ifelse(df$sample%in%names(cla)==TRUE,1,0.5))+guides(fill=guide_legend(title=""))+
xlab("sample")+theme(axis.text.x=element_blank())+ggtitle("Transcript quantitative evaluation"))

	dev.off()
	Single_Cris <- data.frame(class=cla[grep("onlyCRIS",cla)],sample=names(cla)[grep("onlyCRIS",cla)])
	Primary<- sapply(strsplit(cla,"secondary"),"[[",1)
	Primary_Cris <- data.frame(class=Primary,sample=names(Primary))
	Secondary<-lapply(c("CRISA","CRISB","CRISC","CRISD","CRISE"),function(x){
		sec <- data.frame(names(cla)[setdiff(grep(x,cla),grep(x,Primary_Cris))],rep(paste("secondary",x,sep=""),length(setdiff(grep(x,cla),grep(x,Primary_Cris)))))
	})
	do.call(rbind,Secondary)-> Secondary_Cris
	colnames(Secondary_Cris) <- c("sample","class")
	Primary_Cris <- Primary_Cris[which(Primary_Cris$class!="NA"),]
	Single_Cris_summary <- merge(Single_Cris,summary,"sample")
	Primary_Cris_summary <- merge(Primary_Cris,summary,"sample")
	Secondary_Cris_summary <- merge(Secondary_Cris,summary,"sample")
	tot_summary <- rbind(Primary_Cris_summary,Secondary_Cris_summary)

	pdf("out/transcript_distribution.pdf")
 plot(ggplot(Single_Cris_summary,aes(x=class,y=total_transcripts,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Number of total transcripts for sample in each class"))
 plot(ggplot(tot_summary,aes(x=class,y=total_transcripts,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Number of total transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90)))
	
#boxplot(as.numeric(summary[names(Single_Cris),"total_transcripts"]) ~ Single_Cris ,main="Number of total transcript per sample",xlab="CRIS classes",ylab="total transcript number",col=c("red","darkgreen","lightblue","orange","blue"))

        Single_Cris_coverage <- merge(Single_Cris,coverage,"sample")
        Primary_Cris_coverage <- merge(Primary_Cris,coverage,"sample")
        Secondary_Cris_coverage <- merge(Secondary_Cris,coverage,"sample")
        tot_coverage <- rbind(Primary_Cris_coverage,Secondary_Cris_coverage)
 plot(ggplot(Single_Cris_coverage,aes(x=class,y=coverage,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Number of total reads for sample in each class"))
 plot(ggplot(tot_coverage,aes(x=class,y=coverage,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Number of total reads for sample in each class")+theme(axis.text.x = element_text(angle = 90)))

#boxplot(as.numeric(coverage[names(Single_Cris),2]) ~ Single_Cris ,main="Number of uniquely\nassigned reads per sample",xlab="CRIS classes",ylab="reads number",col=c("red","darkgreen","lightblue","orange","blue"))

 plot(ggplot(Single_Cris_summary,aes(x=class,y=total_annotated,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Number of total annotated transcripts for sample in each class"))
 plot(ggplot(tot_summary,aes(x=class,y=total_annotated,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Number of total annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90)))


#boxplot(as.numeric(summary[names(Single_Cris),"total_annotated"]) ~ Single_Cris ,main="Number of annotated transcript per sample",xlab="CRIS classes",ylab="total annotatedtranscript number",col=c("red","darkgreen","lightblue","orange","blue"))

 plot(ggplot(Single_Cris_summary,aes(x=class,y=total_unannotated,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Number of total not annotated transcripts for sample in each class"))
 plot(ggplot(tot_summary,aes(x=class,y=total_unannotated,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Number of total not annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90)))


#boxplot(as.numeric(summary[names(Single_Cris),"total_unannotated"]) ~ Single_Cris ,main="Number of  not annotated transcript per sample",xlab="CRIS classes",ylab="total not annotatedtranscript number",col=c("red","darkgreen","lightblue","orange","blue"))

 plot(ggplot(Single_Cris_summary,aes(x=class,y=assigned_annotated,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Number of assigned annotated transcripts for sample in each class"))
 plot(ggplot(tot_summary,aes(x=class,y=assigned_annotated,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Number of assigned annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90)))

#boxplot(as.numeric(summary[names(Single_Cris),"assigned_annotated"]) ~ Single_Cris ,main="Number of assigned annotated transcript per sample",xlab="CRIS classes",ylab="total transcript number",col=c("red","darkgreen","lightblue","orange","blue"))

 plot(ggplot(Single_Cris_summary,aes(x=class,y=assigned_unannotated,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Number of assigned not annotated transcripts for sample in each class"))
 plot(ggplot(tot_summary,aes(x=class,y=assigned_unannotated,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Number of assigned not annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90)))


#boxplot(as.numeric(summary[names(Single_Cris),"assigned_unannotated"]) ~ Single_Cris ,main="Number of assigned unannotated transcript per sample",xlab="CRIS classes",ylab="total transcript number",col=c("red","darkgreen","lightblue","orange","blue"))

		 plot(ggplot(Single_Cris_summary,aes(x=class,y=(total_transcripts-assigned_transcript),fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Number of unassigned transcripts for sample in each class"))
 plot(ggplot(tot_summary,aes(x=class,y=(total_transcripts-assigned_transcript),fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Number of unassigned transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90)))


#boxplot(as.numeric(summary[names(Single_Cris),"total_transcripts"])-as.numeric(summary[names(Single_Cris),"assigned_transcript"]) ~ Single_Cris ,main="Number of unassigned transcript per sample",xlab="CRIS classes",ylab="total transcript number",col=c("red","darkgreen","lightblue","orange","blue"))

		 plot(ggplot(Single_Cris_summary,aes(x=class,y=((assigned_unannotated*100)/assigned_transcript),fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("% of assigned not annotated transcripts for sample in each class"))
 plot(ggplot(tot_summary,aes(x=class,y=((assigned_unannotated*100)/assigned_transcript),fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("% of assigned not annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90)))


#boxplot(as.numeric(summary[names(Single_Cris),"assigned_unannotated"])*100/as.numeric(summary[names(Single_Cris),"assigned_transcript"]) ~ Single_Cris ,main="% of assigned unannotated transcripts \nover total assigned transcripts",xlab="CRIS classes",ylab="transcript %",col=c("red","darkgreen","lightblue","orange","blue"))

		 plot(ggplot(Single_Cris_summary,aes(x=class,y=((assigned_unannotated*100)/total_unannotated),fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Rate of assignment for not annotated transcripts"))
 plot(ggplot(tot_summary,aes(x=class,y=((assigned_unannotated*100)/total_unannotated),fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Rate of assignment for not annotated transcripts")+theme(axis.text.x = element_text(angle = 90)))


#boxplot(as.numeric(summary[names(Single_Cris),"assigned_unannotated"])*100/as.numeric(summary[names(Single_Cris),"total_unannotated"]) ~ Single_Cris ,main="Rate of assignment for unannotated transcripts",xlab="CRIS classes",ylab="transcript %",col=c("red","darkgreen","lightblue","orange","blue"))

		 plot(ggplot(Single_Cris_summary,aes(x=class,y=((unannotated_IR_containing_transcripts*100)/total_unannotated),fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("% IR containing not annotated transcripts for sample in each class"))
 plot(ggplot(tot_summary,aes(x=class,y=((unannotated_IR_containing_transcripts*100)/total_unannotated),fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("% of IR containing not annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90)))


#boxplot(as.numeric(summary[names(Single_Cris),"unannotated_IR_containing_transcripts"])*100/as.numeric(summary[names(Single_Cris),"total_unannotated"]) ~ Single_Cris ,main="% of IR containing transcripts within unannotated transcripts",xlab="CRIS classes",ylab="transcript %",col=c("red","darkgreen","lightblue","orange","blue"))
	
	#vioplot(as.numeric(summary[names(Single_Cris),"total_transcripts"]) ~ Single_Cris ,main="Number of total transcript per sample",xlab="CRIS classes",ylab="total transcript number",col=c("red","darkgreen","lightblue","orange","blue"))
	
	
plot(log10(as.numeric(summary$total_transcripts)),log10(coverage[summary$sample,"coverage"]),main=paste("Pearson correlation between transcript number\n and assigned reads:",cor(as.numeric(summary$total_transcripts),coverage[summary$sample,"coverage"],method="pearson"),sep=" "),xlab="total transcripts per sample",ylab="assigned reads",pch=20,col=ifelse(summary$sample%in%names(cla)==TRUE,"#00AFBB","black"))

	A <- vector()
	for(i in Single_Cris$sample[which(Single_Cris$class=="onlyCRISA")]){A<-c(A,unname(total_isoforms[,i]))}
	dist_A <- table(A)*100/(sum(table(A)))
	B <- vector()
	for(i in Single_Cris$sample[which(Single_Cris$class=="onlyCRISB")]){B<-c(B,unname(total_isoforms[,i]))}
	dist_B <- table(B)*100/(sum(table(B)))
	C <- vector()
	for(i in Single_Cris$sample[which(Single_Cris$class=="onlyCRISC")]){C<-c(C,unname(total_isoforms[,i]))}
	dist_C <- table(C)*100/(sum(table(C)))
	D <- vector()
	for(i in Single_Cris$sample[which(Single_Cris$class=="onlyCRISD")]){D<-c(D,unname(total_isoforms[,i]))}
	dist_D <- table(D)*100/(sum(table(D)))
	E <- vector()
	for(i in Single_Cris$sample[which(Single_Cris$class=="onlyCRISE")]){E<-c(E,unname(total_isoforms[,i]))}
	dist_E <- table(E)*100/(sum(table(E)))
	plot(as.numeric(names(table(A))),as.numeric(dist_A),main="Distribution of number of isoforms per gene per sample",col="red",type="l",
xlim=c(0,max(c(as.numeric(names(table(A)),names(table(B)),names(table(C)),names(table(D)),names(table(E)))))),xlab="Number of isoforms per gene", ylab="%") 
	lines(as.numeric(names(table(B))),as.numeric(dist_B),col="darkgreen")
	lines(as.numeric(names(table(C))),as.numeric(dist_C),col="lightblue")
	lines(as.numeric(names(table(D))),as.numeric(dist_D),col="yellow")
	lines(as.numeric(names(table(E))),as.numeric(dist_E),col="blue")


	dev.off()
