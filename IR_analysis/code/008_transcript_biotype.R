###library
	library(vioplot)
	library(ggplot2)
###loading
	load("data/transcripts_info.rda")
	load("data/transcripts_biotype_tpm.rda")
	load("out/CRIS_class.rda")
###code
        Single_Cris <- data.frame(class=cla[grep("onlyCRIS",cla)],sample=names(cla)[grep("onlyCRIS",cla)])
        Primary<- sapply(strsplit(cla,"secondary"),"[[",1)
        Primary_Cris <- data.frame(class=Primary,sample=names(Primary))
        Secondary<-lapply(c("CRISA","CRISB","CRISC","CRISD","CRISE"),function(x){
                sec <- data.frame(names(cla)[setdiff(grep(x,cla),grep(x,Primary_Cris))],rep(paste("secondary",x,sep=""),length(setdiff(grep(x,cla),grep(x,Primary_Cris)))))
        })
        do.call(rbind,Secondary)-> Secondary_Cris
        colnames(Secondary_Cris) <- c("sample","class")
        Primary_Cris <- Primary_Cris[which(Primary_Cris$class!="NA"),]
	
	i <-data.frame(retained_intron=sapply(transcripts_info,"[[","retained_intron"),sample=names(sapply(transcripts_info,"[[","retained_intron")))
	i_tpm <- sapply(transcripts_biotype_tpm,"[[","retained_intron")
	pc <- sapply(transcripts_info,"[[","protein_coding")
	pc_tpm <- sapply(transcripts_biotype_tpm,"[[","protein_coding")
	tot <- unlist(lapply(transcripts_info,sum))
	perc <- data.frame(retained_intron=(i$retained_intron/tot)*100,sample=i$sample)
	perc_pc <- data.frame(protein_coding=(pc*100)/tot,sample=i$sample)
	perc_pc_Single <- merge(perc_pc,Single_Cris,"sample")
	perc_pc_Primary <- merge(perc_pc,Primary_Cris,"sample")
	perc_pc_Secondary <- merge(perc_pc,Secondary_Cris,"sample")
	perc_pc_tot <- rbind(perc_pc_Primary,perc_pc_Secondary)
	perc_Single <- merge(perc,Single_Cris,"sample")
	perc_Primary <- merge(perc,Primary_Cris,"sample")
	perc_Secondary <- merge(perc,Secondary_Cris,"sample")
	perc_tot <- rbind(perc_Primary,perc_Secondary)
	i_Single <- merge(i,Single_Cris,"sample")
	i_Primary <- merge(i,Primary_Cris,"sample")
	i_Secondary <- merge(i,Secondary_Cris,"sample")
	i_tot <- rbind(i_Primary,i_Secondary)
        pdf("out/transcript_biotype.pdf")

 plot(ggplot(perc_Single,aes(x=class,y=retained_intron,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Retained intron annotated transcripts for sample in each class")+ stat_compare_means(label = "p.signif", hide.ns = TRUE))
 plot(ggplot(perc_tot,aes(x=class,y=retained_intron,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Retained intron annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90))+stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", hide.ns = TRUE))

#boxplot(as.numeric(perc[names(Single_Cris)]) ~ Single_Cris ,main="Retained intron annotated transcripts per sample",xlab="CRIS classes",ylab="transcript %",col=c("red","darkgreen","lightblue","orange","blue"))

 plot(ggplot(perc_pc_Single,aes(x=class,y=protein_coding,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Protein coding annotated transcripts for sample in each class")+stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", hide.ns = TRUE))
 plot(ggplot(perc_pc_tot,aes(x=class,y=protein_coding,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Protein coding annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90))+ stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", hide.ns = TRUE))

#boxplot(as.numeric(perc_pc[names(Single_Cris)]) ~ Single_Cris ,main="Protein coding annotated transcripts per sample",xlab="CRIS classes",ylab="transcript %",col=c("red","darkgreen","lightblue","orange","blue"))
#vioplot(as.numeric(perc[names(Single_Cris)]) ~ Single_Cris ,main="Retained intron annotated transcripts per sample",xlab="CRIS classes",ylab="transcript %",col=c("red","darkgreen","lightblue","orange","blue"))

 plot(ggplot(i_Single,aes(x=class,y=retained_intron,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                scale_fill_manual(name = "CRIS",values=c("red","darkgreen","lightblue","orange","blue"))+
                ggtitle("Retained intron annotated transcripts for sample in each class")+ stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", hide.ns = TRUE))
 plot(ggplot(i_tot,aes(x=class,y=retained_intron,fill=class))+
                geom_boxplot(alpha=0.5)+geom_jitter(width=0.25,aes(color=class))+
                scale_color_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                scale_fill_manual(name = "CRIS",values=rep(c("red","darkgreen","lightblue","orange","blue"),3))+
                ggtitle("Retained intron annotated transcripts for sample in each class")+theme(axis.text.x = element_text(angle = 90))+ stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", hide.ns = TRUE))


#boxplot(as.numeric(i[names(Single_Cris)]) ~ Single_Cris ,main="Retained intron annotated transcripts per sample",xlab="CRIS classes",ylab="transcript number",col=c("red","darkgreen","lightblue","orange","blue"))
#vioplot(as.numeric(i[names(Single_Cris)]) ~ Single_Cris ,main="Retained intron annotated transcripts per sample",xlab="CRIS classes",ylab="transcript number",col=c("red","darkgreen","lightblue","orange","blue"))

#boxplot(as.numeric(i_tpm[names(Single_Cris)]) ~ Single_Cris ,main="Retained intron annotated transcripts per sample",xlab="CRIS classes",ylab="mean tpm",col=c("red","darkgreen","lightblue","orange","blue"))
#vioplot(as.numeric(i_tpm[names(Single_Cris)]) ~ Single_Cris ,main="Retained intron annotated transcripts per sample",xlab="CRIS classes",ylab="mean tpm",col=c("red","darkgreen","lightblue","orange","blue"))

	dev.off()

