#### Library
library(Biobase)
library(ROCR)
library(ggplot2)

#### Load dataset
load("out/predictive_genes.rda")
load("out/predictive_isoform_proportions.rda")
load("out/signature_genes.rda")

#### Analysis
#ROC on genes
aucs_genes_A <- vector()
aucs_genes_B <- vector()
aucs_genes_C <- vector()
aucs_genes_D <- vector()
aucs_genes_E <- vector()
for (i in 1:length(fData(predictive_genes)$class_specificity))
{
  cris_class <- fData(predictive_genes)$class_specificity[i]
  event <- rep(0,length(pData(predictive_genes)$CRIS_class))
  event[which(pData(predictive_genes)$CRIS_class==cris_class)] <- 1
  vals <- exprs(predictive_genes)[i,]
  pred <- prediction(vals, event)
  perf <- performance(pred,"tpr", "fpr")
  au <- performance(pred,"auc")
  #pdf(file=paste("out/ROC_",names(Scores)[[i]],"times",tempo, "months_new.pdf", sep=""), width=10, heigh=10 )
  #par(mar=c(5,7,4,2), cex.axis=2.5, las=1)
  #plot(perf, lwd=3.5, main=rownames(pData_predictive_genes)[i], col="blue", 
  #     xlab="", ylab="", cex.main=3)
  #par(cex.axis=2.5, las=0)
  #mtext("True positive rate", side=2, line=4.5, cex=2.5, font = 2)
  #mtext("False positive rate", side=1, line=3, cex=2.5, font = 2)
  #abline(a=0, b=1, col="black", lwd=1)
  #legend("bottomright", legend=c(paste("AUC = ",round(au@y.values[[1]], 
  #                                                    digits=3))),
  #       text.col=c("black"), cex=1.5)
  #dev.off()
  if(fData(predictive_genes)$class_specificity[i]=="onlyCRISA")
  {aucs_genes_A <- c(aucs_genes_A, au@y.values[[1]])}
  else if(fData(predictive_genes)$class_specificity[i]=="onlyCRISB")
  {aucs_genes_B <- c(aucs_genes_B, au@y.values[[1]])}
  else if(fData(predictive_genes)$class_specificity[i]=="onlyCRISC")
  {aucs_genes_C <- c(aucs_genes_C, au@y.values[[1]])}
  else if(fData(predictive_genes)$class_specificity[i]=="onlyCRISD")
  {aucs_genes_D <- c(aucs_genes_D, au@y.values[[1]])}
  else if(fData(predictive_genes)$class_specificity[i]=="onlyCRISE")
  {aucs_genes_E <- c(aucs_genes_E, au@y.values[[1]])}
}


aucs_iso_A <- vector()
aucs_iso_B <- vector()
aucs_iso_C <- vector()
aucs_iso_D <- vector()
aucs_iso_E <- vector()
for (i in 1:length(fData(predictive_isoform_proportions)$class_specificity))
{
  cris_class <- fData(predictive_isoform_proportions)$class_specificity[i]
  event <- rep(0,length(pData(predictive_isoform_proportions)$CRIS_class))
  event[which(pData(predictive_isoform_proportions)$CRIS_class==cris_class)] <- 1
  vals <- exprs(predictive_isoform_proportions)[i,]
  pred <- prediction(vals, event)
  perf <- performance(pred,"tpr", "fpr")
  au <- performance(pred,"auc")
  if(fData(predictive_isoform_proportions)$class_specificity[i]=="onlyCRISA")
  {aucs_iso_A <- c(aucs_iso_A, au@y.values[[1]])}
  else if(fData(predictive_isoform_proportions)$class_specificity[i]=="onlyCRISB")
  {aucs_iso_B <- c(aucs_iso_B, au@y.values[[1]])}
  else if(fData(predictive_isoform_proportions)$class_specificity[i]=="onlyCRISC")
  {aucs_iso_C <- c(aucs_iso_C, au@y.values[[1]])}
  else if(fData(predictive_isoform_proportions)$class_specificity[i]=="onlyCRISD")
  {aucs_iso_D <- c(aucs_iso_D, au@y.values[[1]])}
  else if(fData(predictive_isoform_proportions)$class_specificity[i]=="onlyCRISE")
  {aucs_iso_E <- c(aucs_iso_E, au@y.values[[1]])}
}

#Desnity plot
aucs_A <- cbind(c(aucs_genes_A, aucs_iso_A), 
                c(rep("genes", length(aucs_genes_A)), 
                  rep("isoform_proportions", length(aucs_iso_A))))
aucs_A <- data.frame(aucs_A)
colnames(aucs_A) <- c("AUCs", "Predictor")
aucs_A$AUCs <- as.numeric(aucs_A$AUCs)
ggplot_aucs_A <- ggplot(aucs_A, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_B <- cbind(c(aucs_genes_B, aucs_iso_B), 
                c(rep("genes", length(aucs_genes_B)),
                  rep("isoform_proportions", length(aucs_iso_B))))
aucs_B <- data.frame(aucs_B)
colnames(aucs_B) <- c("AUCs", "Predictor")
aucs_B$AUCs <- as.numeric(aucs_B$AUCs)
ggplot_aucs_B <- ggplot(aucs_B, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_C <- cbind(c(aucs_genes_C, aucs_iso_C), 
                c(rep("genes", length(aucs_genes_C)), 
                  rep("isoform_proportions", length(aucs_iso_C))))
aucs_C <- data.frame(aucs_C)
colnames(aucs_C) <- c("AUCs", "Predictor")
aucs_C$AUCs <- as.numeric(aucs_C$AUCs)
ggplot_aucs_C <- ggplot(aucs_C, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_D <- cbind(c(aucs_genes_D, aucs_iso_D), 
                c(rep("genes", length(aucs_genes_D)), 
                  rep("isoform_proportions", length(aucs_iso_D))))
aucs_D <- data.frame(aucs_D)
colnames(aucs_D) <- c("AUCs", "Predictor")
aucs_D$AUCs <- as.numeric(aucs_D$AUCs)
ggplot_aucs_D <- ggplot(aucs_D, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_E <- cbind(c(aucs_genes_E, aucs_iso_E), 
                c(rep("genes", length(aucs_genes_E)), 
                  rep("isoform_proportions", length(aucs_iso_E))))
aucs_E <- data.frame(aucs_E)
colnames(aucs_E) <- c("AUCs", "Predictor")
aucs_E$AUCs <- as.numeric(aucs_E$AUCs)
ggplot_aucs_E <- ggplot(aucs_E, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

pdf("out/density_plot_genes_iso.pdf", width=10, heigh=10)
print(ggplot_aucs_A)
print(ggplot_aucs_B)
print(ggplot_aucs_C)
print(ggplot_aucs_D)
print(ggplot_aucs_E)
dev.off()





### comparing signature genes
aucs_signgenes_A <- vector()
aucs_signgenes_B <- vector()
aucs_signgenes_C <- vector()
aucs_signgenes_D <- vector()
aucs_signgenes_E <- vector()
for (i in 1:length(fData(signature_genes)$class_specificity))
{
  cris_class <- fData(signature_genes)$class_specificity[i]
  event <- rep(0,length(pData(signature_genes)$CRIS_class))
  event[which(pData(signature_genes)$CRIS_class==cris_class)] <- 1
  vals <- exprs(signature_genes)[i,]
  pred <- prediction(vals, event)
  perf <- performance(pred,"tpr", "fpr")
  au <- performance(pred,"auc")
  if(fData(signature_genes)$class_specificity[i]=="onlyCRISA")
  {aucs_signgenes_A <- c(aucs_signgenes_A, au@y.values[[1]])}
  else if(fData(signature_genes)$class_specificity[i]=="onlyCRISB")
  {aucs_signgenes_B <- c(aucs_signgenes_B, au@y.values[[1]])}
  else if(fData(signature_genes)$class_specificity[i]=="onlyCRISC")
  {aucs_signgenes_C <- c(aucs_signgenes_C, au@y.values[[1]])}
  else if(fData(signature_genes)$class_specificity[i]=="onlyCRISD")
  {aucs_signgenes_D <- c(aucs_signgenes_D, au@y.values[[1]])}
  else if(fData(signature_genes)$class_specificity[i]=="onlyCRISE")
  {aucs_signgenes_E <- c(aucs_signgenes_E, au@y.values[[1]])}
}


#Desnity plot
aucs_A <- cbind(c(aucs_signgenes_A, aucs_iso_A), c(rep("sign_genes", 
                                                       length(aucs_signgenes_A)), 
                                                   rep("isoform_proportions", 
                                                       length(aucs_iso_A))))
aucs_A <- data.frame(aucs_A)
colnames(aucs_A) <- c("AUCs", "Predictor")
aucs_A$AUCs <- as.numeric(aucs_A$AUCs)
ggplot_aucs_A <- ggplot(aucs_A, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_B <- cbind(c(aucs_signgenes_B, aucs_iso_B), 
                c(rep("sign_genes", length(aucs_signgenes_B)),
                  rep("isoform_proportions", length(aucs_iso_B))))
aucs_B <- data.frame(aucs_B)
colnames(aucs_B) <- c("AUCs", "Predictor")
aucs_B$AUCs <- as.numeric(aucs_B$AUCs)
ggplot_aucs_B <- ggplot(aucs_B, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_C <- cbind(c(aucs_signgenes_C, aucs_iso_C), 
                c(rep("sign_genes", length(aucs_signgenes_C)), 
                  rep("isoform_proportions", length(aucs_iso_C))))
aucs_C <- data.frame(aucs_C)
colnames(aucs_C) <- c("AUCs", "Predictor")
aucs_C$AUCs <- as.numeric(aucs_C$AUCs)
ggplot_aucs_C <- ggplot(aucs_C, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_D <- cbind(c(aucs_signgenes_D, aucs_iso_D), 
                c(rep("sign_genes", length(aucs_signgenes_D)), 
                  rep("isoform_proportions", length(aucs_iso_D))))
aucs_D <- data.frame(aucs_D)
colnames(aucs_D) <- c("AUCs", "Predictor")
aucs_D$AUCs <- as.numeric(aucs_D$AUCs)
ggplot_aucs_D <- ggplot(aucs_D, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_E <- cbind(c(aucs_signgenes_E, aucs_iso_E), 
                c(rep("sign_genes", length(aucs_signgenes_E)), 
                  rep("isoform_proportions", length(aucs_iso_E))))
aucs_E <- data.frame(aucs_E)
colnames(aucs_E) <- c("AUCs", "Predictor")
aucs_E$AUCs <- as.numeric(aucs_E$AUCs)
ggplot_aucs_E <- ggplot(aucs_E, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

pdf("out/density_plot_signatures_iso.pdf", width=10, heigh=10)
print(ggplot_aucs_A)
print(ggplot_aucs_B)
print(ggplot_aucs_C)
print(ggplot_aucs_D)
print(ggplot_aucs_E)
dev.off()



#Desnity plot sign genes or genes cinzia
aucs_A <- cbind(c(aucs_signgenes_A, aucs_genes_A), c(rep("sign_genes", 
                                                       length(aucs_signgenes_A)), 
                                                   rep("genes", 
                                                       length(aucs_genes_A))))
aucs_A <- data.frame(aucs_A)
colnames(aucs_A) <- c("AUCs", "Predictor")
aucs_A$AUCs <- as.numeric(aucs_A$AUCs)
ggplot_aucs_A <- ggplot(aucs_A, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_B <- cbind(c(aucs_signgenes_B, aucs_genes_B), 
                c(rep("sign_genes", length(aucs_signgenes_B)),
                  rep("genes", length(aucs_genes_B))))
aucs_B <- data.frame(aucs_B)
colnames(aucs_B) <- c("AUCs", "Predictor")
aucs_B$AUCs <- as.numeric(aucs_B$AUCs)
ggplot_aucs_B <- ggplot(aucs_B, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_C <- cbind(c(aucs_signgenes_C, aucs_genes_C), 
                c(rep("sign_genes", length(aucs_signgenes_C)), 
                  rep("genes", length(aucs_genes_C))))
aucs_C <- data.frame(aucs_C)
colnames(aucs_C) <- c("AUCs", "Predictor")
aucs_C$AUCs <- as.numeric(aucs_C$AUCs)
ggplot_aucs_C <- ggplot(aucs_C, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_D <- cbind(c(aucs_signgenes_D, aucs_genes_D), 
                c(rep("sign_genes", length(aucs_signgenes_D)), 
                  rep("genes", length(aucs_genes_D))))
aucs_D <- data.frame(aucs_D)
colnames(aucs_D) <- c("AUCs", "Predictor")
aucs_D$AUCs <- as.numeric(aucs_D$AUCs)
ggplot_aucs_D <- ggplot(aucs_D, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

aucs_E <- cbind(c(aucs_signgenes_E, aucs_genes_E), 
                c(rep("sign_genes", length(aucs_signgenes_E)), 
                  rep("genes", length(aucs_genes_E))))
aucs_E <- data.frame(aucs_E)
colnames(aucs_E) <- c("AUCs", "Predictor")
aucs_E$AUCs <- as.numeric(aucs_E$AUCs)
ggplot_aucs_E <- ggplot(aucs_E, aes(x=AUCs, fill=Predictor), trim=FALSE) + 
  geom_density(alpha=.3, trim=FALSE) +
  theme(panel.grid.major = element_line(size = .7, color = "white"), 
        axis.line = element_line(size=.7, color = "black"),panel.grid.minor = element_blank(), 
        text = element_text(size=25)) + ylab("Density") + xlab("AUCs")

pdf("out/density_plot_signatures_genes.pdf", width=10, heigh=10)
print(ggplot_aucs_A)
print(ggplot_aucs_B)
print(ggplot_aucs_C)
print(ggplot_aucs_D)
print(ggplot_aucs_E)
dev.off()
