###library&loading
	class <- read.table("data/CRIS5/cell_bank_normTMMF_prediction_result.xls",header=T, row.names=1)

###code
	class_assign <- t(apply(class[,c(16:20)],1,function(x){ifelse(as.numeric(x)<0.2,1,0)}))
        colnames(class_assign) <- colnames(class)[c(16:20)]
        su <- apply(class_assign,1,sum)
        cla <- rep("NA",length(su))
        cla[which(su==5)] <- "all"
        cla[which(su==1)] <-apply(class_assign[which(su==1),],1,function(x){paste(c("onlyCRIS",strsplit(names(x)[which(x==1)],"[.]")[[1]][4]),collapse="")})
        cla[which(su>1)] <-apply(class[which(su>1),c(16:20)],1,function(x){paste(c("primaryCRIS",strsplit(names(x)[which(x==min(x))],"[.]")[[1]][4]),collapse="")})
        secondary_class <- class[which(su>1),c(16:20)]
        colnames(secondary_class) <- unlist(lapply(colnames(secondary_class),function(x){paste(c("CRIS",strsplit(x,"[.]")[[1]][4]),collapse="")}))
        sec <- apply(secondary_class,1,function(x){paste(names(x)[which(x<0.2&x>min(x))],collapse="_")})
        cla[which(su>1)] <- paste(cla[which(su>1)],sec,sep="secondary")
	names(cla) <- paste(rownames(class),".STAR",sep="")
        classes <- ifelse(cla=="onlyCRISD","OnlyD","NotD")
        classes[grep("primaryCRISD",cla)] <- "PrimaryD"
        classes[setdiff(grep("CRISD",cla),which(classes!="NotD"))] <- "SecondaryD"
        classes <- factor(classes,levels=c("OnlyD","PrimaryD","SecondaryD","NotD"))
###
	save(cla, file="out/CRIS_class.rda")
	save(classes, file="out/CRIS_D_classification.rda")
