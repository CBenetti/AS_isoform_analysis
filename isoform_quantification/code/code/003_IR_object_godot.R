###Matrix of IRFinder results
        out <-  system("ls out/IRFinder/out/",intern=TRUE) 
        out     <- out[grep(".STAR",out)]
	setwd("out/IRFinder/out")
        for(i in 1:length(out)){ 
                tab <- read.table(paste(out[i],"/IRFinder-IR-nondir-val.txt",sep=""),header=TRUE)
                nam_tab <- paste(tab[,1],tab[,2],tab[,3],tab[,4],tab[,6],sep=":")
		nam_tab <- gsub("/",":",nam_tab)
                rownames(tab) <- nam_tab
                if(i==1){nam<- nam_tab}else{nam <- unique(c(nam,nam_tab))}
                save(tab, file=paste(out[i],"/tab.rda",sep=""))
        }
        IR_ratio <- matrix(0,ncol=length(out),nrow=length(nam))
        IR_score <- matrix(0,ncol=length(out),nrow=length(nam))
        rownames(IR_ratio) <- nam
        rownames(IR_score) <- nam
        for(i in 1:length(out)){
                load(paste(out[i],"/tab.rda",sep=""))
                IR_ratio[rownames(tab),i] <- tab$IRratio
                IR_score[rownames(tab),i] <- tab$CNN_IRscore
        }
        #colnames(IR_ratio) <- unlist(lapply(out,function(x){strsplit(x,"/")[[1]][3]}))
        #colnames(IR_score) <- unlist(lapply(out,function(x){strsplit(x,"/")[[1]][3]}))
        colnames(IR_ratio) <- out
	colnames(IR_score) <- out
	setwd("../")
	save(IR_score,file="IR_score.rda")
        save(IR_ratio,file="IR_ratio.rda")
        write.table(IR_score,"IR_score.txt",sep="/t",col.names=NA)
        write.table(IR_ratio,"IR_ratio.txt",sep="/t",col.names=NA)
	setwd("../../")
