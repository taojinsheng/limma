#!/usr/bin/env Rscript
# CopyRight AnchorDx,All rights reserved
# CopyRight Jinsheng Tao <jinsheng_tao@anchordx.com>
args <- commandArgs(T)

comb=function(x){
	y='';
	for (i in 1:length(x)){
		for (j in 1:length(x)){
			if (i == j) next;
			y=paste(x[i],"-",x[j]," ",y,sep="")
		}
	}
	y
}

exprSet=read.delim(args[1],row.names=1, comment.char ="!",stringsAsFactors=F,header=T,check.names=F)

groupdata <- read.table(args[2],header=TRUE,stringsAsFactors=F)
design <- model.matrix(~0+factor(groupdata$group))

colnames(design)=levels(factor(groupdata$group))
rownames(design)=groupdata$geo_accession

library("limma")

fit <- lmFit(exprSet,design)

# add gene annotation to topTable output
library(annotate)
#fit$genes$Symbol <- getSYMBOL(fit$genes$ID,"hgu133plus2.db")
#library(hgu133plus2.db)
#fit$genes$Symbol <- getSYMBOL(rownames(exprSet),"hgu133plus2.db")
fit$genes$Symbol <- rownames(exprSet)

contrastsCommand=unlist(strsplit(comb(unique(groupdata$group)), split=" "))
cont.matrix <- makeContrasts(contrasts=contrastsCommand, levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#TODO : Add study name as acolumn on ouput
for(i in 1:length(contrastsCommand)){
	tempOutFile <- paste("diffexp.", contrastsCommand[i],".txt", sep="");#tempOutFile <- file.path(dir_meta,tempOutFile)
	tempOutput = topTable(fit2, coef=i, n=Inf)
	write.table(tempOutput,tempOutFile,sep="\t",quote=FALSE,row.names=F)
}

cat("done success !\n")
