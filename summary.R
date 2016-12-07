rm(list=ls())

#Summary for snp result

load('/home/zhu/qtl/rawdata/rlt20161016/old/pathway_snp_cca.qtl')
rlt2p <- apply(rlt[,6:7],1,min)
rlt2p_sel <- rlt2p<(22*0.05/length(rlt2p))

setwd('/home/zhu/qtl/rawdata')
gene31 <- readLines('gene31.csv')
gene31_sel <- rlt$genenm%in%gene31
gene31_sel <- (gene31_sel & rlt2p<0.02)
rlt2 <- rlt[rlt2p_sel|gene31_sel,]
snptoget <- unique(rlt2[,c(2,4,5)])

setwd('/home/zhu/qtl/rawdata/1707snp_new')
fs <- paste0('chr',1:22,'.448fpca')
snpdata <- lapply(fs,function(x){
	print(x)
	load(x)
	return(rlt_fpca)
})
snp_sel <- lapply(1:nrow(snptoget),function(k){
	i <- as.numeric(snptoget[k,2])
	j <- as.numeric(snptoget[k,3])
	xj <- snpdata[[i]][[j]]
	xj1 <- xj[[3]][[1]][,1:which(xj[[3]][[2]]>0.9)[1],drop=F]#fpca
	xj1[260,] <- (xj1[260,]+xj1[261,])/2; xj1 <- xj1[-261,,drop=F]
	return(xj1)
})
names(snp_sel) <- snptoget[,1]
selected_snp <- list(
	scores=snp_sel,mapping=rlt2
)

#Summary for methylation result

setwd('/home/zhu/qtl/448samples/result20161021')
fs <- paste0('path_methy_',1:22,'.qtl')
methy_pheno_qtl <- lapply(fs,function(x){
	print(x)
	load(x)
	rlt2 <- t(do.call(cbind,rlt))
	rlt2p <- sapply(5:10,function(i){as.numeric(rlt2[,i])})
	rltsel <- apply(rlt2p,1,min)<(0.05/nrow(rlt2p))
	return(rlt2[rltsel,])
})
	
methy_pheno_qtl <- do.call(rbind,methy_pheno_qtl)
setwd('/home/zhu/qtl/448samples/processed')
load('methylation2.rda')

toget <- unique(methy_pheno_qtl[,c(1,3,4)])
methy_sel <- lapply(1:nrow(toget),function(k){
	i <- as.numeric(toget[k,1])
	j <- as.numeric(toget[k,2])
	methy[[i]][[j]]$qfpca
})
names(methy_sel) <- toget[,3]
selected_methy <- list(
	scores=methy_sel,mapping=methy_pheno_qtl
)
