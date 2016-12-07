rm(list=ls())

#Summary for methylation result

setwd('/home/zhu/qtl/448samples/result20161021')
fs <- paste0('path_methy_',1:22,'.qtl')
methy_pheno_qtl <- lapply(fs,function(x){
	load(x)
	rlt2 <- t(do.call(cbind,rlt))
	rlt2p <- sapply(5:10,function(i){as.numeric(rlt2[,i])})
	rltsel <- apply(rlt2p,1,min)<(0.05/nrow(rlt2p))
	return(rlt2[rltsel,])
})
	
