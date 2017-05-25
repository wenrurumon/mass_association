
rm(list=ls())
library(reshape)
setwd('/home/zhu/rushdata/tacc/rlt')

#get snp2exor result
f <- dir(pattern='snp2expr')
rlt <- lapply(f,function(x){load(x);rlt[[1]]})
rlt.snp2expr <- do.call(cbind,rlt)

#get methylation2expr result
f <- dir(pattern='mfpca2expr')
rlt <- lapply(f,function(x){load(x);rlt[[1]]})
rlt.methy2expr <- do.call(cbind,rlt)

#get snp2methylation result
f <- dir(pattern='sfpca2mfpca_')
rlt <- lapply(f,function(x){load(x);
  do.call(cbind,rlts)
})
rlt.snp2methy <- do.call(rbind,rlt)

#summary

rlt <- list(snp2expr=rlt.snp2expr,methy2expr=rlt.methy2expr,snp2methy=rlt.snp2methy)
rlt2 <- lapply(rlt,function(x){
  print(1)
  x[is.na(x)] <- 1
  x <- (x*length(x))
  r <- lapply(1:ncol(x),function(i){
    y <- x[,i,drop=F]
    y <- y[y<=0.5,,drop=F]
    if(length(y)==0){
      return(NULL)
    }else{
      return(melt(y))
    }
  })
  do.call(rbind,r)
})
