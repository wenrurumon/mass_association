
rm(list=ls())

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

