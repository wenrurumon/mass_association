
rm(list=ls())
library(qfcca)
library(GenABEL)
library(flare)
library(fda)
library(corpcor)
setwd("/Users/wenrurumon/Documents/uthealth/data")

load('rawdata4snpcandidate.rda')
snps <- do.call(c,rlt)
load('D6.rda')
snpid <- as.numeric(gsub('[^[:digit:]]','',rownames(snps[[1]][[1]])))
did <- as.numeric(rownames(D6))
pid <- intersect(did,snpid)
disease <- (D6[match(pid,did),,drop=F]+1)/2
snps <- lapply(snps,function(x){x$snp <- x$snp[match(pid,snpid),];return(x)})

#define a loss function to optimize the nbasis
loss <- function(x){cca(x,disease)[[1]]}

#Fpca macro for fpca and qfpca
fast_fpca <- function (x, pos=1:length(x), nbasis=length(pos)*2-1, lambda=0,prop=0.9){
  nbasis <- min(nrow(x), nbasis)
  x <- x[, order(pos), drop = F]
  pos <- pos[order(pos)]
  pos <- (pos - min(pos))/(max(pos) - min(pos))
  x.expanded <- fourier.expansion(x, pos, nbasis, lambda)
  fcoef <- scale(t(x.expanded$coef - rowMeans(x.expanded$coef))/sqrt(ncol(x)))
  A.svd <- try(svd(fcoef))
  while (!is.list(A.svd)) {
    nbasis <- nbasis - 2
    x.expanded <- fourier.expansion(x, pos, nbasis, lambda)
    fcoef <- scale(t(x.expanded$coef - rowMeans(x.expanded$coef))/sqrt(ncol(x)))
    A.svd <- try(svd(fcoef))
  }
  prop1 <- (A.svd$d)^2
  prop1 <- cumsum(prop1)/sum(prop1)
  r <- which(prop1 > prop)[1]
  d <- A.svd$d - A.svd$d[min(r + 1, dim(fcoef))]
  d <- d[d > 1e-10]
  prop2 <- d^2
  prop2 <- cumsum(prop2)/sum(prop2)
  d <- diag(d, length(d), length(d))
  score1 <- fcoef %*% A.svd$v
  score1 <- score1[, 1:which(prop1 > 0.9999)[1], drop = F]
  prop1 <- prop1[1:which(prop1 > 0.9999)[1]]
  score2 <- A.svd$u[, 1:ncol(d), drop = F] %*% sqrt(d)
  list(fpca=list(score=score1,prop=prop1),qfpca=list(score=score2,prop=prop2),nbasis=nbasis)
}
optim_fpca <- function(x,loss, pos=1:length(x), nbasis=length(pos)*2-1, lambda=0,prop=0.9,i=2, n = 10){n
  n <- floor(min(n,nbasis)/2)-1
  test <- lapply(0:n,function(j){
    nbasis <- nbasis - j * 2
    s <- fast_fpca(x,pos,nbasis,lambda,prop)
    l <- loss(s[[i]]$score[,1:which(s[[i]]$prop>=prop)[1],drop=F])
    list(x=s,loss=l)
  })
  return(test[[which(sapply(test,function(x){x$loss})==min(sapply(test,function(x){x$loss})))[1]]]$x)
}

#Main

i <- 0
test <- lapply(snps[1:10],function(x){
  print(i<<-i+1)
  optim_fpca(x$snp,loss,x$pos)
})


