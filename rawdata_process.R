
rm(list=ls())
#Function
hw_test <- function(X){
  n = length(X);
  n_11 = sum(X==0);
  n_12 = sum(X==1);
  n_22 = sum(X==2);
  p = (2*n_11+n_12)/(2*n);
  q = 1-p;
  t_stat = ((n_11-n*p^2)^2)/(n*(p^2)) + ((n_12-2*n*p*q)^2)/(2*n*p*q) + ((n_22-n*q^2)^2)/(n*(q^2));
  p_value = pchisq(t_stat,1,lower.tail=FALSE);
  return(p_value);
}
ccheck <- function(x){
  p2 <- mean(x==2)+1/2*mean(x==1)
  p2<0.95&p2>0.05
}
#PID
setwd('/Users/wenrurumon/Documents/uthealth/data')
load('expr_residual.rda')
pid <- rownames(expr_residual[[1]])
#getdata
setwd('/Users/wenrurumon/Documents/uthealth/snprawdata')

out <- lapply(1:22,function(i){
  print(i)
  load(paste0('chr',i,'.common'))
  load(paste0('chr',i,'.rare'))
  x1 <- lapply(common[[1]],function(x){
    x[match(pid,rownames(x)),,drop=F]
  })
  x2 <- lapply(rare[[1]],function(x){
    x[match(pid,rownames(x)),,drop=F]
  })
  gene <- unique(c(common[[2]],rare[[2]]))
  x <- lapply(gene,function(g){
    # print(g)
    g1 <- x1[which(names(x1)==g)]
    g1 <- ifelse(length(g1)==0,list(NULL),g1)
    g2 <- x2[which(names(x1)==g)]
    g2 <- ifelse(length(g2)==0,list(NULL),g2)
    if(is.null(g1[[1]])){return(g2[[1]])}
    if(is.null(g2[[1]])){return(g1[[1]])}
    return(cbind(g1[[1]],g2[[1]]))
  })
  gene <- gene[!sapply(x,is.null)]
  x <- x[!sapply(x,is.null)]
  x <- sapply(x,function(snp){
    sel1 <- apply(snp,2,var)!=0
    sel2 <- apply(snp,2,hw_test)>=1e-6;sel2[is.na(sel2)]<-FALSE
    sel3 <- colSums(snp)>3
    snp <- snp[,sel1&sel2&sel3,drop=F]
    c(common=sum(apply(snp,2,ccheck)),rare=sum(1-apply(snp,2,ccheck)))
  })
  x <- t(x)
  rownames(x) <- gene
  x <- cbind(chr=i,x)
  x
})
out <- do.call(rbind,out)
dim(out)
colSums(out)
