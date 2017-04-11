
rm(list=ls())
#setwd('/home/zhu/rushdata/tacc')
setwd('/work/03963/zhu2/rushdata_qtl/')

args = as.numeric(commandArgs(trailingOnly=TRUE))

if(length(args)==0){args <- 1}
fn <- paste0('cut/methylation_fpca_',args,'.rda')
load(fn)
load('expr4cut.rda')

#CCA

p_ginv_sq <- function(X,p){
  X.eigen = eigen(X);
  X.rank = sum(X.eigen$values>1e-8);
  X.value = X.eigen$values[1:X.rank]^(-1*p);
  if (length(X.value)==1){
    D = as.matrix(X.value);
  }else{
    D = diag(X.value);
  }
  rlt = X.eigen$vectors[,1:X.rank] %*% D %*% t(X.eigen$vectors[,1:X.rank]);
  return(rlt);
}
mrank <- function(X){
  X.svd = svd(X);
  X.rank = sum(X.svd$d>1e-6);
  return(X.rank);
}
mrank_sq <- function(X){
  X.eigen = eigen(X);
  X.rank = sum(Re(X.eigen$values)>1e-6);
  return(X.rank);
}
CCA_chisq_test <- function(rho,n,p,q){
  tstat = -1*n*sum(log(1-rho^2));
  p_value = pchisq(tstat,(p*q),lower.tail=FALSE);
  return(p_value);          
}
cca <- function(A,B){
  n = nrow(A);
  p = mrank(A);
  q = mrank(B);
  if (p <= q){
    X = A;
    Y = B;
  }else{
    X = B;
    Y = A;
  }
  R = p_ginv_sq(cov(X),0.5) %*% cov(X,Y) %*% p_ginv_sq(cov(Y),1) %*% cov(Y,X) %*% p_ginv_sq(cov(X),0.5);
  k = mrank_sq(R);
  d = Re(eigen(R)$values);
  rho = d[1:k]^(0.5);
  rho[rho >= 0.9999]=0.9;
  chisq_p = CCA_chisq_test(rho,n,p,q);
  return(c("chisq_p"=chisq_p,"df"=p*q));
}
ccap <- function(A,B){as.numeric(cca(A,B)[1])}
ccap2 <- function(A,B){
	out <- try(ccap(A,B))
	if(is.numeric(out)){
		return(out)
	} else {
		return(NA)
	}}

#Modeling
#i <- 0
time1 <- Sys.time()
system.time(rlt <- sapply(x,function(x){
	#print(Sys.time())
	#print(i<<-i+1)
	sapply(exprs,function(y){
		ccap2(x,y)
		})
	}))
time2 <- Sys.time()
timeu <- time2-time1

rlt <- list(rlt,args,timeu)
save(rlt,file=paste0('rlt/mfpca2expr_',args,'.rda'))
