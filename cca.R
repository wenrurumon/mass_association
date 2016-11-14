
rm(list=ls())

##########################
# Macro
##########################

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
	return(c("rho"=rho,"chisq_p"=chisq_p,"df"=p*q));
}
ccap <- function(A,B){as.numeric(cca(A,B)[2])}
ccaps <- function(As,Bs){
	rlt <- lapply(As,function(A){
		sapply(Bs,function(B){
			ccap(A,B)
		})
	})
	return(unlist(rlt))
}

##########################
# Load data
##########################

setwd('/home/zhu/rushdata/processed')

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){args <- c(20,20)}
i <- as.numeric(args[[1]])
j <- as.numeric(args[[2]])

load(paste0('methylation/processed_methy_chr',i,'.rda'))
load(paste0('snp/processed_snps_chr',j,'.rda'))
loadscore <- function(x){x[[1]][,1:which(x[[2]]>0.9)[1],drop=F]}

methyid <- rownames(pmethy[[1]][[1]]); snpid <- as.numeric(gsub('[^[:digit:]]','',rownames(psnp[[1]][[1]])))
pid <- intersect(methyid,snpid)
methyid <- match(pid,methyid); snpid <- match(pid,snpid)

##########################
# Association
##########################

cca_methyi_snpj <- function(i,j){
	methyi <- lapply(pmethy[[i]][-1:-2],function(x){
		x[[1]] <- x[[1]][methyid,,drop=F]
		loadscore(x)
	})
	snpj <- lapply(psnp[[j]][-1:-2],function(x){
		x[[1]] <- x[[1]][snpid,,drop=F]
		loadscore(x)
	})
	if((length(methyi)*length(snpj))==0){
		return(NULL)
	} else {
		return(c(methyi=i,snpj=j,ccaps(methyi,snpj)))
	}
}
system.time(rlt <- lapply(1:length(pmethy),function(i){
	do.call(rbind,
		lapply(1:length(psnp),function(j){
			cca_methyi_snpj(i,j)
		}))
}))
rlt <- do.call(rbind,rlt))

save(rlt,file=paste0('/home/zhu/rushdata/qtl_methylation_snp/cca_methy_snp_',i,'_',j,'.rda'))
