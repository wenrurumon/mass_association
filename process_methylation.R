
rm(list=ls())

#######################
# Config
#######################

#Setup

args = as.numeric(commandArgs(trailingOnly=TRUE))
if(length(args)==0){args <- 22}#for test

#Import data

datafolder <- '/home/zhu/rushdata/rawdata/methylation_rawdata/methylation_imputed_all_chr/'
mapfolder <- '/home/zhu/rushdata/rawdata/methylation_rawdata/illqc_probes_loc/'

datafile <- paste0(datafolder,'dnaMeth_matrix_chr',args,'_748qc_imputed.txt')
posfile <- paste0(mapfolder,'illqc_probes_chr',args,'_loc.txt')
mapfile <- '/home/zhu/rushdata/rawdata/methylation_rawdata/GPL13534.map'
idfile <- '/home/zhu/rushdata/rawdata/methylation_rawdata/id_delt.txt'

raw <- t(read.table(datafile))
pos <- read.table(posfile,header=T)
map <- read.table(mapfile)
id <- read.table(idfile,header=T)
rownames(raw) <- id$projid

#######################
# Macro
#######################

library(fda)
library(MASS)
library(GenABEL)
library(flare)
library(corpcor)

#FPCA

scale_ivn <- function(x){apply(x,2,rntransform)}
myinv<-function(A){
  A_svd<-fast.svd(A)
  if(length(A_svd$d)==1){
    A_inv<-A_svd$v%*%as.matrix(1/A_svd$d)%*%t(A_svd$u)
  }else{
    A_inv<-A_svd$v%*%diag(1/A_svd$d)%*%t(A_svd$u)
  }
  return(A_inv)
}

fourier.expansion<- function(x,pos,nbasis,lambda){
  frange <- c(pos[1], pos[length(pos)])
  rlt=list();
  rlt$fbasis<-create.fourier.basis(frange,nbasis=nbasis)
  rlt$phi = eval.basis(pos,rlt$fbasis) + eval.basis(pos,rlt$fbasis,2)*lambda
  rlt$coef<-myinv(t(rlt$phi)%*%rlt$phi)%*%t(rlt$phi)%*%t(x)
  return(rlt)
}

fpca <- function(x,pos,nbasis,lambda){
	nbasis <- min(nrow(x),nbasis)
	#process data
	x <- x[,order(pos),drop=F]
	pos <- pos[order(pos)]
	pos <- (pos-min(pos))/(max(pos)-min(pos))
	#fourier expansion
	x.expanded <- fourier.expansion(x,pos,nbasis,lambda)
	fcoef<-scale(t(x.expanded$coef-rowMeans(x.expanded$coef))/sqrt(ncol(x)))
	#PCA
	A.svd <- try(svd(fcoef))
	while(!is.list(A.svd)){
		nbasis <- nbasis - 2
		x.expanded <- fourier.expansion(x,pos,nbasis,lambda)
		fcoef<-scale(t(x.expanded$coef-rowMeans(x.expanded$coef))/sqrt(ncol(x)))
		A.svd <- try(svd(fcoef))
	}
	prop1 <- (A.svd$d)^2; prop1 <- cumsum(prop1)/sum(prop1)
	r <- which(prop1>0.8)[1]
	d <- A.svd$d-A.svd$d[min(r+1,dim(fcoef))]
	d <- d[d>1e-10]
	prop2 <- d^2; prop2 <- cumsum(prop2)/sum(prop2)
	d <- diag(d,length(d),length(d))
	score1 <-  fcoef %*% A.svd$v
		score1 <- score1[,1:which(prop1>0.9999)[1],drop=F]
		prop1 <- prop1[1:which(prop1>0.9999)[1]]
	score2 <- A.svd$u[,1:ncol(d),drop=F] %*% sqrt(d)
	if(ncol(score1)==1){
		score3 <- score1
		prop3 <- prop1
	} else {
		score3 <- qpca(score1,rank=r)
		prop3 <- score3$prop
		score3 <- score3$X
	}
	list(fpca=list(score=score1,prop=prop1),
		qfpca=list(score=score2,prop=prop2),
		q2fpca=list(score=score3,prop=prop3))
}

qpca <- function(A,rank=0){
  A <- scale_ivn(A)
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-10]
  r <- length(d)
  prop <- d^2; prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop)
  return(rlt)
}

#Methylation Process

getmethy <- function(gene){
	print(paste(args,gene))
	xi <- raw[,map$V5==gene,drop=F]
	posi <- pos$MAPINFO[map$V5==gene];names(posi) <- pos$TargetID[map$V5==gene]
	nbasis <- length(posi)*2-1
	lambda <- 0
	sel <- length(posi)
	if( sum(sel) == 0 ){
		rlt <- NULL
	} else if (sum(sel)<=3){
		x.fpca <- list(
			fpca=list(score=xi,prop=c(rep(0,sum(sel)-1),1)),
			qfpca=list(score=xi,prop=c(rep(0,sum(sel)-1),1)),
			q2fpca=list(score=xi,prop=c(rep(0,sum(sel)-1),1))
		)
		rlt <- list(methy=xi,pos=posi)
		rlt[3:5] <- x.fpca; names(rlt)[3:5] <- names(x.fpca)
	} else {
		x.fpca <- fpca(x=xi,pos=posi,nbasis=nbasis,lambda=lambda)
		rlt <- list(methy=xi,pos=posi)
		rlt[3:5] <- x.fpca; names(rlt)[3:5] <- names(x.fpca)
	}
	return(rlt)
}

#######################
# Process
#######################

sel <- which(pos$TargetID%in%map$V4)
	raw <- raw[,sel,drop=F]
	pos <- pos[sel,,drop=F]
sel <- match(pos$TargetID,map$V4)
	map <- map[sel,]
genename <- paste(unique(map$V5))

pmethy <- lapply(genename,getmethy)
names(pmethy) <- genename
pmethy <- pmethy[!is.null(pmethy)]

print(filename <- paste0('/home/zhu/rushdata/processed/methylation/processed_methy_chr',args,'.rda'))
save(pmethy,file=filename)

##########################
# Batch Push
##########################

#R CMD BATCH --no-save --no-restore '--args 1' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 2' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 3' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 4' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 5' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 6' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 7' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 8' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 9' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 10' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 11' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 12' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 13' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 14' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 15' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 16' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 17' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 18' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 19' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 20' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 21' process_methylation.R &
#R CMD BATCH --no-save --no-restore '--args 22' process_methylation.R &
