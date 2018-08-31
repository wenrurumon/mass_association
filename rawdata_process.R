
rm(list=ls())

##################
# CCA
##################

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
ccap <- function(A,B){
  x <- try(cca(A,B))
  if(length(x)==1){NA}else(x[2])
}
model <- function(A,B,np=0){
  if(np==0){
    fun <- function(A,B,np){return(NA)}
  } else {
    fun <- function(A,B,np){permANM(A[,1,drop=F],B,np)$P_no_causation}
  }
  c(as.numeric(ccap(A,B)),
    as.numeric(ccap(A[,1,drop=F],B[,1,drop=F]))
    ,fun(A,B,np))
}
model_cca <- function(A,B,np=0){
  x <- try(as.numeric(cca(A,B)[2]))
  ifelse(is.numeric(x),x,NA)
}
models <- function(As,Bs,np=0){
  rlt <- lapply(As,function(A){
    sapply(Bs,function(B){
      model_cca(A,B,np=np)
    })
  })
  return((rlt))
}
source('/Users/wenrurumon/Documents/uthealth/data/contiANM_2.0.R')

##################
# loaddata
##################

load('/Users/wenrurumon/Documents/uthealth/data/nodes.rda')
cc <- cand <- readLines('/Users/wenrurumon/Documents/uthealth/data/candidates.txt')
cand <- unique(c(paste0('s:',cand),grep('s:',nodes$gg,value=T),grep('s:',nodes$pg,value=T)))
cand2 <- unique(c(paste0('m:',cand),grep('m:',nodes$gg,value=T),grep('m:',nodes$pg,value=T)))
setwd("/Users/wenrurumon/Documents/uthealth/data/cut_snp")
data.snp <- do.call(c,lapply(dir(),function(x){
  load(x)
  g <- names(x)
  g <- paste0('s',substr(g,regexpr(':',g)+1,nchar(g)))
  x[g%in%cand]
}))
setwd("/Users/wenrurumon/Documents/uthealth/data/cut_methy")
data.methy <- do.call(c,lapply(dir(),function(x){
  load(x)
  g <- names(x)
  g <- paste0('m',substr(g,regexpr(':',g)+1,nchar(g)))
  x[g%in%cand2]
}))
load('/Users/wenrurumon/Documents/uthealth/data/expr_residual.rda')
data.expr <- expr_residual[sapply(expr_residual,colnames)%in%nodes$gg]
names(data.expr) <- sapply(data.expr,colnames)
load('/Users/wenrurumon/Documents/uthealth/data/D6.rda')
data.dis <- lapply(1:ncol(D6),function(i){D6[,i,drop=F]})
# save(data.snp,data.expr,data.dis,data.methy,file='/Users/wenrurumon/Documents/uthealth/data/cutdata.rda')

scc <- paste0('s',substr(names(data.snp),regexpr(':',names(data.snp))+1,nchar(names(data.snp))))[(paste0('s',substr(names(data.snp),regexpr(':',names(data.snp))+1,nchar(names(data.snp))))%in%paste0('s:',cc))]
mcc <- paste0('m',substr(names(data.methy),regexpr(':',names(data.methy))+1,nchar(names(data.methy))))[(paste0('m',substr(names(data.methy),regexpr(':',names(data.methy))+1,nchar(names(data.methy))))%in%paste0('m:',cc))]


################
# GO ahead
################
#
# system.time(snp2expr500 <- models(data.snp[1:500],data.expr))
# system.time(snp2expr1000 <- models(data.snp[1:500+500],data.expr))
# system.time(snp2expr1500 <- models(data.snp[1001:length(data.snp)],data.expr))
# rlt.snp2expr <- c(snp2expr500,snp2expr1000,snp2expr1500)
# save(rlt.snp2expr,file='rlt.snp2expr.rda')
# system.time(snp2m <- models(data.snp,data.methy))
#save(snp2methy,file='rlt.snp2methy.rda')

###############
#
###############

setwd('/Users/wenrurumon/Documents/uthealth/path_via_gene')
load('genenet.rda'); load('pathnet.rda')
library(igraph)
library(data.table)
library(dplyr)
g2d <- read.csv('gene2disease.csv',stringsAsFactors=F)
g2d <- as.matrix(g2d)
genenet <- as.matrix(genenet)
genenet <- unique(rbind(genenet,g2d))
pathnet <- unique(rbind(pathnet,g2d))
pathnet <- gsub("ph':","ph'",pathnet)
pathnet <- gsub("ph'AGE","ph'Age",pathnet)
node.p <- names(geneincluster)
node.ph <- unique(grep("ph'",as.vector(pathnet),value=T)); node.ph <- substr(node.ph,4,nchar(node.ph))
node.d <- unique(grep("d'",as.vector(pathnet),value=T)); node.d <- substr(node.d,3,nchar(node.d))
pathnet <- gsub("g'","s:",pathnet)
pathnet <- gsub("m'","m:",pathnet)
pathnet <- gsub("d'|ph'|p'","",pathnet)
genenet <- gsub("m::","m:",genenet)
genenet <- gsub("s::","s:",genenet)
genenet <- rbind(genenet,c('s:CREBBP','CBL'))
gnet <- as.data.table(genenet); colnames(gnet) <- c('from','to')
pnet <- as.data.table(pathnet)
gg <- graph_from_data_frame(gnet)
pg <- graph_from_data_frame(pnet)

setwd('..')
setwd('output')
jr <- lapply(dir(),read.table,header=T,fill=T)
jr <- lapply(jr,function(x){
  filter(x,ANM<=0.01)
})
jrnames <- names(data.snp)[as.numeric(gsub('SNP_','',substr(dir(),1,regexpr('expression',dir())-2)))]
jrnames <- rep(jrnames,sapply(jr,nrow))
jrsel <- data.frame(snp=jrnames,do.call(rbind,jr))
jrsel <- filter(jrsel,ANM==0)

snp2 <- paste(jrsel$snp)
snp2 <- paste0('s',substr(snp2,regexpr(':',snp2)+1,nchar(snp2)))
exp2 <- paste(jrsel$exName)
gnet2 <- data.table(from=snp2,to=exp2)
gnet2 <- rbind(gnet,gnet2)
gg2 <- graph_from_data_frame(gnet2)
snp2 <- unique(snp2[snp2%in%cand]);length(snp2)

jrnames <- paste0('s',substr(jrnames,regexpr(':',jrnames)+1,nchar(jrnames)))

test <- sapply(snp2,function(x){
  c(AD=length(all_shortest_paths(gg2,x,'AD')[[1]]),DIA=length(all_shortest_paths(gg2,x,'Diabete')[[1]]))
})
table(AD=test[1,]>0,T2D=test[2,]>0)

####
cc2 <- unique(substr(paste(jrsel$snp),regexpr(':',paste(jrsel$snp))+2,nchar(paste(jrsel$snp))))
cc2 <- cbind(gene=gsub('s:','',cand),
             intest=gsub('s:','',cand)%in%cc2,
             inchenchao=gsub('s:','',cand)%in%cc)
write.csv(cc2,'cc2.csv')

test <- lapply(dir(),function(x){
  print(x)
  x <- read.table(x,header=T,fill=T)
  table(paste(x[,2]<0.05,x[,4]==0))
})
