
#Check id
rm(list=ls())
setwd('/home/zhu/rushdata/processed/methylation')
load(dir()[1])
setwd('/home/zhu/rushdata/processed/snp')
load(dir()[1])
setwd('/home/zhu/rushdata/expression_clustering')
load('expr_cluster_list.rda')
setwd('/home/zhu/rushdata/expression_clustering')
load('disease.rda')

mid <- rownames(pmethy[[1]][[1]])
sid <- gsub('ROS|MAP','',rownames(psnp[[1]][[1]]))
eid <- lapply(exprincluster,function(x){rownames(x[[1]])})
did <- rownames(disease)

pid <- unlist(list(mid,sid,eid,did))
pid <- names(which(table(pid)==5))
save(pid,file='pid.rda')

#methylation
rm(list=ls())
load('/home/zhu/rushdata/expression_clustering/pid.rda')
setwd('/home/zhu/rushdata/processed/methylation')
raw <- lapply(dir(),function(x){
	print(x)
	load(x)
	x <- gsub('[a-zA-Z|_|.]','',x)
	names(pmethy) <- paste0(x,'::',names(pmethy))
	pmethy
})
raw2 <- do.call(c,raw)
methylation_fpca <- lapply(raw2,function(x){
	x.pid <- rownames(x$methy)
	x <- x$q2fpca$score
	x[match(pid,x.pid),,drop=F]
})
methylation_site <- lapply(raw2,function(x){
	x.pid <- rownames(x$methy)
	x <- x$methy
	x[match(pid,x.pid),,drop=F]
})
setwd('/home/zhu/rushdata/tacc')
save(methylation_site,methylation_fpca,file='methylation4cut.rda')

#SNP
rm(list=ls())
load('/home/zhu/rushdata/expression_clustering/pid.rda')
setwd('/home/zhu/rushdata/processed/snp')
raw <- lapply(dir(),function(x){
	print(x)
	load(x)
	x <- gsub('[a-zA-Z|_|.]','',x)
	names(psnp) <- paste0(x,'::',names(psnp))
	psnp
})
raw2 <- do.call(c,raw)
raw2 <- raw2[sapply(raw2,length)==5]
snp_fpca <- lapply(raw2,function(x){
	x.pid <- rownames(x$snp)
	x <- x$q2fpca$score
	x[match(pid,x.pid),,drop=F]
})
snp_site <- lapply(raw2,function(x){
	x.pid <- rownames(x$snp)
	x <- x$snp
	x[match(pid,x.pid),,drop=F]
})
setwd('/home/zhu/rushdata/tacc')
save(snp_site,snp_fpca,file='snp4cut.rda')

#expression
rm(list=ls())
setwd('/home/zhu/rushdata/expression_clustering')
load('expr_cluster_list.rda')
load('/home/zhu/rushdata/expression_clustering/pid.rda')

raw2 <- do.call(c,exprincluster)
expr <- lapply(raw2,function(x){
	x.pid <- rownames(x)
	x[match(pid,x.pid),,drop=F]
})
expr <- do.call(cbind,expr)
expr <- t(unique(t(expr)))
save(expr,file='expr4cut.rda')
