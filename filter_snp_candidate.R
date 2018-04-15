rm(list=ls())

load('cutdata.rda')
gene <- names(data.snp)
gene <- substr(gene,regexpr('::',gene)+2,nchar(gene))

files <- dir(pattern='processed')

rlt <- lapply(files,function(f){
print(f)
load(f)
x <- lapply(psnp[names(psnp)%in%gene],function(x){
	list(snp=x$snp,pos=x$pos)
})
names(x) <- names(psnp[names(psnp)%in%gene])
return(x)
})

names(rlt) <- files
save(rlt,file='rawdata4snpcandidate.rda')

#scp zhu@129.106.3.110:/home/zhu/rushdata/processed/snp/rawdata4snpcandidate.rda .
