
using DataFrames
using LowRankModels
using Distributions
using MultivariateStats
using GLM

genes = readdlm("genes.txt")

function cca(i,j)
  gene1 = genes[i] #Input data
  gene2 = genes[j]
  methy = readdlm("methy_$gene1.txt")
  exp = readdlm("exp_$gene2.txt")
  if size(methy)[2] > size(exp)[1]
    chipvalue = NA
    fpvalue = NA
  else
    rho = correlations(MultivariateStats.fit(CCA,methy',exp')) #Test
    n = size(methy)[1]
    p = size(methy)[2]
    q = size(exp)[2]
    m = n-(3+p+q)/2
    s = sqrt(((p*q)^2-4)/(p^2+q^2-5))
    df1 = p*q
    df2 = m*s -p*q/2+1
    if df2 > 0
      lambda = prod(1-rho.^2)
      ftest = ((1-lambda^(1/s))*df2)/(lambda^(1/s)*df1)
      fpvalue = ccdf(FDist(df1,df2),ftest)
    else
      f_pvalue = 2
    end
    chitest = -1*n*sum(log(1-rho[1:min(p,q)].^2))
    chipvalue = ccdf(Chisq(p*q),chitest)
  end
  return chipvalue,fpvalue
end

function lrcca(i,j,thes=0.5)
  gene1 = genes[i] #Input data
  gene2 = genes[j]
  X = readdlm("methy_$gene1.txt")
  Y = readdlm("exp_$gene2.txt")
   xsvd = svd(X);
   ysvd = svd(Y);
   n = size(X)[1];
   p = size(X)[2];
   q = size(Y)[2];
  ###estimate the value for r and k for x and y seperately, the default value is that the sigular values that can explain more than 80% of variations##
   kx = sum(cumsum(xsvd[2])/sum(xsvd[2]) .< thes)+1;
   rx = xsvd[2][kx];
  ###fit qpca for x and y seperately
  if p  ==1
     methy = X;
  else
   xglrm = qpca(X,p,scale=rx);
   xgx,xgy,ch= LowRankModels.fit!(xglrm);
   methy = xgx' * xgy;
  end
  exp = Y;  
  if size(methy)[2] > size(exp)[1]
    chipvalue = NA
    fpvalue = NA
  else
    rho = correlations(MultivariateStats.fit(CCA,methy',exp')) #Test
    n = size(methy)[1]
    p = size(methy)[2]
    q = size(exp)[2]
    m = n-(3+p+q)/2
    s = sqrt(((p*q)^2-4)/(p^2+q^2-5))
    df1 = p*q
    df2 = m*s -p*q/2+1
    if df2 > 0
      lambda = prod(1-rho.^2)
      ftest = ((1-lambda^(1/s))*df2)/(lambda^(1/s)*df1)
      fpvalue = ccdf(FDist(df1,df2),ftest)
    else
      f_pvalue = 2
    end
    chitest = -1*n*sum(log(1-rho[1:min(p,q)].^2))
    chipvalue = ccdf(Chisq(p*q),chitest)
  end
  return chipvalue,fpvalue
end


