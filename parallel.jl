
require("/home/nlin/schizaphrenia/lrpca.screening/lrcca/new/source/ridge_cca.jl");
p_value = SharedArray(Float64,(18026,2));

@parallel for i = 1:18026
   println(i)
   A = readdlm(string(i,".txt"));
   B = readdlm("pheno.txt");
   p_chi,p_f = ridge_cca(A,B,0.1);
   p_value[i,1] = p_chi;
   p_value[i,2] = p_f;  
end
