Sys.setenv(RETICULATE_PYTHON = "/nas/longleaf/home/wancen/.conda/envs/scdali/bin/python")
library(reticulate)
library(SingleCellExperiment)
library(airpart)
library(pbapply)
library(tidyverse)
# use_python("C:/Users/wancenmu/anaconda3/envs/scdali/python.exe", required = TRUE)
source_python("/proj/milovelab/mu/SC-ASE/simulation/scdali/test.py")
# import numpy and specify no automatic Python to R conversion
np <- import("numpy", convert = FALSE)
nct <-6
n <- 40
set.seed(37)
p.vec1 <- rep(c(0.5, 0.5,0.5,0.5,0.6,0.7),each=2)
# p.vecc <-vctrs::vec_c(p.vec1, p.vec2)

# sig_true <- c(rep(1,100),rep(c(2,3),50))
sce <- makeSimulatedData(
  mu1 = 2,
  mu2 = 10,
  nct = nct,
  n = n,
  ngenecl = 400,
  theta = 20,
  ncl =2,
  p.vec = p.vec1)
sce <- preprocess(sce)
rowData(sce)$cluster <- seq_len(dim(sce)[1])
# sig_true <-(rowData(sce)[,1:nct]) %>% as.data.frame()


# scdali
ase <- assays(sce)[["a1"]]
cts <- assays(sce)[["counts"]]
x<-rep(seq_len(nct),each=n)
np.x<-np_array(x,dtype = "float64")
np.ase<-np_array(t(ase),dtype = "float64")
np.cts<-np_array(t(cts),dtype = "float64")
t_dali <- system.time(scdali <-testsig(np.ase,np.cts,np.x,1))[[3]]
scdali_res <- ifelse(scdali<0.1,TRUE,FALSE)
true <- rep(c(1,2),each=400)
table(scdali_res,true)/400


#fused lasso
t_fl_binomial <- system.time(
  fl_binomial <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      sce_sub <- fusedLasso(sce,model="binomial",genecluster = i,ncores=1,se.rule.mult = 0.5)
      index <- nlevels(metadata(sce_sub)$partition$part)>1
      return(index)
    }, error = function(e) {
      sce_sub <- fusedLasso(sce,model="gaussian",genecluster = i,ncores=1,se.rule.mult = 0.5)
      index <- nlevels(metadata(sce_sub)$partition$part)>1
      return(index)
    })
  },cl=15))[[3]]
table(fl_binomial,true)/400
#
# # fl_binomial <- do.call(rbind, fl)
# #fused lasso
t_fl_gaussian <- system.time(
  fl_gaussian <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      sce_sub <- fusedLasso(sce,model="gaussian",genecluster = i,ncores=1)
      index <- nlevels(metadata(sce_sub)$partition$part)>1
      return(index)
    }, error = function(e) {return(NA)})
  },cl=15))[[3]]
#
# # wilcoxon
t_wilcoxon <- system.time(
  wilcoxon <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      sce_sub <- wilcoxExt(sce,genecluster = i)
      index <- nlevels(metadata(sce_sub)$partition$part)>1
      return(index)
    }, error = function(e) {return(NA)})
  },cl=15))[[3]]

save(p.vec1,sce,fl_gaussian,t_fl_gaussian,fl_binomial,t_fl_binomial,
     scdali,scdali_res,t_dali,wilcoxon,t_wilcoxon,file = "/proj/milovelab/mu/SC-ASE/simulation/scdali/fp.rda")
