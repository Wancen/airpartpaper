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
true <- rep(c(FALSE,TRUE),each=400)
 table(scdali_res,true)/400


#fused lasso
t_fl_binomial <- system.time(
  fl_binomial <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      sce_sub <- fusedLasso(sce,model="binomial",genecluster = i,ncores=1,se.rule.mult = 0.5)
      sce_sub <- allelicRatio(sce_sub, DAItest = TRUE)
      p <- rowData(sce_sub)$p.value
      # index <- nlevels(metadata(sce_sub)$partition$part)>1
      return(p)
    }, error = function(e) {
      sce_sub <- fusedLasso(sce,model="gaussian",genecluster = i,ncores=1,se.rule.mult = 0.5)
      sce_sub <- allelicRatio(sce_sub, DAItest = TRUE)
      p <- rowData(sce_sub)$p.value 
      # index <- nlevels(metadata(sce_sub)$partition$part)>1
      return(p)
    })
  },cl=20))[[3]]
fl_binomial.adj <- p.adjust(fl_binomial,"fdr") <0.05
table(true,fl_binomial.adj)/400

# no group
t_fl <- system.time(
  fl <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      sce_sub <- sce[i,]
      sce_sub <- allelicRatio(sce_sub, nogroup = T, DAItest = T)
      p <- rowData(sce_sub)$p.value
      return(p)
    }, error = function(e) {return(NA)})
  },cl=20))[[3]]
fl.adj <- p.adjust(fl,"fdr") <0.05
table(true,fl.adj)/400


# #fused lasso
t_fl_gaussian <- system.time(
  fl_gaussian <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      sce_sub <- fusedLasso(sce,model="gaussian",genecluster = i,ncores=1)
      sce_sub <- allelicRatio(sce_sub, DAItest = TRUE)
      p <- rowData(sce_sub)$p.value 
      return(p)
    }, error = function(e) {return(NA)})
  },cl=20))[[3]]
fl_gaussian.adj <- p.adjust(fl_gaussian,"fdr") <0.05
table(true,fl_gaussian.adj)/400
#
# # wilcoxon
t_wilcoxon <- system.time(
  wilcoxon <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      sce_sub <- wilcoxExt(sce,genecluster = i)
      sce_sub <- allelicRatio(sce_sub, DAItest = TRUE)
      p <- rowData(sce_sub)$p.value 
      return(p)
    }, error = function(e) {return(NA)})
  },cl=20))[[3]]
wilcoxon.adj <- p.adjust(wilcoxon,"fdr") <0.05
table(true,wilcoxon.adj)/400
# save(p.vec1,sce,fl_gaussian,t_fl_gaussian,fl_binomial,t_fl_binomial,
#      scdali,scdali_res,t_dali,wilcoxon.adj,t_wilcoxon,file = "/proj/milovelab/mu/SC-ASE/simulation/scdali/fp100.rda")
save(p.vec1,sce,fl_gaussian.adj,t_fl_gaussian,fl_binomial.adj,t_fl_binomial,
     wilcoxon.adj,t_wilcoxon,file = "/proj/milovelab/mu/SC-ASE/simulation/scdali/fp40_DAI.rda")
