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
nct <-8
n <- 60
set.seed(37)
ardiff<- c(0.05,0.1,0.2,0.3)
p.vec <- lapply(1:length(ardiff), function(i) {
  rep(c(0.3,0.3+ardiff[i],0.3+2*ardiff[i],0.3+ardiff[i]),each=2)
})
p.vec1 <- do.call(c,p.vec)
# p.vec2 <- rep(rep(c(0.5,0.7),each=8), 5)
# p.vecc <-vctrs::vec_c(p.vec1, p.vec2)

# sig_true <- c(rep(1,100),rep(c(2,3),50))
sce <- makeSimulatedData(
  mu1 = 2,
  mu2 = 10,
  nct = nct,
  n = n,
  ngenecl = 400,
  theta = 20,
  ncl = length(ardiff),
  p.vec = p.vec1)
sce <- preprocess(sce)
rowData(sce)$cluster <- seq_len(dim(sce)[1])
# sig_true <-(rowData(sce)[,1:nct]) %>% as.data.frame()


# scdali
ase <- assays(sce)[["a1"]]
cts <- assays(sce)[["counts"]]
x<-rep(seq_len(nct),each=n)
np.x<-np_array(x,dtype = "int")
np.ase<-np_array(t(ase),dtype = "float64")
np.cts<-np_array(t(cts),dtype = "float64")
loc <-testsig(np.ase,np.cts,np.x,1)
# x_start<-seq(1,240,30)
t_dali <- system.time(
      scdali <-practice(np.ase,np.cts,np.x,ncores=1))[[3]]
# sig<-obj[[1]]+1
# obj2<-obj[[2]][x_start,]
# diff_scdali <- sig_true-t(obj2)
# rmse_scdali<-sqrt((rowSums(diff_scdali)^2)/nct)

#fused lasso
t_fl_binomial <- system.time(
  fl_binomial <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      if(i<=400){
        sce_sub <- fusedLasso(sce,model="binomial",genecluster = i,ncores=1,se.rule.mult = 0.5)
      }else{
        if(i>400 & i<=800){
          sce_sub <- fusedLasso(sce,model="binomial",genecluster = i,ncores=1,se.rule.mult = 0.5)
        }else{
          sce_sub <- fusedLasso(sce,model="binomial",genecluster = i,ncores=1)
        }
      }
      # index <- nlevels(metadata(sce_sub)$partition$part)>1
      sce_sub <- allelicRatio(sce_sub)
      est <- extractResult(sce_sub)%>% as.matrix()
      up <- extractResult(sce_sub,estimates = "upper")%>% as.matrix()
      low <- extractResult(sce_sub,estimates = "lower")%>% as.matrix()
      return(c(est,up,low))
    }, error = function(e) {return(rep(NA,3*nct)) })
  },cl=15))[[3]]
#
t_fl <- system.time(
  fl <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      # sce_sub <- fusedLasso(sce,model="binomial",genecluster = i,ncores=1)
      sce_sub <- sce[i,]
      sce_sub$part <- factor(rep(seq_len(nct),each=n))
      # index <- nlevels(metadata(sce_sub)$partition$part)>1
      sce_sub <- allelicRatio(sce_sub)
      est <- extractResult(sce_sub)%>% as.matrix()
      up <- extractResult(sce_sub,estimates = "upper")%>% as.matrix()
      low <- extractResult(sce_sub,estimates = "lower")%>% as.matrix()
      return(c(est,up,low))
    }, error = function(e) {return(rep(NA,3*nct)) })
  },cl=15))[[3]]
# # fl_binomial <- do.call(rbind, fl)
# #fused lasso
t_fl_gaussian <- system.time(
  fl_gaussian <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      if(i<=1200){
          sce_sub <- fusedLasso(sce,model="gaussian",genecluster = i,ncores=1,se.rule.mult = 0.5)
        }else{
          sce_sub <- fusedLasso(sce,model="gaussian",genecluster = i,ncores=1)
        }
      # index <- nlevels(metadata(sce_sub)$partition$part)>1
      sce_sub <- allelicRatio(sce_sub)
      est <- extractResult(sce_sub)%>% as.matrix()
      up <- extractResult(sce_sub,estimates = "upper")%>% as.matrix()
      low <- extractResult(sce_sub,estimates = "lower")%>% as.matrix()
      return(c(est,up,low))
    }, error = function(e) {return(rep(NA,3*nct)) })
  },cl=15))[[3]]
#
# # wilcoxon
t_wilcoxon <- system.time(
  wilcoxon <- pbsapply(1:nrow(sce), function(i) {
    res <- tryCatch({
      sce_sub <- wilcoxExt(sce,genecluster = i)
      sce_sub <- allelicRatio(sce_sub)
      est <- extractResult(sce_sub)%>% as.matrix()
      up <- extractResult(sce_sub,estimates = "upper")%>% as.matrix()
      low <- extractResult(sce_sub,estimates = "lower")%>% as.matrix()
      return(c(est,up,low))
    }, error = function(e) {return(rep(NA,3*nct)) })
  },cl=15))[[3]]

# save(ardiff,loc,sce,fl_gaussian,t_fl_gaussian,fl_binomial,t_fl_binomial,
#      scdali,t_dali,wilcoxon,t_wilcoxon,file = "/proj/milovelab/mu/SC-ASE/simulation/scdali/n200_n40_rbf.rda")
# load("/proj/milovelab/mu/SC-ASE/simulation/scdali/scdali.rda")
save(ardiff,loc,sce,fl_gaussian,t_fl_gaussian,fl_binomial,t_fl_binomial,fl,t_fl,
     scdali,t_dali,wilcoxon,t_wilcoxon,file = "/proj/milovelab/mu/SC-ASE/simulation/scdali/default_n60_rbf.rda")
