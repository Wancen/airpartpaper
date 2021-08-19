library(reticulate)
library(airpart)
library(SummarizedExperiment)
library(pbapply)
library(tidyverse)
use_python("C:/Users/wancenmu/anaconda3/envs/scdali/python.exe", required = TRUE)
source_python('C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/test.py')
# import numpy and specify no automatic Python to R conversion
np <- import("numpy", convert = FALSE)
nct <-9
set.seed(1)
ardiff<- seq(0.05,0.3,0.05)
p.vec <- lapply(1:length(ardiff), function(i) {
  rep(seq(from = 0.2, by=ardiff[i],length.out = 3),each=2)
})
p.vec1 <- do.call(c,p.vec)
# p.vec2 <- rep(rep(c(0.5,0.7),each=8), 5)
# p.vecc <-vctrs::vec_c(p.vec1, p.vec2)

sce <- makeSimulatedData(
  mu1 = 2,
  mu2 = 10,
  nct = nct,
  n = 30,
  ngenecl = 200,
  theta = 20,
  ncl = length(ardiff),
  p.vec = p.vec1)
sce <- preprocess(sce)
rowData(sce)$cluster <- seq_len(dim(sce)[1])
# sig_true <-(rowData(sce)[,1:nct]) %>% as.data.frame()


# scdali
ase <- assays(sce)[["a1"]]
cts <- assays(sce)[["counts"]]
x<-rep(seq_len(nct),each=30)
x_start<-seq(1,240,30)
np.ase<-np_array(t(ase),dtype = "float64")
np.cts<-np_array(t(cts),dtype = "float64")
np.x<-np_array(x,dtype = "float64")
t_dali <- system.time(obj <-practice(np.ase,np.cts,np.x))[[3]]
# sig<-obj[[1]]+1
obj2<-obj[[1]][x_start,]
diff_scdali <- sig_true-t(obj2)
rmse_scdali<-sqrt((rowSums(diff_scdali)^2)/nct)

#fused lasso
t_fl <- system.time(
fl <- pbsapply(1:nrow(sce), function(i) {
  res <- tryCatch({
    sce_sub <- fusedLasso(sce,model="binomial",genecluster = i,ncores=1)
  # index <- nlevels(metadata(sce_sub)$partition$part)>1
    sce_sub <- allelicRatio(sce_sub)
    est <- metadata(sce_sub)$estimator
    ar<- est$estimate
    se <- est$std.error
    return(c(ar,se))
  }, error = function(e) {return(rep(NA,2*nct)) })
},cl=8))[[3]]
ar <-fl[1:nct,]
se <-fl[(nct+1):(2*nct),]
diff_fl <- sig_true-t(ar)
rmse_fl<-sqrt((rowSums(diff_fl,na.rm = T)^2)/nct)

dat <- data.frame(rmse_fl=rmse_fl,
                     rmse_scdali=rmse_scdali,
                     type=factor(c(rep(round(ardiff,2),each=100))))
dat2<-dat %>% gather(key=method,value=rmse,rmse_fl:rmse_scdali)
#
library(ggplot2)
ggplot(dat2,mapping = aes(x=type,y=rmse,fill=method)) +
  geom_boxplot()

save(sce,fl,t_fl, t_dali,obj,file = "./data/n100_2group.rda")
load(file = "C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/scdali/n200_2group.rda")
