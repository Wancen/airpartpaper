cmd_args=commandArgs(TRUE)

ngenecl <- as.numeric(cmd_args[1]) # cells per cell type
mu2 <- as.numeric(cmd_args[2]) # the mean "high" total count
out <- cmd_args[3]

n=40
library(pbapply)
library(airpart)
library(mclust)
library(tidyverse)
library(smurf)
library(SingleCellExperiment)
source("/proj/milovelab/mu/SC-ASE/simulation/fusedlasso.R")
ans <- pbsapply(1:200, function(i) {
  set.seed(i)
  p.vec <- c(0.95, 0.9, 0.85, 0.85, 0.7, 0.7, 0.65, 0.6, 0.5, 0.5)
  sce <- makeSimulatedData(mu1=2,mu2=mu2,nct=length(p.vec),n=n,ngenecl=ngenecl,theta=20,ncl=1,p.vec =p.vec)
  # hist(assays(se)[["counts"]],xlim = c(0,60),breaks = 50)
  # se <- makeSimulatedData2(mu1=1,mu2=mu2,nct=length(p.vec),n=n,ngenecl=ngenecl,theta=20,ncl=1,p.vec =p.vec)
  sce<-preprocess(sce)
  rowdata<-cbind(rowData(sce),cluster=rep(1,ngenecl))
  rowData(sce)<-rowdata
  f <- ratio ~ p(x, pen="gflasso")
  t1 <- system.time(sce_sub <- airpart::fusedLasso(sce,formula=f,model="binomial",genecluster = 1,ncores=1,se.rule.nct = 10))[[3]]
  a1 <- adjustedRandIndex(factor(p.vec), metadata(sce_sub)$partition$part)
  t2 <- system.time(sce_sub <- airpart::fusedLasso(sce,formula=f,model="gaussian",genecluster = 1,ncores=1,se.rule.nct = 10))[[3]]
  a2 <- adjustedRandIndex(factor(p.vec), metadata(sce_sub)$partition$part)
  t3 <- system.time(sce_sub <- airpart::wilcoxExt(sce, genecluster=1 ))[[3]]
  a3 <- adjustedRandIndex(factor(p.vec), metadata(sce_sub)$partition$part)
  # t4 <- system.time(se_sub <- wilcoxExt2(se, genecluster=1))[[3]]
  # a4 <- adjustedRandIndex(factor(p.vec), metadata(se_sub)$partition$part)
  out <- c(a1,a2,a3,t1,t2,t3)
  out
},cl=15)

# dealing with the tryCatch errors...
if (!is.matrix(ans)) {
  ans <- do.call(cbind, ans)
}
# save the results as a data.frame
dat <- data.frame(type=rep(c("fl_binomial","fl_gaussian","WilcoxExt"),each=ncol(ans)),
                  ARI=as.vector(t(ans[1:3,])),
                  time=as.vector(t(ans[4:6,])),
                  n=ngenecl,
                  cnt=mu2)

# write out as a table
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
