n <- 50 # cells per cell type
mu2 <- 10 # the mean "high" total count
library(pbapply)
library(airpart)
library(mclust)
library(SummarizedExperiment)
ans <- pbsapply(1:100, function(i) {
  set.seed(i)
  p.vec <- c(0.95, 0.9, 0.85, 0.85, 0.7, 0.7, 0.65, 0.6, 0.5, 0.5)
  ngenecl=5
  se <- makeSimulatedData(mu1=1,mu2=mu2,nct=length(p.vec),n=n,ngenecl=ngenecl,theta=20,ncl=1,p.vec =p.vec)
  # se <- makeSimulatedData(mu1=1,mu2=10,nct=4,n=10,ngenecl=25,theta=20,ncl=1,p.vec =p.vec)
  makeRatioHeatmap(se)
  se<-preprocess(se)
  rowdata<-cbind(rowData(se),cluster=rep(1,ngenecl))
  rowData(se)<-rowdata
  f <- ratio ~ p(x, pen="gflasso") # formula
  t1 <- system.time(se_sub <- fusedLasso(formula=f,model="binomial",se,genecluster = 1))[[3]]
  a1 <- adjustedRandIndex(factor(p.vec), metadata(se_sub)$partition$part)
  t2 <- system.time(se_sub <- fusedLasso(formula=f,model="gaussian",se,genecluster = 1))[[3]]
  a2 <- adjustedRandIndex(factor(p.vec), metadata(se_sub)$partition$part)
  t3 <- system.time(se_sub <- wilcoxExt(se, genecluster=1, threshold =10^seq(from=-2,to=-0.4,by=0.1)))[[3]]
  a3 <- adjustedRandIndex(factor(p.vec), metadata(se_sub)$partition$part)
  t4 <- system.time(se_sub <- wilcoxExt2(se, genecluster=1, threshold =10^seq(from=-2,to=-0.4,by=0.1)))[[3]]
  a4 <- adjustedRandIndex(factor(p.vec), metadata(se_sub)$partition$part)
  out <- c(a1,a2,a3,a4,t1,t2,t3,t4)
  out
},cl=10)
# }

# dealing with the tryCatch errors...
if (!is.matrix(ans)) {
  ans <- do.call(cbind, ans)
}
# save the results as a data.frame
dat <- data.frame(type=rep(c("fl_binomial","fl_gaussian","WilcoxExt_mean","WilcoxExt_wmean"),each=ncol(ans)),
                  ARI=as.vector(t(ans[1:4,])),
                  n=n,
                  cnt=mu2)
