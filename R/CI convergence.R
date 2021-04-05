n=40
ngenecl=5
# source("/proj/milovelab/mu/SC-ASE/simulation/fusedlasso.R")
source("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/scripts/fusedlasso.R")
library(smurf)
library(pbapply)
library(mclust)
library(emdbook)
library(fastmatch)
library(clue)
library(airpart)

ans <- pbsapply(1:2, function(i) {
  set.seed(i)
  p.vec <- c(0.95, 0.9, 0.85, 0.85, 0.7, 0.7, 0.65, 0.6, 0.5, 0.5)
  sce <- makeSimulatedData(mu1=2,mu2=mu2,nct=length(p.vec),n=n,ngenecl=ngenecl,theta=20,ncl=1,p.vec =p.vec)
  sce <- preprocess(sce)

  rowdata<-cbind(rowData(sce),cluster=rep(1,ngenecl))
  rowData(sce)<-rowdata
  f <- ratio ~ p(x, pen="gflasso")
  sce$part<-factor(rep(c(1,2,3,3,4,4,5,6,7,7),each=n))
  cl <- data.frame(part=c(1,2,3,3,4,4,5,6,7,7), x = levels(sce$x))
  metadata(sce)$partition <- cl
  # Normal approximation
  sce <- allelicRatio(sce)
  ci<-metadata(sce)$estimator[,4:5]
  index_bb<-ifelse(sapply(1:10, function(j)
    any(ci[j,1] <= p.vec[j] & ci[j,2] >= p.vec[j])),1, 0)
  # Bootstrap
  sce2 <- allelicRatio(sce,method = "bootstrap", R=200)
  ci<-metadata(sce2)$estimator[,4:5]
  index_boot<-ifelse(sapply(1:10, function(j)
    any(ci[j,1] <= p.vec[j] & ci[j,2] >= p.vec[j])),1, 0)

  out <- list(index_bb,index_boot)
  out
},cl=15)
# }

# save the results as a data.frame
normal<-do.call(rbind, ans[seq(1,length(ans), by = 2)])
boot<-do.call(rbind, ans[seq(2,length(ans), by = 2)])


colSums(normal)
# [1] 192 188 188 188 186 186 192 193 189 189
colSums(boot)
# save the results as a data.frame
dat <- data.frame(normal=colSums(normal),
                  boot=colSums(boot),
                  n=ngenecl,
                  cnt=mu2)

# write out as a table
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
