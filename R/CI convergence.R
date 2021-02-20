cmd_args=commandArgs(TRUE)

n <- as.numeric(cmd_args[1]) # cells per cell type
cnt <- as.numeric(cmd_args[2]) # the mean "high" total count
out <- cmd_args[3]


k<-10
# source("/proj/milovelab/mu/SC-ASE/simulation/fusedlasso.R")
source("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/scripts/fusedlasso.R")
library(smurf)
library(pbapply)
library(mclust)
library(emdbook)
library(fastmatch)
library(clue)

ans <- pbsapply(1:200, function(i) {
  set.seed(i)
  k <- k # number of cell types
  low_count <- 2 # the mean "low" total count
  mean_total_count <- rep(rep(c(low_count, cnt),each=n/2), times=k) # total count
  size <- rpois(n * k, mean_total_count)
  # size[size == 0] <- 1
  p.vec <- c(0.95, 0.9, 0.85, 0.85, 0.7, 0.7, 0.65, 0.6, 0.5, 0.5)
  # p.vec <- (5 + rep(seq(from=3.5,to=4.5,length.out=k/2),each=2))/10
  p <- rep(p.vec, each=n) # true prob
  #y <- rbinom(k*n, prob=p, size=size) # obs counts
  y <- rbetabinom(k*n, prob=p, size=size, theta=20) # obs counts
  ratio <- y/size # ratio
  x <- factor(rep(1:k,each=n)) # cell type dummy
  f <- ratio ~ p(x, pen="gflasso", refcat="1") # formula
  dat=data.frame(x=x,ratio=ratio,cts=size)
  # t <- system.time(fit<-fusedlasso(formula=f,model="quasibinomial",data=dat,ncores=1))[[3]] # quasibinomial
  # t5<-system.time(fit2<-wilcox_adj(dat,nct=k,k=5,threshold =10^seq(from=-2,to=-0.4,by=0.1)))[[3]]
  
  # co <- coef(fit)
  # co <- co + c(0,rep(co[1],k-1))
  class<-fmatch(p.vec,unique(p.vec))
  label<-tibble(type=factor(seq_len(10)),par_wilcoxon=factor(class))
  dat2<-dat %>% left_join(label,by=c("x"="type"))
  
  # test<-tryCatch({
  #   vglm(cbind(ratio*cts, cts-ratio*cts) ~par_wilcoxon, betabinomial, data = dat2, trace = F)}, error=function(e) {
  # vglm(cbind(ratio*cts, cts-ratio*cts) ~ 1, betabinomial, data = dat2, trace = F)}) # beta binomial confidence interval
  #subset = cts>1
  estimator<-sapply (1:max(class), function(m){
  bb<-vglm(cbind(ratio*cts, cts-ratio*cts) ~1, betabinomial, data = dat2[which(dat2$par_wilcoxon==m),], trace = F)
  coef_bb<-Coef(bb)[-2] # betabinomial estimator
  rho<-Coef(bb)[2]
  confint_bb<-confintvglm(bb,matrix=T)[-2,]
  confint_wilcoxon<-1/(1+exp(-confint_bb))
  return(list(coef_bb,confint_wilcoxon,rho))
  })
  
  # allelic ratio estimator
  coef<-as.vector(do.call(rbind, estimator[seq(1,length(estimator), by = 3)]))
  est<-tibble(par_wilcoxon=factor(seq_len(length(coef))),estimator=coef)
  est2<-merge(label, est, by="par_wilcoxon", all = T) #combine with partation label
  print(knitr::kable(est2))
  
  # Normal approximation
  confint<-matrix(do.call(rbind, estimator[seq(2,length(estimator), by = 3)]),ncol=2) %>%as_tibble() %>% setNames(names(estimator[[2]]))
  normal<-confint[as.numeric(class),]
  print(knitr::kable(normal))
  
  # Bootstrap
  

  out <- list(index_vgam,index_mm,index_mle)
  out
},cl=15)
# }

# save the results as a data.frame
vgam<-do.call(rbind, ans[seq(1,length(ans), by = 3)])
mm<-do.call(rbind, ans[seq(2,length(ans), by = 3)])
mle<-do.call(rbind, ans[seq(3,length(ans), by = 3)])

colSums(vgam)
# [1] 192 188 188 188 186 186 192 193 189 189
colSums(mm)
# [1] 176 173 191 191 189 190 168 173 196 197
colSums(mle)
# [1] 173 172 190 192 190 191 170 167 191 195


# write out as a table
# write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
