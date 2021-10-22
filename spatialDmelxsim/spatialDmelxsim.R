# The way we pre-filter the genes by QC and LRT and save the result in "./data/spatial.RData"

# library(spatialDmelxsim)
# sce <- spatialDmelxsim()
# # Find correlation
# library(SingleCellExperiment)
# library(airpart)
# library(scater)
# sce <- as(spatialDmelxsim, "SingleCellExperiment")
# assay(sce,"a1")[is.na(assay(sce,"a1"))] <- 0
# assay(sce,"a2")[is.na(assay(sce,"a2"))] <- 0
# sce <- preprocess(sce)
# ## QC on genes with 10 ASE counts in at least half of the samples
# keep_feature2<-nexprs(assay(sce,"counts"),byrow = TRUE, detection_limit = 20) >= ncol(sce)/2
# table(keep_feature2) ## Keep 5026 genes
#
# sce<-sce[keep_feature2,] ## keep 5026
#
# ## perform LRT
# library(lmtest)
# library(pbapply)
# library(tidyverse)
# selectgene <- pbsapply(seq_len(nrow(sce)),function(i){
#   data <- data.frame(ratio = assay(sce,"ratio")[i,] %>% unlist(),
#                      slice = sce$slice,
#                      strain = sce$strain)
#   nested <- glm(ratio~slice,data=data)
#   complex <- glm(ratio~slice + strain,data=data)
#   fit<-lrtest(nested, complex)
#   psig <- fit$`Pr(>Chisq)`[2] > 0.01
#   return(psig)
# })
# table(selectgene) # 299 gene
# selectgene[which(is.na(selectgene))]<-FALSE
# sce <- sce[selectgene,]

## Intergrating slices by correlation
library(plyr)
library(SingleCellExperiment)
library(airpart)
library(scater)
library(tidyverse)
## rename x levels
# sce <- sce[,-c(42,74,115)]
# sce$x <- factor(c(1,1,2,3,3,4,4,5,5,6,6,6,7,7,8,8,9,10,11,12,13,14,15,16,17,18,19,
#                   1,1,2,3,4,4,4,5,6,6,6,7,7,8,9,10,11,12,13,14,15,16,17,18,19,
#                   1,1,2,3,4,4,5,6,6,7,7,8,9,10,11,11,11,11,12,13,14,15,16,17,18,19,
#                   1,1,1,2,3,3,4,5,6,6,7,7,8,8,9:19,
#                   1,1,2,3,3,4,4,4,5,6,6,7,7,8,8,9:19))
# # se$x <- mapvalues(se$x, from = levels(se$x), to = seq_len(12))
# cell_meta_unique <- seq_len(19)
# order_of_development <- order(match(sce$x,cell_meta_unique))
# sce <- sce[,order_of_development]
#
# cellQCmetrics<-cellQC(sce,mad_detected=3)
# cellQCmetrics
# keep_cell <- (
#   cellQCmetrics$filter_sum  # sufficient features (genes)
# )
# table(keep_cell)
# ## keep 126 out of 129 cells
# sce<-sce[,keep_cell]
#
# ## QC on genes with 25% of slice expressed within each slice group
# featureQCmetric<-featureQC(sce,sd = 0.025,detection_limit = 20)
# keep_feature<-(featureQCmetric$filter_celltype &
#                  featureQCmetric$filter_sd)
# table(keep_feature) ## Keep 222 genes
# sce<-sce[keep_feature,] ## keep 222

## helper function
plotGene <- function(gene) {
  # x <- sce$normSlice
  y <- assay(sce, "ratio")[gene,]
  Q <- quantile(y, probs=c(.25, .75), na.rm = T)
  iqr <- IQR(y, na.rm = T)
  up <-  min(Q[2]+1.5*iqr,1) # Upper Range
  low<- max(0,Q[1]-1.5*iqr) # Lower Rangeÿ
  dat <- data.frame(
    ratio = as.vector(unlist(y)),
    x = rep(sce$x, each = nrow(y)) %>% as.numeric(),
    feat = factor(rep(gene, ncol(y)), levels = gene),
    rep = factor(rep(sce$rep, each = nrow(y)))
  )
  # dat <- data.frame(slice = x, rep = as.factor(sce$rep),ratio=y)
  ggplot(dat,aes(x=x,y=ratio)) +
    geom_point(aes(color=rep)) +
    geom_smooth(method = "loess", se = FALSE) +
    theme_minimal() +
    geom_hline(yintercept=0.5,linetype=2) +
    ylim(low,up) +
    scale_color_brewer(palette = "Paired")+
    facet_wrap(~feat) +
    theme(strip.text = element_text(size=20),
          legend.position = "none",panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = 'transparent'),
          axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
    labs(x = "grouped slices", y = "allelic ratio") +
    scale_x_continuous(breaks = seq_len(19))
}

## airpart analysis
load(file="./data/spatial.RData")
sce<-geneCluster(sce, G=60)
# sce_hclust<-geneCluster(sce, method = "hierarchical")
# ## Choose gene cluster 2,10,12,17
makeHeatmap(sce,genecluster = 10,show_row_names = T)

# clover5<-which(metadata(sce)$geneCluster>=5)
estDisp(sce,genecluster = 10)
library(smurf)
## construct graph \Gamma
f <- ratio ~ p(x, pen = "ggflasso")
nct <- nlevels(sce$x)
b <- matrix(0, nct, nct)
a <- diag(nct-1)
a[lower.tri(a,diag = T)]
b[lower.tri(b, diag = FALSE)] <- a[lower.tri(a,diag = T)]
b2 <- b + t(b)
colnames(b2) <- rownames(b2) <- levels(sce$x)

sce_sub <- fusedLasso(sce, formula=f,adj.matrix = list(x = b2),
                        model = "binomial", niter = 1,
                        genecluster = 10, ncores = 4,lambda.max=100
  )
metadata(sce_sub)
sce_sub <- allelicRatio(sce_sub)
## for visualization, remove gene "ICA69"
p <- makeStep(sce_sub,xlab = "grouped slices")
gene <- rownames(sce_sub)
q <- plotGene(gene[-6])

jpeg(file="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/plot/60_cl10airpart.jpeg",width = 8, height = 6,units = "in",res=360)
p
dev.off()
metadata(sce_sub)
sce_sub <- allelicRatio(sce_sub)
s<-extractResult(sce_sub,"svalue")
lower<-extractResult(sce_sub,"lower")
upper<-extractResult(sce_sub,"upper")
makeHeatmap(sce_sub, order_by_group = F,show_row_names = T)

## scdali
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

ase <- assays(sce_sub)[["a1"]]
cts <- assays(sce_sub)[["counts"]]
x<-sce_sub$x %>% as.numeric()
np.x<-np_array(x,dtype = "float64")
np.ase<-np_array(t(ase),dtype = "float64")
np.cts<-np_array(t(cts),dtype = "float64")
t_dali <- system.time(
  scdali <- pbsapply(1:nrow(sce_sub), function(i) {
    np.ase<-np_array(t(ase[i,]),dtype = "float64")
    np.ase <- array_reshape(np.ase, c(ncol(sce_sub),1))
    np.cts<-np_array(t(cts[i,]),dtype = "float64")
    np.cts <- array_reshape(np.cts, c(ncol(sce_sub),1))
    res <- tryCatch({
      obj <-practice(np.ase,np.cts,np.x,ncores=1)
      obj <- do.call(c,obj)
      return(obj)
    }, error = function(e) {return(rep(NA,2*ncol(sce_sub))) })
  }))[[3]]

values<-rle(x)$lengths
x_start<-cumsum(c(1,values))[1:nct]
ar_scdali<-scdali[1:ncol(sce_sub),]
sd_scdali<-scdali[(ncol(sce_sub)+1):(2*ncol(sce_sub)),]
mu<-ar_scdali[x_start,] %>% t() %>% `rownames<-`(rownames(sce_sub))
library(reshape2)
ar_long <- melt(mu, id = "slice", value.name = "ar") %>% `colnames<-`(c("gene","slice","ar"))
ar_long$slice <- factor(ar_long$slice)
raw <- assay(sce_sub, "ratio")
Q <- quantile(raw, probs=c(.25, .75), na.rm = T)
iqr <- IQR(raw, na.rm = T)
up <-  Q[2]+1.5*iqr # Upper Range
low<- Q[1]-1.5*iqr # Lower RangeÃ¿
scdali<-ggplot(ar_long %>% filter(gene!="ICA69"), aes(x = slice,y = ar)) +
  geom_point() +
  # geom_step(group = 1) +
  theme_minimal()+
  ylim(low,up) +
  ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
  theme(strip.text = element_text(size=20),
        legend.position="none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  labs(x = "grouped slices", y = "allelic ratio") +
  facet_wrap(~gene)


## Overall trend
sce <- sce_sub
part <- metadata(sce)$part
part$coef <- 1/(1+exp(-part$coef))
part$x <- factor(part$x,levels = seq_len(19))
raw <- assay(sce, "ratio")
Q <- quantile(raw, probs=c(.25, .75), na.rm = T)
iqr <- IQR(raw, na.rm = T)
up <-  min(Q[2]+1.5*iqr,1) # Upper Range
low<- max(0,Q[1]-1.5*iqr) # Lower Rangeÿ
dat <- data.frame(
  ratio = as.vector(unlist(ar)),
  x = factor(rep(unique(sce$x), each = length(sce))),
  part = factor(rep(metadata(sce)$partition$part, each = length(sce))),
  feat = factor(rep(row.names(sce), nlevels(sce$x)), levels = row.names(sce))
)
p <- ggplot2::ggplot(part, aes(x = .data$x, y = .data$coef)) +
  ggplot2::geom_point(aes(color = .data$part)) +
  ggplot2::geom_step(group = 1) +
  ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
  ylim(low,up) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = "transparent"),
    axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)
  ) +
  labs(x = "grouped slices", y = "allelic ratio")
p

## no group step, smoothing splines ##############################
sce_sub <- sce[rowData(sce)$cluster==10]
f <- ratio ~ ns(x, df = 5)
sce_sub <- allelicRatio(sce_sub, f)
ar <-extractResult(sce_sub) %>% as.matrix()

library(reshape2)
ar_long <- melt(ar, id = "slice", value.name = "ar") %>% `colnames<-`(c("gene","slice","ar"))
ar_long$slice <- factor(ar_long$slice)
nogroup<-ggplot(ar_long %>% filter(gene!="ICA69"), aes(x = slice,y = ar)) +
  geom_point() +
  # geom_step(group = 1) +
  theme_minimal()+
  ylim(low,up) +
  ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
  theme(strip.text = element_text(size=13),
        legend.position="none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  labs(x = "grouped slices", y = "allelic ratio") +
  facet_wrap(~gene)
jpeg(file="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/plot/suppspatial_nogroup.jpg",width = 8, height = 5.2,units = "in",res=450)
nogroup
dev.off()
# save(sce,file="C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/Github/airpartpaper/data/spatial.RData")


