# Find correlation
sce <- as(se, "SingleCellExperiment")
assay(sce,"a1")[is.na(assay(sce,"a1"))] <- 0
assay(sce,"a2")[is.na(assay(sce,"a2"))] <- 0
sce <- preprocess(sce)
## QC on genes with 10 ASE counts in at least half of the samples
keep_feature2<-nexprs(assay(sce,"counts"),byrow = TRUE, detection_limit = 20) >= ncol(sce)/2
table(keep_feature2) ## Keep 5026 genes

sce<-sce[keep_feature2,] ## keep 5026

## perform LRT
library(lmtest)
library(pbapply)
selectgene <- pbsapply(seq_len(nrow(sce)),function(i){
  data <- data.frame(ratio = assay(sce,"ratio")[i,],
                     slice = sce$slice,
                     strain = sce$strain)
  nested <- glm(ratio~slice,data=data)
  complex <- glm(ratio~slice + strain,data=data)
  fit<-lrtest(nested, complex)
  psig <- fit$`Pr(>Chisq)`[2] > 0.01
  return(psig)
})
table(selectgene) # 299 gene
selectgene[which(is.na(selectgene))]<-FALSE
sce <- sce[selectgene,]

ratio_pseudo <- assay(sce2,"ratio_pseudo")
index <- order(rowVars(ratio_pseudo), decreasing=TRUE)
sce_select <- sce2[,]
sce_1 <- assay(sce_select[,which(sce_select$rep==1)],"ratio_pseudo")
sce_2 <- assay(sce_select[,which(sce_select$rep==2)],"ratio_pseudo")
sce_3 <- assay(sce_select[,which(sce_select$rep==3)],"ratio")
sce_4 <- assay(sce_select[,which(sce_select$rep==4)],"ratio_pseudo")
sce_5 <- assay(sce_select[,which(sce_select$rep==5)],"ratio_pseudo")
cor45<-cor(sce_4,sce_5)
cor34<-cor(sce_3,sce_4)
cor24 <- cor(sce_2,sce_4)
cor23<-cor(sce_2,sce_3)
cor13<-cor(sce_1,sce_3)
cor25<-cor(sce_2,sce_5)
cor15<-cor(sce_1,sce_5)



cor15<-cor(sce_1,sce_5)
diag(cor12)

sf6 <- c("CG43110", "veil", "hb", "CG4500", "CG4594", "Cyp4p2",
         "dream", "CG13384", "CG34266", "l(1)sc", "Adgf-A", "Cht3", "scw",
         "CG3502", "pxb", "CG15628", "CG8960", "CG43085", "path", "CG8147",
         "lea", "Bsg25D", "CG14915", "CG10035", "bmm", "Elba2", "CG14767",
         "prd", "CG15480", "bnb", "comm", "Pino", "fz2", "NimC4", "CG15479",
         "slp1", "mira", "CG14937", "CG14317", "unc-5", "ps", "CG12581", "rst",
         "Surf1", "mas", "Esp", "cnc", "uif", "CG12730", "Ance", "Bsg25A",
         "pcs", "CG13454", "CG43394", "CG9775", "sro", "CG2930",
         "l(2)08717", "scb", "srw", "CG30015", "sprt", "edl", "Mipp1",
         "CG9863", "Pepck")
plotGene(sf6[21])

# miss <- sf6[!sf6 %in% rownames(sce)]
# se$x <- cut(se$normSlice,12)
library(plyr)
library(SingleCellExperiment)
library(airpart)
library(scater)
library(tidyverse)
## rename x levels
sce <- sce[,-c(42,74,115)]
sce$x <- factor(c(1,1,2,3,3,4,4,5,5,6,6,6,7,7,8,8,9,10,11,12,13,14,15,16,17,18,19,
                  1,1,2,3,4,4,4,5,6,6,6,7,7,8,9,10,11,12,13,14,15,16,17,18,19,
                  1,1,2,3,4,4,5,6,6,7,7,8,9,10,11,11,11,11,12,13,14,15,16,17,18,19,
                  1,1,1,2,3,3,4,5,6,6,7,7,8,8,9:19,
                  1,1,2,3,3,4,4,4,5,6,6,7,7,8,8,9:19))
# se$x <- mapvalues(se$x, from = levels(se$x), to = seq_len(12))
cell_meta_unique <- seq_len(19)
order_of_development <- order(match(sce$x,cell_meta_unique))
sce <- sce[,order_of_development]

cellQCmetrics<-cellQC(sce,mad_detected=3)
cellQCmetrics
keep_cell <- (
  cellQCmetrics$filter_sum  # sufficient features (genes)
)
table(keep_cell)
## keep 126 out of 129 cells
sce<-sce[,keep_cell]

## QC on genes with 25% of slice expressed within each slice group
featureQCmetric<-featureQC(sce,sd = 0.025,detection_limit = 20)
keep_feature<-(featureQCmetric$filter_celltype &
                 featureQCmetric$filter_sd)
table(keep_feature) ## Keep 222 genes
sce<-sce[keep_feature,] ## keep 222

load(file="C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/Github/airpartpaper/data/spatial.RData")
sce<-geneCluster(sce, G=60)
sce_hclust<-geneCluster(sce, method = "hierarchical")
metadata(sce)$geneCluster # 3 gene cluster
metadata(sce[index,])$geneCluster # 3 gene cluster
table(sf6 %in% rownames(sce))
index <- which(rownames(sce) %in% sf6)
# rowData(sce[index,])
# rownames(sce)[index]
# rowData(sce[index,])$cluster[9]
table(rowData(sce[index,])$cluster)
# ## Choose gene cluster 10 and 36
makeHeatmap(sce,genecluster = 10,show_row_names = T)
# summary<-summaryAllelicRatio(sce)
# sapply(1:length(summary), function(i) {
#   inst <- summary[[i]]
#   inst_order <- inst[order(inst$weighted.mean), ]
#   max(diff(inst_order$weighted.mean)) > 0.05
# })

clover5<-which(metadata(sce)$geneCluster>=5)
estDisp(sce,genecluster = 27)
library(smurf)
f <- ratio ~ p(x, pen = "ggflasso")
nct <- nlevels(sce$x)
b <- matrix(0, nct, nct)
a <- diag(nct-1)
a[lower.tri(a,diag = T)]
b[lower.tri(b, diag = FALSE)] <- a[lower.tri(a,diag = T)]
b2 <- b + t(b)
colnames(b2) <- rownames(b2) <- levels(sce$x)
q<-list()
p<-list()
for (j in 1:length(clover5)) {
  sce_sub <- fusedLasso(sce, formula=f,adj.matrix = list(x = b2),
                        model = "binomial", niter = 1,
                        genecluster = 10, ncores = 4,lambda.max=100
  )
  metadata(sce_sub)
  sce_sub <- allelicRatio(sce_sub)
  p[[j]] <- makeStep(sce_sub[-6,],xlab = "grouped slices")
  gene <- rownames(sce_sub[-6,])
  q[[j]] <- plotGene(gene)
}

jpeg(file="C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/Github/airpartpaper/plot/60_cl10orig.jpeg",width = 8, height = 7,units = "in",res=360)
orig
dev.off()
metadata(sce_sub)
sce_sub <- allelicRatio(sce_sub)
s<-extractResult(sce_sub,"svalue")
lower<-extractResult(sce_sub,"lower")
upper<-extractResult(sce_sub,"upper")
makeHeatmap(sce_sub, order_by_group = F,show_row_names = T)
makeStep <- function(sce, xlab = "cell type") {
  ar <- rowData(sce)[, c(grep("ar", colnames(rowData(sce)), value = TRUE))] %>%
    `colnames<-`(levels(sce$x))
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
  p <- ggplot2::ggplot(dat, aes(x = .data$x, y = .data$ratio)) +
    ggplot2::geom_point(aes(color = .data$part)) +
    ggplot2::geom_step(group = 1) +
    ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
    ggplot2::facet_wrap(~feat) +
    ylim(low,up) +
    theme_minimal() +
    theme(
      strip.text = element_text(size=20),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = "transparent")
    ) +
    labs(x = xlab, y = "allelic ratio")
  p
}

makeStep(sce_sub,xlab = "grouped slices")

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
          panel.border = element_rect(fill = 'transparent')) +
    labs(x = "grouped slices", y = "allelic ratio") +
    scale_x_continuous(breaks = seq_len(19))
}
plotGene("CG4500")
plotGene("DOR")
plotGene("uif")
plotGene("bmm")
plotGene("Ndfip")
plotGene("CG13868")
plotGene("lea")
plotGene("ICA69")


gene <- rownames(sce_sub)
plotGene(gene)


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
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = "transparent")
  ) +
  labs(x = "grouped slices", y = "allelic ratio")
p

save(sce,file="C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/Github/airpartpaper/data/spatial.RData")
load(file="C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/Github/airpartpaper/data/spatial.RData")
