load(file="./data/dynamicASE.rda")

## Construct sce and QC
remove <- which(names(index)=="Intergenic")
a1_new<-a1_new[order_of_development,]
a2_new<-a2_new[order_of_development,]
assay.list<-list(a1=t(a1_new)[-remove,],a2=t(a2_new)[-remove,])
coldata<-data.frame(x=factor(id$V2[order_of_development],levels=cell_meta_unique),
                    id = factor(id$V1[order_of_development]))

rowdata<-merge(data.frame(V2=names(index),snp=snpname),geneGuess[which(!duplicated(geneGuess[,3])),2:3] %>% as.data.frame(),
               ) %>% `colnames<-`(c("geneName","snp","geneID"))

sce<-SingleCellExperiment(assays=assay.list, colData=coldata, rowData= DataFrame(rowdata[-remove,]))
sce <- preprocess(sce)

cellQCmetrics<-cellQC(sce,mad_detected=3,mad_sum = 4,threshold = 1,spike = "pseudo")
keep_cell <- (
  cellQCmetrics$filter_sum & # sufficient features (genes)
    cellQCmetrics$filter_detected & cellQCmetrics$filter_spike)
table(keep_cell)

sce<-sce[,keep_cell]

featureQCmetric<-featureQC(sce,sd = 0.02,detection_limit = 5,threshold = 0.25,
                           spike = "pseudo")
keep_feature<-(featureQCmetric$filter_celltype &
                 featureQCmetric$filter_sd & featureQCmetric$filter_spike)
table(keep_feature)
## selecting 1932 genes
sce<-sce[keep_feature,]

## selecting genes with high variance of weighted mean
sce_select <- sce
ratio_pseudo <- assay(sce_select,"ratio_pseudo")
ratio <- assay(sce_select,"ratio")
### Remove batch effect
y <- sce_select$id
ratio_pseudo_select <- limma::removeBatchEffect(ratio_pseudo, batch=y)
ratio_select <- limma::removeBatchEffect(ratio, batch=y)

assay(sce_select,"ratio_pseudo") <- ratio_pseudo_select
assay(sce_select,"ratio") <- ratio_select

mcols(sce_select)$cluster <- seq_len(nrow(sce_select))
sce_select<-summaryAllelicRatio(sce_select)
summary_select2 <- do.call(cbind,sce_select)
weighted_mean <- summary_select2[,seq(2,ncol(summary_select2),4)] %>% as.matrix()
varindex <- order(colVars(weighted_mean), decreasing=TRUE)
# index <- head(order(colVars(weighted_mean), decreasing=TRUE),50)
top50_name <- row.names(sce)[varindex]
top50<- weighted_mean[,varindex] %>% `colnames<-`(top50_name)

### See dynamic genes in paper
interest_gene <- c("HLA-DQB","DDX11","CASP8","PITRM1","F11R","UBASH3A","IL10","GNLY","CD3G")
for (gene in interest_gene) {
  print(which(str_detect(top50_name,gene)))
}

sce <- sce[varindex[1:50],]
mcols(sce)$cluster <- seq_len(nrow(sce))

### Modeling each gene with airpart
nct <- nlevels(sce$x)
b <- matrix(0, nct, nct)
a <- diag(nct-1)
a[lower.tri(a,diag = T)]
b[lower.tri(b, diag = FALSE)] <- a[lower.tri(a,diag = T)]
b2 <- b + t(b)
colnames(b2) <- rownames(b2) <- levels(sce$x)
f <- ratio ~ p(x, pen = "ggflasso") + id
library(smurf)
res <- pbsapply(seq_len(nrow(sce)),function(g){
  sce_sub <- fusedLasso(sce,formula = f,lambda = "is.bic",
                            model = "binomial",niter = 1,
                            genecluster = g, ncores = 1,
                            adj.matrix = list(x = b2),
                            lambda.max = 8000)
  if(nlevels(metadata(sce_sub)$partition$part)==1){
    ar <- rep(NA,nct)
  } else{
    sce_sub <- allelicRatio(sce_sub, formula = f)
    ar <- extractResult(sce_sub)
    ar %>% as.matrix()
  }
},cl=8)
colnames(res) <-rownames(sce)

save(res,file="./data/dynamicASE_top50.rda")
load("./data/dynamicASE_top50.rda")

## Cluster genes with same transcriptional cis regulatory programs: Late-Spike, Constant-Low, and Fluctuating
res <- res[ , apply(res, 2, function(x) !any(is.na(x)))]
geneVar_index <- order(colVars(res), decreasing=TRUE)
geneVar_name <- colnames(res)[geneVar_index]
for (gene in interest_gene) {
  print(which(str_detect(geneVar_name,gene)))
}
sce_cluster <- sce[geneVar_name[1:ncol(res)],]

#
library(tidyverse)
res2 <- res[,geneVar_name] %>% t()
# res2 <- res2 - rowMeans(res2)
#
# makeHeatmap(sce_cluster,genecluster = 1,show_row_names = T)
# makeHeatmap(sce_cluster,genecluster = 2,show_row_names = T)
# makeHeatmap(sce_cluster,genecluster = 3,show_row_names = T)
# makeHeatmap(sce_cluster,genecluster = 4,show_row_names = T)
## Combine all plot together
library("reshape2")
data <- list()
# test<-kmeans(res2,4)
sce_cluster <- sce_cluster[-17,]
mcols(sce_cluster)$cluster <- c(1,2,4,1,2,1,2,4,2,1,1,3,2,1,2,4,2,4,1,1,2,2,1,2,1,4,2,1,2,2,1,2,4,1,1,1,2,2,1,2,1,1,2)
# my_clusters <- c(3,5,1,1,3,1,1,4,1,3,1,3,5,5,5,4,1,6,5,5,3,1,3,3,3,1,4,3,5,6,4,3,4,2,1,6,1,1,1,6,2,2,2,5,3,5,2,6,2,5,1,6,2,2,1,4,3,4,1,3)
sce_cluster <- sce_cluster[c(str_which(rownames(sce_cluster),"F11R"),
                             str_which(rownames(sce_cluster),"HLA-DQB1"),
                             str_which(rownames(sce_cluster),"GNLY"),
                             str_which(rownames(sce_cluster),"HMGN4"),
                             str_which(rownames(sce_cluster),"DDX11"),
                             str_which(rownames(sce_cluster),"SP110"),
                             str_which(rownames(sce_cluster),"PTPRCAP"),
                             str_which(rownames(sce_cluster),"IGHMBP2"),
                             str_which(rownames(sce_cluster),"CD3G"),
                             str_which(rownames(sce_cluster),"HEATR3"),
                             str_which(rownames(sce_cluster),"LGALS3"),
                             str_which(rownames(sce_cluster),"ATP6V0E2"))]
for (i in seq_len(4)) {
  cluster3 <- rownames(sce_cluster[mcols(sce_cluster)$cluster == i])
  # cluster3 <- rownames(sce_cluster)
  ar_cluster3 <- res[,cluster3] %>% as.data.frame() %>% `colnames<-`(cluster3)
  part_cluster3 <- sapply(seq_len(ncol(ar_cluster3)),
                          function(j)match(ar_cluster3[,j],unique(ar_cluster3[,j]))) %>% as.data.frame() %>% `colnames<-`(cluster3)
  ar_cluster3$time <- cell_meta_unique
  part_cluster3$time <- cell_meta_unique
  ar_long <- melt(ar_cluster3, id = "time", value.name = "ar")    # Convert data to long format
  part_long <- melt(part_cluster3, id = "time", value.name = "partition")
  data_long <- merge(ar_long,part_long)
  data_long$partition <- factor(data_long$partition)
  data_long$time <- factor(data_long$time)
  data_long$variable <- sub("\\-protein.*", "", data_long$variable)
  data_long$variable <- sub("\\-antisense.*", "", data_long$variable)
  data_long$variable <- sub("\\-lincRNA.*", "", data_long$variable)
  # data_long <- data_long[!duplicated(data_long[,c("time","variable")]),]
  data[[i]]<-data_long
}

g1 <- ggplot(data[[1]], aes(x = time,y = ar)) +
  geom_point(aes(color = partition)) +
  geom_step(group = 1) +
  theme_minimal()+
  ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  labs(x = "Time points", y = "allelic ratio") +
  facet_wrap(~variable,nrow = 2)
g1

g2 <- ggplot(data[[2]], aes(x = time,y = ar)) +
  geom_point(aes(color = partition)) +
  geom_step(group = 1) +
  theme_minimal()+
  ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  labs(x = "Time points", y = "allelic ratio") +
  facet_wrap(~variable,nrow = 2)
g2

g3 <- ggplot(data[[3]], aes(x = time,y = ar)) +
  geom_point(aes(color = partition)) +
  geom_step(group = 1) +
  theme_minimal()+
  ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  labs(x = "Time points", y = "allelic ratio") +
  facet_wrap(~variable,nrow = 3)
g3

g4 <- ggplot(data[[4]], aes(x = time,y = ar)) +
  geom_point(aes(color = partition)) +
  geom_step(group = 1) +
  theme_minimal()+
  ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  labs(x = "Time points", y = "allelic ratio") +
  facet_wrap(~variable,nrow = 2)
g4
jpeg(file="C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/Github/airpartpaper/plot/dynamicmain3.jpeg",width = 1.65, height = 1.8,units = "in",res=360)
g3
dev.off()

## specific gene s value
sce_cluster<-sce[geneVar_name[9],]
mcols(sce_cluster)$cluster
sce_cluster <- fusedLasso(sce_cluster,formula = f,lambda = "is.bic",
                          model = "binomial",niter = 1,
                          genecluster = 223, ncores = 1,
                          adj.matrix = list(x = b2),
                          lambda.max = 8000)
sce_cluster <- allelicRatio(sce_cluster, formula = f)
s <- extractResult(sce_cluster,"svalue")

ggdraw() +
  draw_plot(g1, x = 0, y = 0, width = 0.29 , height = 1) +
  draw_plot(g2, x = 0.3, y = 0, width = 0.34 , height = 1) +
  draw_plot(g3, x = 0.65, y = 0.5, width = 0.34, height = 0.5) +
  draw_plot(g4, x = 0.65, y = 0, width = 0.34, height = 0.5) +
  draw_plot_label(label = c("A","B","C","D"), x= c(0, 0.3, 0.65,0.65), y = c(1,1,1,0.5), size = 10)

