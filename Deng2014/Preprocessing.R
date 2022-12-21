library(tidyverse)
library(scater)
library(pheatmap)
library(factoextra)
library(pbapply)
library(mclust)
library(parallel)
nct=10

## Preprocess ############################################
#C57BL/6J(mother) x CAST/EiJ(father)
B6<-read_delim("./data/B6.txt",delim =" ", skip = 1,col_names = F)
CAST<-read_delim("./data/CAST.txt",skip = 1,delim=" ",col_names = F)
genename<-read_delim("./data/genename.txt",skip = 1,delim=" ",col_names = F)[,1]
sampleinfo<-scan("./data/sample.txt", character(), quote = "")
celltype<-sub("_.*", "", sampleinfo)


## Check categories of cells, remove unuseful cells and outlier cells indicated in Deng's Supplyment###################
ratio1<-(CAST)/(CAST+B6)

celltype[grep("zy",celltype)]="zy";
sampleinfo[grep("8cell_2pooled", sampleinfo)]="8cell_nd";
sampleinfo[grep("8cell_split", sampleinfo)]="8cell_nd";
sampleinfo[grep("16cell_2pooled", sampleinfo)]="16cell_nd";
sampleinfo[grep("16cell_split", sampleinfo)]="16cell_nd";
sampleinfo[grep("8cell_1", sampleinfo)]="8cell_nd";
sampleinfo[grep("8cell_5", sampleinfo)]="8cell_nd";
sampleinfo[grep("lateblast_2", sampleinfo)]="lateblast_nd";
midblast<-colMeans(ratio1[,grep("midblast", sampleinfo)],na.rm = T)
hist(midblast) #one extremely high midblast
sampleinfo[order(midblast,decreasing = T)[1]]="midblast_nd";
eicell<-colMeans(ratio1[,grep("8cell", sampleinfo)],na.rm = T)
hist(eicell) #one extremely low 8cell
sampleinfo[order(eicell)[1]]="8cell_nd";
sixtcell<-colMeans(ratio1[,grep("16cell", sampleinfo)],na.rm = T)
hist(sixtcell) #one extremely high 16cell
sampleinfo[order(sixtcell,decreasing = T)[1]]="16cell_nd";
indices_not_reqd <- which(celltype=="BXC"   | celltype=="C57twocell" | celltype=="fibroblast" | sampleinfo =="8cell_nd"| sampleinfo =="lateblast_nd" |sampleinfo =="midblast_nd" | sampleinfo == "16cell_nd")
celltype<-celltype[-indices_not_reqd]
sampleinfo<-sampleinfo[-indices_not_reqd]
table(celltype)

#return the position where have duplicate gene name
genedui<-which(duplicated(genename) | duplicated(genename, fromLast = TRUE))

## remove unrelated genes and cells and order count matrix ####
b6<-B6[-genedui,-indices_not_reqd] %>% as.matrix()
rownames(b6)<-unname(as.character(genename$X1[-genedui]))
colnames(b6)<-sampleinfo
cast<-CAST[-genedui,-indices_not_reqd] %>% as.matrix()
row.names(cast)<-unlist(genename[-genedui,1])
colnames(cast)<-sampleinfo

cell_meta_unique <- c("zy","early2cell","mid2cell","late2cell","4cell","8cell","16cell","earlyblast","midblast","lateblast") ;

save(b6,cast,celltype,file = "./data/Deng.rda")

## Change the order based on development order ############################
order_of_development <- order(match(celltype,cell_meta_unique))
sampleinfo<-sampleinfo[order_of_development]
celltype<- celltype[order_of_development]



b6_uq<-b6_uq[,order_of_development]
colnames(b6_uq)<-sampleinfo

cast_uq<-cast_uq[,order_of_development]
colnames(cast_uq)<-sampleinfo
total<-b6_uq+cast_uq


## Define embryo level for futher research ##############
cell_embryo<-sub("-.*", "", sampleinfo)
cell_embryo[grep("zy",cell_embryo)]="zy";
cell_embryo_unique<-unique(cell_embryo)
## create annotation df for heatmap ##########################
anno_df = data.frame(CellType = factor(celltype,levels = cell_meta_unique),row.names = sampleinfo) #cell type level
anno_df2 = data.frame(CellType = factor(cell_embryo,levels = cell_embryo_unique),row.names = sampleinfo) #embryo level

## specific to Deng's data ####
## define cell types except 3 mainly expressed cell types for filtering ###############
celltype_short<-celltype[cellQC]
indices_specific<-which(celltype_short%in%c("early2cell","late2cell","zy"))
celltype_short[-indices_specific]<-"others"

## filter out genes mainly expressed in three cell types #################
celltype_prop<- pbsapply(1:nrow(total_fil),function(i){tapply(unlist(total_fil[i,]),as.factor(celltype_short), FUN=sum)})
other_prop<-celltype_prop["others",]/colSums(celltype_prop)
h2<-hist(other_prop,main="4597 genes",
         xlab="prop of other cell types")
text(h2$mids,h2$counts,labels=h2$counts, adj=c(0.5, -0.5))
abline(v=0.8,col="red")
check2<-which(other_prop>0.8)

save(cast_uq,keep_feature, keep_feature2,check2,total,cellQC,cell_meta_unique,celltype,anno_df,file = "./data/postprocess.rda")
