library(airpart)
library(tidyverse)
library(SingleCellExperiment)
# read ASE metrics for all sites and samples
# both alleles seen per site (0 or 1)
# bas <- read.table("C:/Users/wancenmu/Test Dropbox/Wancen Mu/ASEsummaries_cd4tcLS/summary_both_alleles_seen_swappedT4T8.txt", header = T, stringsAsFactors = F, row.names = 1)
# reference allele counts per site
ref <- read.table("C:/Users/wancenmu/Test Dropbox/Wancen Mu/ASEsummaries_cd4tcLS/summary_refCount_swappedT4T8.txt", header = T, stringsAsFactors = F, row.names = 1)
# alternative allele counts per site
alt <- read.table("C:/Users/wancenmu/Test Dropbox/Wancen Mu/ASEsummaries_cd4tcLS/summary_altCount_swappedT4T8.txt", header = T, stringsAsFactors = F, row.names = 1)
# gene snp mapping file
snp2gene <- read.table("C:/Users/wancenmu/Test Dropbox/Wancen Mu/ASEsummaries_cd4tcLS/gene_snp_mapping_file.txt", stringsAsFactors = F) %>% `colnames<-`(c("chr","pos","ID", "intergenic"))

sample <- colnames(ref)
id<-do.call(rbind,str_split(sample,"_")) %>% as.data.frame()
ref[is.na(ref)] <- 0
alt[is.na(alt)] <- 0
total <- ref + alt
cell_meta_unique <- c(0,2,4,8,12,24,48,72)
order_of_development <- order(match(id$V2,cell_meta_unique))

sigsnp<- rownames(ref)

snp2gene_sig <- snp2gene[which(snp2gene$ID %in% unique(sigsnp)),]

geneGuess <- NULL
for(i in 1:length(sigsnp)){
    snp <- sigsnp[i]
    gene <- snp2gene_sig[which(snp2gene_sig[,3] %in% snp),4]
    if(length(gene) == 0){
      print("ERROR: SNP not in snp2gene file")
      print(snp)
      gene <- "Not found"
      geneID <- "Intergenic"
      geneName <- "Intergenic"
    }else{
      #sp <- strsplit(gene, ":")[[1]]
      sp <- strsplit(gene, ";")[[1]]
      exon <- sp[grep("-exon", sp)]
      if(length(exon) > 0){
        pc <- exon[grep("protein", exon)]
        if(length(pc) > 0){
          parts <- strsplit(pc[1], ":")[[1]]
          geneID <- parts[1]
          geneName <- parts[4]
        }
        else{
          lr <- exon[grep("lincRNA", exon)]
          if(length(lr) > 0){
            parts <- strsplit(lr[1], ":")[[1]]
            geneID <- parts[1]
            geneName <- parts[4]
          }
          else{
            parts <- strsplit(exon[1], ":")[[1]]
            geneID <- parts[1]
            geneName <- parts[4]
          }
        }
      }
      else{
        intron <- sp[grep("-intron", sp)]
        if(length(intron) > 0){
          pc <- intron[grep("protein", intron)]
          if(length(pc) > 0){
            parts <- strsplit(pc[1], ":")[[1]]
            geneID <- parts[1]
            geneName <- parts[4]
          }
          else{
            lr <- intron[grep("lincRNA", intron)]
            if(length(lr) > 0){
              parts <- strsplit(lr[1], ":")[[1]]
              geneID <- parts[1]
              geneName <- parts[4]
            }
            else{
              parts <- strsplit(intron[1], ":")[[1]]
              geneID <- parts[1]
              geneName <- parts[4]
            }
          }
        }
        else{
          parts <- strsplit(sp[1], ":")[[1]]
          geneID <- parts[1]
          geneName <- parts[4]
        }
      }
    }
    geneGuess <- rbind(geneGuess, c(snp,geneID, geneName))
  }

index <- split(1:nrow(geneGuess), geneGuess[,3])

library(pbapply)
snpindex <- pbsapply(1:length(index), function(i){
  poi <- which.max(rowSums(total[index[[i]],]))
  poi
})
snpname <- names(snpindex)
a1_new <- pbsapply(1:length(index), function(i){
  ref[index[[i]][snpindex[i]],] %>% as.numeric()
})
colnames(a1_new)<-names(index)
rownames(a1_new)<-sample
a2_new <- pbsapply(1:length(index), function(i){
  alt[index[[i]][snpindex[i]],] %>% as.numeric()
})
colnames(a2_new)<-names(index)
rownames(a2_new)<-sample
save(geneGuess,index,snpindex,snpname,a1_new,a2_new,id,cell_meta_unique,order_of_development,file="C:/Users/wancenmu/Test Dropbox/Wancen Mu/ASEsummaries_cd4tcLS/processed.RData")


