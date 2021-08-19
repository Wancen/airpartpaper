
# script to test for dynamic ASE

library(aod)
library(lme4)

wd <- "/medpop/srlab/cd4tc_LS/rnaseq/both/waspGatkASE/TiDeASE/mixedEffects"
setwd(wd)

# arguments are timecourse ID (individual ID + replicate)
args <- commandArgs(trailingOnly = TRUE)
tc <- args[1]
print(tc)

# read ASE metrics for all sites and samples
 # both alleles seen per site (0 or 1)
bas <- read.table("C:/Users/wancenmu/Test Dropbox/Wancen Mu/ASEsummaries_cd4tcLS/summary_both_alleles_seen_swappedT4T8.txt", header = T, stringsAsFactors = F, row.names = 1)
 # reference allele counts per site
ref <- read.table("C:/Users/wancenmu/Test Dropbox/Wancen Mu/ASEsummaries_cd4tcLS/summary_refCount_swappedT4T8.txt", header = T, stringsAsFactors = F, row.names = 1)
  # alternative allele counts per site
alt <- read.table("C:/Users/wancenmu/Test Dropbox/Wancen Mu/ASEsummaries_cd4tcLS/summary_altCount_swappedT4T8.txt", header = T, stringsAsFactors = F, row.names = 1)
  # total coverage per site
cov <- read.table("/medpop/srlab/cd4tc_LS/rnaseq/both/waspGatkASE/summary_total_count_swappedT4T8.txt", header = T, stringsAsFactors = F, row.names = 1)
  # reference fraction per site
refRatio <- read.table("/medpop/srlab/cd4tc_LS/rnaseq/both/waspGatkASE/summary_refRatio_swappedT4T8.txt", header = T, stringsAsFactors = F, row.names = 1)

ref[is.na(ref)] <- 0
alt[is.na(alt)] <- 0
assay.list<-list(a1=ref,a2=alt)
# coldata<-data.frame(x=factor(celltype$x[order_of_development],levels=cell_meta_unique))
sce<-SingleCellExperiment(assays=assay.list)
sce <- preprocess(sce)
keep_feature2<-nexprs(counts(sce),byrow = TRUE, detection_limit = 20) >= ncol(sce)/2
table(keep_feature2)

# Select only those columns corresponding to the timecourse in question
indiv <- strsplit(tc, "_")[[1]][1]
rep <- strsplit(tc, "_")[[1]][2]
pos1 <- grep(indiv, colnames(refRatio))
pos2 <- grep(paste("_", rep, sep = ""), colnames(refRatio))
posIn <- pos1[which(pos1 %in% pos2)]
cnams <- colnames(refRatio)[posIn]
tpts <- sapply(cnams, function(sample){
			strsplit(sample, "_")[[1]][2]
		})
colsIn <- cnams[order(as.numeric(tpts), decreasing = F)]
sub.refRatio <- refRatio[,colsIn]
sub.cov <- cov[,colsIn]
sub.ref <- ref[,colsIn]
sub.alt <- alt[,colsIn]
sub.bas <- bas[,colsIn]

# Select sites with at least 10 reads in at least 4 samples (the mrs20 filter is further down)
mrs10 <- apply(sub.cov, 1, function(snp){
			length(which(snp >= 10))
  		})

posOK <- which(mrs10 >= 4)

rr.all <- sub.refRatio[posOK,]
cv.all <- sub.cov[posOK,]
rc.all <- sub.ref[posOK,]
ac.all <- sub.alt[posOK,]
bs.all <- sub.bas[posOK,]

# Select sites with both alleles seen in at least one time point (BASall filter further down)
bs <- apply(bs.all, 1, function(snp){
			length(which(snp == 1))
  		})
n <- length(colsIn) - 4
if(n < 1){
	print(paste("n = ", n, ". Please discard timecourse: ", tc, sep = ""))
    n <- 1
}

rr.all.bs <- rr.all[which(bs >= 1),]
cv.all.bs <- cv.all[which(bs >= 1),]
rc.all.bs <- rc.all[which(bs >= 1),]
ac.all.bs <- ac.all[which(bs >= 1),]

# Collect time info from colnames (for random effect t is overwritten to be 1,2,3,4...)
t <- as.numeric(sapply(colnames(rr.all.bs), function(sample){ strsplit(sample, "_")[[1]][2]}))

#################################
#### RUN LOGISTIC REGRESSION ####
#################################

lmt <- matrix(ncol = 41)
colnames(lmt) <- c("SNP",
                   "intercept", "intercept.stderr", "intercept.pval",
                   "timeCoef.randPoly", "timeStdErr.randPoly","timeP.randPoly",
                   "tSqCoef.randPoly", "tSqStdErr.randPoly", "tSqP.randPoly",
                   "Pall.randPoly",
                   "flagConv.rp", "flagRescale.rp", "flagOther.rp"
                   )
cntNotTested <- 0

for(i in 1:nrow(rc.all.bs)){
  #print(paste("Processing SNP", i))
  tot <- as.numeric(cv.all.bs[i,])

  posCov <- which(tot >= 20) #changed here to require mrs of 20

  # also require both alleles seen in all
  refC <- as.numeric(rc.all.bs[i,])
  altC <- as.numeric(ac.all.bs[i,])
  posBas <- which(refC > 0 & altC > 0)

  posIn <- posCov[which(posCov %in% posBas)]

  if(length(posIn) < 4){
    #print(paste("NOT TESTING SNP BC LESS THAN 4 TMPTS AFTER FILTERS, LINE:", i))
    cntNotTested <- cntNotTested + 1
  }
  else{
    #print(i)
    #print(paste("Positions to be included are", paste(posIn, collapse = " ")))
    tab <- rbind(rc.all.bs[i,posIn], ac.all.bs[i,posIn])

    ###############

    ##### logistic regression models

    # set time to go from 1 to 8 (for SNP2 you can change to 1:5)
    sub.t <- scale(1:ncol(tab))

    # format data to run logistic regression model
    x <- NULL
    y <- NULL
    for(j in 1:length(sub.t)){
      if(length(which(is.na(tab[,j]))) == 0){
        x <- c(x, rep(sub.t[j], (tab[1,j] + tab[2,j])))
        y <- c(y, c(rep(1, tab[1,j]), rep(0, tab[2,j])))
      }
    }

    mydata <- as.data.frame(cbind(x,y))
    colnames(mydata) <- c("time", "refAllele")
    mydata$timepoint <- as.factor(mydata$time)
    mydata$tsq <- mydata$time^2
    #head(mydata)
    mydata$refAllele <- as.factor(mydata$refAllele)

    # check if intercept significantly different from 0
    lm0 <- glm(refAllele ~ 1, data = mydata, family = "binomial")
    #summary(lm0)
    intercept <- coef(summary(lm0))[1,1]
    intercept.stderr <- coef(summary(lm0))[1,2]
    intercept.pval <- coef(summary(lm0))[1,4] #R uses standard normal to test if intercept diff from 0, 2*pnorm(-1*abs(estimate/std.error))

    # null model with intercept and sample as random effect
    lm1 <- glmer(refAllele ~ 1 + (1 | timepoint), data = mydata, family = "binomial")

    # alternative model adding time and time squared terms
    lm2.rand <- glmer(refAllele ~ 1 + time + tsq + (1 | timepoint), data = mydata, family = "binomial")

    timeCoefs.rand <- coef(summary(lm2.rand))[2,c(1,2,4)]
    timeSqCoefs.rand <- coef(summary(lm2.rand))[3,c(1,2,4)]

    # likelihood ratio test to test for significance of alternative model
    anv.rand <- anova(lm2.rand, lm1)
    p.rand <- anv.rand$P[2]

    # record any warnings
    s <- summary(lm2.rand)
    warn <- s$optinfo$conv$lme4$messages
    flagConv.rp <- 0 #
    flagRescale.rp <- 0
    flagOther.rp <- 0
    if(length(warn) > 0){
      #print(paste("i =",i,":", warn))
      if(length(grep("failed to converge", warn)) > 0){
        flagConv.rp <- 1
      }
      else{
        if(length(grep("Rescale variables", warn)) > 0){
          flagRescale.rp <- 1
        }
        else{
          flagOther.rp <- 1
        }
      }
    }

    res <- c(rownames(rc.all.bs)[i],
             intercept, intercept.stderr, intercept.pval,
             timeCoefs.rand,
             timeSqCoefs.rand,
             p.rand,
             flagConv.rp, flagRescale.rp, flagOther.rp,
             )

    #################

    if(nrow(mydata) > 0){
      lmt <- rbind(lmt, res)
    }
    else{
      lmt <- rbind(lmt, c(rownames(rc.all.bs)[i], rep(NA, (ncol(lmt) - 1))))
      print("Something went wrong...")
    }
  }
}
lmt <- lmt[-1,]
write.table(lmt, file = paste("results_dynASEtest_", tc, sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

warnings()
print(paste("Did not test", cntNotTested, "SNPs because less than 4 timepoints with mrs20 and basAll, for timecourse", tc))
