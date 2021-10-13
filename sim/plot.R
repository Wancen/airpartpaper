## main1
library(patchwork)
library(tidyr)
library(dplyr)
library(ggsci)
load("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/scdali/default_n40.rda")
nct <-8
n=40
true <- (rowData(sce)[,1:nct]) %>% as.data.frame()
x_start<-seq(1,n*nct,n)
ar_scdali<-scdali[1:320,]
sd_scdali<-scdali[321:640,]
mu<-ar_scdali[x_start,]
sum(is.na(colSums(mu)))
sd<-sd_scdali[x_start,]
diff_scdali <- true-t(mu)
rmse_scdali<-sqrt(rowMeans(diff_scdali^2))

fl_mean <- fl[1:nct,]
sum(is.na(colSums(fl_mean)))
upper <-fl[(nct+1):(2*nct),]
lower <-fl[(2*nct+1):(3*nct),]
diff_fl <- true-t(fl_mean)
rmse_fl<-sqrt(rowMeans(diff_fl^2))


fl_binomial_mean <- fl_binomial[1:nct,]
sum(is.na(colSums(fl_binomial_mean)))
upper_binomial <-fl_binomial[(nct+1):(2*nct),]
lower_binomial <-fl_binomial[(2*nct+1):(3*nct),]
diff_fl_binom <- true-t(fl_binomial_mean)
rmse_fl_binom<-sqrt(rowMeans(diff_fl_binom^2))
##
fl_gaussian_mean <- fl_gaussian[1:nct,]
upper_gaussian <-fl_gaussian[(nct+1):(2*nct),]
lower_gaussian <-fl_gaussian[(2*nct+1):(3*nct),]
diff_fl_gaussian <- true-t(fl_gaussian_mean)
rmse_fl_gaussian<-sqrt(rowMeans(diff_fl_gaussian^2))

##
wilcoxon_mean <- wilcoxon[1:nct,]
upper_wilcoxon <-wilcoxon[(nct+1):(2*nct),]
lower_wilcoxon <-wilcoxon[(2*nct+1):(3*nct),]
diff_wilcoxon <- true-t(wilcoxon_mean)
rmse_wilcoxon<-sqrt(rowMeans(diff_wilcoxon^2))


dat <- data.frame(rmse_fl_binom=rmse_fl_binom,
                  rmse_fl_gaussian=rmse_fl_gaussian,
                  rmse_wilcoxon=rmse_wilcoxon,
                  rmse_scdali=rmse_scdali,
                  rmse_fl=rmse_fl,
                  # rmse_scdali=rmse_scdali,
                  type=factor(c(rep(round(ardiff,3),each=400))))

# dat <- data.frame(rmse_fl=rmse_fl,
#                   rmse_scdali=rmse_scdali,
#                   type=factor(c(rep(round(ardiff,3),each=400))))
dat2<-dat %>% gather(key=method,value=rmse,rmse_fl_binom:rmse_fl)
dat2  %>% group_by(method,type) %>% summarise(fail=sum(is.na(rmse))/4)
# dat2$method <- factor(dat2$method,ordered = F, levels = c("rmse_fl_binom","rmse_scdali"),
#                       labels  = c("airpart","scdali"))
dat2$method <- factor(dat2$method,ordered = T, levels = c("rmse_fl_binom","rmse_fl_gaussian","rmse_wilcoxon","rmse_fl","rmse_scdali"),
                      labels  = c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"))

mu_0.3 <- t(mu[,801:1200]) %>% `colnames<-`(paste("scdali",paste0("ct",1:8)))
binomal_0.3 <- t(fl_binomial_mean[,801:1200])%>% `colnames<-`(paste("binomial",paste0("ct",1:8)))
guassian_0.3 <- t(fl_gaussian_mean[,801:1200])%>% `colnames<-`(paste("gaussian",paste0("ct",1:8)))
wilcoxon_0.3 <- t(wilcoxon_mean[,801:1200])%>% `colnames<-`(paste("wilcoxon",paste0("ct",1:8)))
nogroup_0.3 <- t(fl_mean[,801:1200])%>% `colnames<-`(paste("nogroup",paste0("ct",1:8)))

dat_0.3 <- data.frame(cbind(binomal_0.3,guassian_0.3,mu_0.3,nogroup_0.3,wilcoxon_0.3))
dat2_0.3<-dat_0.3 %>% gather(key=type,value=ar,binomial.ct1:wilcoxon.ct8)
dat2_0.3$type <- factor(dat2_0.3$type,ordered = F, levels = paste(rep(c("binomial","gaussian","wilcoxon","nogroup","scdali"),each=nct),rep(paste0("ct",1:nct),4),sep = "."),
                        labels  = paste(rep(c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"),each=nct),rep(paste0("ct",1:nct),4),sep = "."))
dat2_0.3$method <- rep(c(rep("airpart.bin",nct),rep("airpart.gau",nct),rep("scdali",nct),rep("airpart.nogroup",nct),rep("airpart.np",nct)),each=400)
dat2_0.3$method <- factor(dat2_0.3$method,ordered = T, levels = c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"),
                          labels  = c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"))
# dat2_0.3$true <- rep(rep(true[1601,],each=400),4)
dat_true <- data.frame(type=levels(dat2_0.3$type),true=rep(as.numeric(true[801,]),5))

p1<-ggplot(dat2,mapping = aes(x=type,y=rmse,fill=method)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_igv()+
  # scale_fill_jco()+
  theme(legend.position = "right",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent')) +
  xlab("Allelic ratio difference")+
  ylab("RMSE")
# scale_fill_discrete(name="methods", labels=c("airpart.bin","airpart.gau","airpart.np","scdali"))
p1
p1.2<- ggplot() +
  geom_boxplot(data=dat2_0.3[6401:12800,],mapping = aes(x=type,y=ar,fill=method))+
  geom_point(dat_true[25:40,],mapping = aes(x=type,y=true,col="red")) +
  theme_classic() +
  scale_fill_manual(
    values = pal_igv("default")(5)[4:5])+
  theme(legend.position = "none",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  xlab("Cell type") +
  ylab("Allelic Ratio")+
  # facet_wrap(~ar) +
  scale_x_discrete(breaks=paste(rep(c("airpart.nogroup","scdali"),each=nct),rep(paste0("ct",1:nct),2),sep = "."),
                   labels=rep(c(paste0("ct",1:nct)),2))
p1.2
# pcompare <- p1 + inset_element(p1.2, left = 0.4, bottom = 0.55, right = 0.995, top = 0.995)


l5 <- lapply(list.files(path="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/scdali/csv", full=TRUE), function(f) {
  read.csv(f, header=FALSE)
})
dat5 <- do.call(rbind, l5)
colnames(dat5) <- c("type","ARI","time","g","cnt")
dat5.2<-dat5%>% filter(cnt%in%c(5,10,20)&type%in%c("fl_binomial","fl_gaussian","WilcoxExt")&g%in%c(5,10,20))
str(dat5.2)
dat5.2$type <- ifelse(dat5.2$type=="fl_binomial","airpart.bin",ifelse(dat5.2$type=="fl_gaussian","airpart.gau","airpart.np"))

colnames(dat5.2) <- c("method","ARI","time","g","cnt")
dat5.2$method<-factor(dat5.2$method)
# p2<- ggplot(dat5.2, aes(x=ARI,fill=method)) +
p2<-ggplot(dat5.2, aes(x=method,y=ARI,fill=method)) +
  geom_boxplot()+
  # geom_bar(aes(y = ..prop..), position="dodge", width=0.15) +scale_x_continuous(limits=c(0, 1.08), breaks=c(0, 0.25, 0.5, 0.75,1))+
  # scale_y_continuous(labels = scales::percent,limits=c(0, 1)) +
  facet_grid(g ~ cnt, labeller = label_both)+
  # scale_fill_brewer(palette = "Paired")+
  scale_fill_igv()+
  theme_classic() +
  ylab("ARI")+
  theme(legend.position = "none",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=8))

p3<-ggplot(dat5.2, aes(x=method,y=time,fill=method)) +
  geom_boxplot(outlier.color=NA)+
  geom_jitter(width=.1,size=0.2)+
  facet_grid(g ~ cnt, labeller = label_both)+
  # scale_fill_brewer(palette = "Paired")+
  scale_fill_igv()+
  theme_classic() +
  ylab("time(s)")+ ylim(0,50)+
  theme(legend.position = "none",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=8))

jpeg(file="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/plots/main1.jpg",width = 12, height = 7,units = "in",res=450)
p2+p1/p1.2+ plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect')
dev.off()

## theta=3 in supp
l5 <- lapply(list.files(path="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/scdali/csv", full=TRUE), function(f) {
  read.csv(f, header=FALSE)
})
dat5 <- do.call(rbind, l5)
colnames(dat5) <- c("type","ARI","time","g","cnt")
dat5.2<-dat5%>% filter(cnt%in%c(5,10,20)&type%in%c("bin","gau","Wilcoxon")&g%in%c(200,500,1000))
str(dat5.2)
dat5.2$type <- ifelse(dat5.2$type=="bin","airpart.bin",ifelse(dat5.2$type=="gau","airpart.gau","airpart.np"))
dat5.2$g <- ifelse(dat5.2$g==200,5,ifelse(dat5.2$g==500,10,25))

colnames(dat5.2) <- c("method","ARI","time","g","cnt")
dat5.2$method<-factor(dat5.2$method)
# p2<- ggplot(dat5.2, aes(x=ARI,fill=method)) +
p2<-ggplot(dat5.2, aes(x=method,y=ARI,fill=method)) +
  geom_boxplot()+
  # geom_jitter(width=.1,size=0.2)+
  # geom_bar(aes(y = ..prop..), position="dodge", width=0.15) +scale_x_continuous(limits=c(0, 1.08), breaks=c(0, 0.25, 0.5, 0.75,1))+
  # scale_y_continuous(labels = scales::percent,limits=c(0, 1)) +
  facet_grid(g ~ cnt, labeller = label_both)+
  # scale_fill_brewer(palette = "Paired")+
  scale_fill_igv()+
  theme_classic() +
  ylab("ARI")+
  theme(legend.position = "none",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=8))


jpeg(file="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/plot/suppsim1.jpg",width = 12, height = 8,units = "in",res=450)
p3+p2+ plot_annotation(tag_levels = 'A')
dev.off()

## supp2 n=100
load("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/scdali/default_n100.rda")
nct <-8
n=100
true <- (rowData(sce)[,1:nct]) %>% as.data.frame()
x_start<-seq(1,n*nct,n)
ar_scdali<-scdali[1:(n*8),]
sd_scdali<-scdali[(n*8+1):(n*16),]
mu<-ar_scdali[x_start,]
sum(is.na(colSums(mu)))
sd<-sd_scdali[x_start,]
diff_scdali <- true-t(mu)
rmse_scdali<-sqrt(rowMeans(diff_scdali^2))

fl_mean <- fl[1:nct,]
sum(is.na(colSums(fl_mean)))
upper <-fl[(nct+1):(2*nct),]
lower <-fl[(2*nct+1):(3*nct),]
diff_fl <- true-t(fl_mean)
rmse_fl<-sqrt(rowMeans(diff_fl^2))

fl_binomial_mean <- fl_binomial[1:nct,]
sum(is.na(colSums(fl_binomial_mean)))
upper_binomial <-fl_binomial[(nct+1):(2*nct),]
lower_binomial <-fl_binomial[(2*nct+1):(3*nct),]
diff_fl_binom <- true-t(fl_binomial_mean)
rmse_fl_binom<-sqrt(rowMeans(diff_fl_binom^2))
##
fl_gaussian_mean <- fl_gaussian[1:nct,]
upper_gaussian <-fl_gaussian[(nct+1):(2*nct),]
lower_gaussian <-fl_gaussian[(2*nct+1):(3*nct),]
diff_fl_gaussian <- true-t(fl_gaussian_mean)
rmse_fl_gaussian<-sqrt(rowMeans(diff_fl_gaussian^2))

##
wilcoxon_mean <- wilcoxon[1:nct,]
upper_wilcoxon <-wilcoxon[(nct+1):(2*nct),]
lower_wilcoxon <-wilcoxon[(2*nct+1):(3*nct),]
diff_wilcoxon <- true-t(wilcoxon_mean)
rmse_wilcoxon<-sqrt(rowMeans(diff_wilcoxon^2))

dat <- data.frame(rmse_fl_binom=rmse_fl_binom,
                  rmse_fl_gaussian=rmse_fl_gaussian,
                  rmse_wilcoxon=rmse_wilcoxon,
                  rmse_scdali=rmse_scdali,
                  rmse_fl=rmse_fl,
                  # rmse_scdali=rmse_scdali,
                  type=factor(c(rep(round(ardiff,3),each=400))))

# dat <- data.frame(rmse_fl=rmse_fl,
#                   rmse_scdali=rmse_scdali,
#                   type=factor(c(rep(round(ardiff,3),each=400))))
dat2<-dat %>% gather(key=method,value=rmse,rmse_fl_binom:rmse_fl)
dat2  %>% group_by(method,type) %>% summarise(fail=sum(is.na(rmse))/4)
# dat2$method <- factor(dat2$method,ordered = F, levels = c("rmse_fl_binom","rmse_scdali"),
#                       labels  = c("airpart","scdali"))
dat2$method <- factor(dat2$method,ordered = T, levels = c("rmse_fl_binom","rmse_fl_gaussian","rmse_wilcoxon","rmse_fl","rmse_scdali"),
                      labels  = c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"))
p1<-ggplot(dat2,mapping = aes(x=type,y=rmse,fill=method)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_igv()+
  # scale_fill_jco()+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent')) +
  xlab("Allelic ratio difference")+
  ylab("RMSE")
# scale_fill_discrete(name="methods", labels=c("airpart.bin","airpart.gau","airpart.np","scdali"))
jpeg(file="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/plot/sim60.jpg",width = 12, height = 7,units = "in",res=450)
p1
dev.off()

## CI coverage
scdali_lower <-mu-qnorm(0.95)*sd
scdali_upper <- mu+qnorm(0.95)*sd
ngene = 1600
poi_scdali <- which(is.na(colSums(mu)))
ci_scdali <- matrix(0,nrow = ngene,ncol = nct) %>% `colnames<-`(paste("scdali",paste0("ct",1:nct)))
ci_scdali[which(true<t(scdali_upper) & true>t(scdali_lower))]<-1
ci_scdali[poi_scdali,]<-NA

poi_binomial <- which(is.na(colSums(lower_binomial)))
ci_binomial <- matrix(0,nrow = ngene,ncol = nct)%>% `colnames<-`(paste("binomial",paste0("ct",1:nct)))
ci_binomial[which(true<t(upper_binomial) & true>t(lower_binomial))]<-1
ci_binomial[poi_binomial,]<-NA

ci <- matrix(0,nrow = ngene,ncol = nct)%>% `colnames<-`(paste("nogroup",paste0("ct",1:nct)))
ci[which(true<t(upper) & true>t(lower))]<-1

ci_gaussian <- matrix(0,nrow = ngene,ncol = nct)%>% `colnames<-`(paste("gaussian",paste0("ct",1:nct)))
ci_gaussian[which(true<t(upper_gaussian) & true>t(lower_gaussian))]<-1

ci_wilcoxon <- matrix(0,nrow = ngene,ncol = nct)%>% `colnames<-`(paste("wilcoxon",paste0("ct",1:nct)))
ci_wilcoxon[which(true<t(upper_wilcoxon) & true>t(lower_wilcoxon))]<-1

datci <- data.frame(cbind(ci_scdali,ci,ci_binomial,ci_gaussian,ci_wilcoxon),
                    ar=factor(c(rep(round(ardiff,3),each=400))))
datci2<-datci %>% gather(key=type,value=coverage,scdali.ct1:wilcoxon.ct8)
datci2 <- datci2[!is.na(datci2$coverage),]
library(dplyr)
detach("package:plyranges", unload=TRUE)
datci2<-datci2%>% group_by(ar,type) %>% dplyr::summarise(coverage=sum(coverage,na.rm=T)/n())
datci2$type <- factor(datci2$type,ordered = F, levels = paste(rep(c("binomial","gaussian","wilcoxon","nogroup","scdali"),each=nct),rep(paste0("ct",1:nct),5),sep = "."),
                      labels  = paste(rep(c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"),each=nct),rep(paste0("ct",1:nct),5),sep = "."))
datci2$method <- rep(c(rep("airpart.bin",nct),rep("airpart.gau",nct),rep("airpart.nogroup",nct),rep("scdali",nct),rep("airpart.np",nct)),length(ardiff))

datci2$method <- factor(datci2$method,ordered = T, levels = c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"),
                      labels  = c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"))
library(ggplot2)
library(ggsci)
pci40<-ggplot(datci2,mapping = aes(x=type,y=coverage,fill=method)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_igv()+
  geom_hline(yintercept=0.95, linetype="dashed",
             color = "red", size=1)+
  theme(legend.position = "none",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  xlab("Allelic ratio difference")+
  facet_wrap(~ar) +
  scale_x_discrete(breaks=paste(rep(c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"),each=nct),rep(paste0("ct",1:nct),4),sep = "."),
                   labels=rep(c(paste0("ct",1:nct)),5))
pci40

pci100<-ggplot(datci2,mapping = aes(x=type,y=coverage,fill=method)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_igv()+
  geom_hline(yintercept=0.95, linetype="dashed",
             color = "red", size=1)+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=8)) +
  xlab("Allelic ratio difference")+
  facet_wrap(~ar)+
  scale_x_discrete(breaks=paste(rep(c("airpart.bin","airpart.gau","airpart.np","airpart.nogroup","scdali"),each=nct),rep(paste0("ct",1:nct),5),sep = "."),
                   labels=rep(c(paste0("ct",1:nct)),5))
jpeg(file="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/plot/suppsim2_new.jpg",width = 14, height = 10,units = "in",res=450)
p1+ pci40 / pci100 + plot_annotation(tag_levels = 'A')
dev.off()


## DAI test ##################
load("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/scdali/fp.rda")
nct <-6
n=40

true <- rep(c(FALSE,TRUE),each=400)
table(true,scdali_res[1:800])/400
table(true,fl_binomial[1:800])/400
table(true,fl_gaussian[1:800])/400
table(true,wilcoxon[1:800])/400
