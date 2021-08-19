library(patchwork)
patchwork <- p2 + p1 + inset_element(p1.2, left = 0.5, bottom = 0.6, right = 1, top = 1)
patchwork + plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size = 10))

library(cowplot)
fp1 <- grid.grabExpr(p1)
ggdraw() +
  draw_plot(p1, x = 0, y = 0, width = 1 , height = 0.5) +
  draw_plot(p, x = 0, y = 0.5, width = 1 , height = 0.5) +
  draw_plot_label(label = c("A","B"), x= c(0, 0.5), y = c(1,1), size = 10)

l11 <- lapply(list.files(path="C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1/data/csv", full=TRUE), function(f) {
  read.csv(f, header=FALSE)
})
dat11 <- do.call(rbind, l11)
colnames(dat11) <- c("type","ARI","time","n","cnt")
dat11$n <- factor(dat11$n, levels = c(5, 10, 20,40))
dat11$cnt <- factor(dat11$cnt, levels = c(5, 10, 20))
dat11$ARI<-as.numeric(dat11$ARI)
dat11$time<-as.numeric(dat11$time)
dat11<-dat11%>% filter(type%in%c("fl_binomial","fl_gaussian","WilcoxExt")&n %in%c(5,10,20,40))
dat11$type <- factor(dat11$type,ordered = F, levels = c("fl_binomial","fl_gaussian","wilcoxExt"),
                      labels  = c("airpart.bin","airpart.gau","airpart.np","scdali"))

p2 <- ggplot(dat11 , aes(x=ARI,fill=type)) +
  geom_bar(stat="count", position="dodge", width=0.1) +
  scale_x_continuous(limits=c(-0.05, 1.08), breaks=c(0, 0.25, 0.5, 0.75,1))+
  facet_grid(n ~cnt, labeller = label_both) +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired")+
  theme(legend.position = "none",
    panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent')) +
  labs( y = "count")
p2


p3 <- ggplot(dat11, aes(x=type, y=time,fill=type)) +
  geom_boxplot(outlier.color=NA) +
  geom_jitter(width=.1,size=0.1) +
  facet_grid(n ~ cnt, labeller = label_both)+ylim(0,50) +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired")+
  theme(legend.position = "none",panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        axis.text.x = element_blank()) +
  labs(x = "methods", y = "time(s)")
p3

patchwork <- (p2 + p3 )/ p4
patchwork + plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size = 10))

sce <- makeSimulatedData(ncl=1,p.vec = c(0.25,0.25,0.26,0.26),ngenecl = 10)
sce <-preprocess(sce)
sce<-geneCluster(sce,G=1)
sce_sub <- fusedLasso(sce,genecluster = 1)
metadata(sce_sub)
sce_sub <- allelicRatio(sce_sub)

p4+p5

sce_sub <- fusedLasso(sce,model = "binomial", genecluster = 1,ncores = 2, se.rule.mult = 0)
metadata(sce_sub)$partition
