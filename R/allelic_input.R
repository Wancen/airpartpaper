cell_cycle <- data.frame(read.table("./data/cell_cycle_annotation.csv", header=TRUE, sep=",", row.names=1))

mESC_c57 <- data.frame(read.table("./data/SS3_c57_UMIs_mESC.csv", header=TRUE, sep=",", row.names=1))
mESC_cast <- data.frame(read.table("./data/SS3_cast_UMIs_mESC.csv", header=TRUE, sep=",", row.names=1))

fib_c57 <- data.frame(read.table("./data/SS3_c57_UMIs_concat.csv", header=TRUE, sep=",", row.names=1))
fib_cast <- data.frame(read.table("./data/SS3_cast_UMIs_concat.csv", header=TRUE, sep=",", row.names=1))

c57<-merge(mESC_c57,fib_c57,by="row.names")
rownames(c57)<-c57[,1]
c57<-c57[,-1]
cast<-merge(mESC_cast,fib_cast,by="row.names")
rownames(cast)<-cast[,1]
cast<-cast[,-1]

mESC_anno<-data.frame(x=rep("mESC",ncol(mESC_c57)),row.names = colnames(mESC_c57))
celltype<-rbind(mESC_anno,cell_cycle)

save(c57,cast,celltype,file = "./data/larsson.rda")


