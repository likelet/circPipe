library(dplyr)
library(ggplot2)
library(readr)

options(stringsAsFactors=F)
setwd("./")
args <-commandArgs(T)
source(args[1])
#source("/home/wqj/code/circPipe-master/R_function.R")
###########mean expression of circRNA and co-linear
##########data loading
#args[2]<-"/data2/wqj/singleend/gencode.rsem.fpkm_m6Astatus_11_29.mat"
linear_FPKM0 =  read_delim(args[2],delim = "\t")
#读入circRNA的表达矩阵
#args[3]<-"/data2/wqj/singleend/find_circ_merge.matrix"
circ_RPM_mat0 = read_delim(args[3],delim = "\t")
circ_RPM_mat0 <- circ_RPM_mat0[!duplicated(circ_RPM_mat0$id), ]

#读入注释矩阵
#args[4]<-"/data2/wqj/singleend/find_circ_annotation_annote.txt"
circ_feature=read.delim(args[4],header = FALSE,sep="\t",fill = T,col.names = c("chr","start","end","circ_type","exon_number","strand","strand2","ensemble_id","symbol","transcript","gene_feature","type1","type2","m6Astatus"))
circ_name = paste(circ_feature$chr,circ_feature$start,circ_feature$end,circ_feature$strand,sep = "_")
circ_feature$id = circ_name
circ_feature <- circ_feature[!duplicated(circ_feature$id), ]

#合并注释数据框和表达矩阵
circ_input_RPM_feature = inner_join(circ_feature,circ_RPM_mat0,by="id") #####需要测试

#求出每个circRNA的平均值
circ_RPM_mat<-circ_input_RPM_feature[,-(1:15)]
circ_RPM_mean = data.frame(gene = circ_input_RPM_feature$ensemble_id,
                           circ_name=circ_input_RPM_feature$id,mean_circ=rowMeans(circ_RPM_mat))
linear_FPKM<-linear_FPKM0[,-1]
linear_FPKM<-linear_FPKM[,-ncol(linear_FPKM)]
linear_FPKM_mean = data.frame(gene=as.character(linear_FPKM0$id),mean_linear=rowMeans(linear_FPKM))
circ_linear_mean = inner_join(circ_RPM_mean,linear_FPKM_mean,by="gene")

###############cor: circ and co-linear with mean expression level

# cor.test(circ_linear_mean$mean_circ,circ_linear_mean$mean_linear,method = "kendall")
# cormethod="spearman"  #"pearson", "kendall", "spearman"
# pcutoff=0.05
# 
# for (i in c("pearson", "kendall","spearman")) {
#   print(i)
#   print(cor.test(circ_linear_mean_over1$mean_circ,circ_linear_mean_over1$mean_linear,method = i))
#   #plot 
#   CoxExpPoint = Cor_plot(CoxExpPlotData = circ_linear_mean,gene1 =  "mean_circ",gene2 = "mean_linear", cormethod= i)+
#     ggtitle("Correlation of expression level of circRNAs and co-linear")
#   print(CoxExpPoint)
# }

CoxExpPoint1 = Cor_plot(CoxExpPlotData = circ_linear_mean,gene1 =  "mean_circ",gene2 = "mean_linear", cormethod= "pearson")+
  ggtitle("Correlation of expression level of circRNAs and co-linear")
png("pearson.png", width=2400, height=1200, res = 300)
plot(CoxExpPoint1)
dev.off()

CoxExpPoint2 = Cor_plot(CoxExpPlotData = circ_linear_mean,gene1 =  "mean_circ",gene2 = "mean_linear", cormethod= "kendall")+
  ggtitle("Correlation of expression level of circRNAs and co-linear")
png("kendall.png", width=2400, height=1200, res = 300)
plot(CoxExpPoint2)
dev.off()

CoxExpPoint3 = Cor_plot(CoxExpPlotData = circ_linear_mean,gene1 =  "mean_circ",gene2 = "mean_linear", cormethod= "spearman")+
  ggtitle("Correlation of expression level of circRNAs and co-linear")
png("spearman.png", width=2400, height=1200, res = 300)
plot(CoxExpPoint3)
dev.off()

pdf("correlation.pdf")
plot(CoxExpPoint1)
plot(CoxExpPoint2)
plot(CoxExpPoint3)
dev.off()