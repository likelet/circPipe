library("readr")
library(dplyr)
#source("/home/wqj/code/circPipe-master/R_function.R")

options(stringsAsFactors=F)

#setwd("/home/wqj/code/circPipe-master")
library(Biobase)
library(edgeR)
########cut off init
p_value <- 0.05 # p value of FDR
lfc <- 0.58  #logFC
args <-commandArgs(T)

source(args[1])

########load data
countData1 <- as.data.frame(read_delim(args[2],delim = "\t"))
dim(countData1)
rownames(countData1) <- countData1$id
countData=countData1[,-1]

###########pheno data loading
colData <- get_pheno(x = colnames(countData),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal")
group1 = "Tumor"
group2 = "Normal"
label1 = "T"
label2 = "N"
colnames(countData)=sub("X","S",as.character(colnames(countData)))
countData <- countData[,as.character(rownames(colData))]
colData$Type <- as.factor(colData$Type)
type_level <- levels(as.factor(colData$Type))
comb <- combn(type_level,2)
Circ_norm_edgeR=cpm(countData)
#########filter low-abandance circRNA ; the step has been done in node1:filter_circ
 countData <- countData[which(rowSums(countData > 0) >= 2),] #a Count>0 in at least 2 samples
 dim(countData)
  ###norm expres_mat


##############edgeR
#source("./R_function.R")

#sharedCirc_edgeR_tmp <- edgeR_test(expre_mat = shared_circ,group_mat = colData,test_method = "LRT" )
if( sum(colData$Type==group1)< 2 || sum(colData$Type==group2)< 2 ){
  mean.mat = data.frame(group1_mean = rowMeans(Circ_norm_edgeR[,grep(label1,colnames(Circ_norm_edgeR)), drop = FALSE]),
                        group2_mean =rowMeans(Circ_norm_edgeR[,grep(label2,colnames(Circ_norm_edgeR)),drop = FALSE]),
                        all_mean = rowMeans(Circ_norm_edgeR))
  rownames(mean.mat) = rownames(Circ_norm_edgeR)
  mean.mat = subset(mean.mat,group1_mean > 0)
  
  sharedCirc_edgeR = data.frame(logFC=log2(mean.mat$group2_mean/mean.mat$group1_mean+1),
                        logCPM = log2(mean.mat$all_mean),
                    `F`= rep("_",nrow(mean.mat)),PValue = rep("_",nrow(mean.mat)),
                    FDR =rep(0,nrow(mean.mat)))
  rownames(sharedCirc_edgeR) = rownames(mean.mat)
}else {
  sharedCirc_edgeR <-edgeR_test(expre_mat = countData,group_mat = colData)
}

DE_sharedCirc_res = subset(sharedCirc_edgeR,FDR <= p_value & abs(logFC) >= lfc)
select <- DE_sharedCirc_res[order(DE_sharedCirc_res$logFC, decreasing = TRUE),]

rownames(Circ_norm_edgeR)=countData1$id
sharedCirc_DE_edgeR=Circ_norm_edgeR[rownames(select),]
#all DE
DE_list <- c(rownames(DE_sharedCirc_res))

################plot
#pdf(file = paste0("DE/",comb[2,1], "_vs_", comb[1,1], "_DE_10_17.edgeR.pdf"))
########Vocano plot
png(args[3], width=1800, height=1800, res = 300)
if(sum(colData$Type==group1)< 2 || sum(colData$Type==group2)< 2){
  print("cannot plot volcano plot.")
}else{
  volcano_plot(sharedCirc_edgeR,c(group1,group2))
}
dev.off()


########heatmap
png(args[4], width=1800, height=1800, res = 300)
heatmap_house(Circ_norm_edgeR[DE_list,],get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),title_hp = "Heatmap of Different expressed all DE circRNA")
dev.off()
png(args[5], width=1800, height=1800, res = 300)
heatmap_house(Circ_norm_edgeR[DE_list,],get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),
              title_hp = "Heatmap of Different expressed all DE circRNA",cluster_rule = "row")
dev.off()
png(args[6], width=1800, height=1800, res = 300)
heatmap_house(sharedCirc_DE_edgeR,get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),title_hp = "Heatmap of Different expressed shared_circRNA")
dev.off()
########PCA
png(args[7], width=1800, height=1800, res = 300)
PCA_plot(DE_list,Circ_norm_edgeR,colData,withtext=FALSE)
dev.off()
png(args[8], width=1800, height=1800, res = 300)
PCA_plot(DE_list,Circ_norm_edgeR,colData,withtext=TRUE)
dev.off()

pdf(args[9])
if(sum(colData$Type==group1)< 2 || sum(colData$Type==group2)< 2){
  print("cannot plot volcano plot.")
}else{
  volcano_plot(sharedCirc_edgeR,c(group1,group2))
}
heatmap_house(Circ_norm_edgeR[DE_list,],get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),title_hp = "Heatmap of Different expressed all DE circRNA")
heatmap_house(Circ_norm_edgeR[DE_list,],get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),
              title_hp = "Heatmap of Different expressed all DE circRNA",cluster_rule = "row")
heatmap_house(sharedCirc_DE_edgeR,get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),title_hp = "Heatmap of Different expressed shared_circRNA")
PCA_plot(DE_list,Circ_norm_edgeR,colData,withtext=FALSE)
PCA_plot(DE_list,Circ_norm_edgeR,colData,withtext=TRUE)
dev.off()
