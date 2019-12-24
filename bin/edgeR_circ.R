library("readr")
library(dplyr)
library(ggplot2)

#source("/home/wqj/code/circPipe-master/R_function.R")

options(stringsAsFactors=F)

setwd("./")
library(Biobase)
library(edgeR)
########cut off init
p_value <- 0.05 # p value of FDR
lfc <- 0.58  #logFC
args <-commandArgs(T)

source(args[1])

########load data
# args[2]<-"/home/wqj/code/segemehl_merge.matrix"
# args[2]<-"/home/wqj/code/circPipe-master/onlychrx.matrix"
# countData1 <- as.data.frame(read.table(args[2],sep = "\t"))
countData1 <- as.data.frame(read_delim(args[2],delim = "\t"))
countData1 <- countData1[!duplicated(countData1$id), ]
dim(countData1)
# colnames(countData1) <- countData1[1,]
# countData1<-countData1[-1,]
rownames(countData1) <- countData1$id
#countData1=countData1[,-(2:6)]
countData=countData1[,-1]
# countData<-sapply(countData, as.numeric)
# countData<-data.frame(countData)
# rownames(countData)<-countData1$id
# countData[countData==1]=0
###########pheno data loading
# args[3]<-"/data1/wqj/database/data/12_21/design.file"
colData1 <- as.data.frame(read_delim(args[3],delim = "\t"))
colData<-data.frame(Type = colData1$Type)
rownames(colData)<-colData1[,1]
# colData <- get_pheno(x = colnames(countData),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal")
group1 = "Tumor"
group2 = "Normal"
label1 = "T"
label2 = "N"
# colnames(countData)=sub("X","S",as.character(colnames(countData)))
countData <- countData[,as.character(rownames(colData))]
colData$Type <- as.factor(colData$Type)
# args[4]<-"/home/wqj/code/compare.file"
temp<-read.table(args[4])
type_level <- rev(strsplit((temp[1,1]),split = "_vs_")[[1]])
#type_level <- levels(as.factor(colData$Type))
comb <- combn(type_level,2)

Circ_norm_edgeR=cpm(countData)
#########filter low-abandance circRNA ; the step has been done in node1:filter_circ
 countData <- countData[which(rowSums(countData > 0) >= 2),] #a Count>0 in at least 2 samples
 dim(countData)

######file for KEGG
 circ_feature=read.delim(args[5],header = FALSE,sep="\t",fill = T,col.names = c("chr","start","end","circ_type","exon_number","strand","strand2","ensemble_id","symbol","transcript","gene_feature","type1","type2","m6Astatus"))
 
##############edgeR

#sharedCirc_edgeR_tmp <- edgeR_test(expre_mat = shared_circ,group_mat = colData,test_method = "LRT" )
if( sum(colData$Type==group1)< 3 || sum(colData$Type==group2)< 3 ){
  # mean.mat = data.frame(group1_mean = rowMeans(Circ_norm_edgeR[,grep(label1,colnames(Circ_norm_edgeR)), drop = FALSE]),
  #                       group2_mean =rowMeans(Circ_norm_edgeR[,grep(label2,colnames(Circ_norm_edgeR)),drop = FALSE]),
  #                       all_mean = rowMeans(Circ_norm_edgeR))
  # rownames(mean.mat) = rownames(Circ_norm_edgeR)
  # mean.mat = subset(mean.mat,group1_mean > 0)
  # 
  # sharedCirc_edgeR = data.frame(logFC=log2(mean.mat$group2_mean/mean.mat$group1_mean+1),
  #                       logCPM = log2(mean.mat$all_mean),
  #                   `F`= rep("_",nrow(mean.mat)),PValue = rep("_",nrow(mean.mat)),
  #                   FDR =rep(0,nrow(mean.mat)))
  # rownames(sharedCirc_edgeR) = rownames(mean.mat)
  notice="No enough samples!"
  write.table(notice,file = "Message.txt",sep = ' ')
}else {
  sharedCirc_edgeR <-edgeR_test(expre_mat = countData,group_mat = colData)
  DE_sharedCirc_res = subset(sharedCirc_edgeR,FDR <= p_value & abs(logFC) >= lfc)
  if(nrow(DE_sharedCirc_res)==0){
    notice="No differential genes!"
    write.table(notice,file = "Message.txt",sep = ' ')   
  }else{
    select <- DE_sharedCirc_res[order(DE_sharedCirc_res$logFC, decreasing = TRUE),]
    
    rownames(Circ_norm_edgeR)=countData1$id
    sharedCirc_DE_edgeR=Circ_norm_edgeR[rownames(select),]
    #all DE
    DE_list <- c(rownames(DE_sharedCirc_res))
    
    ################plot
    #######Vocano plot
    png("volcaon.png", width=1200, height=1200, res = 300)
    volcano_plot(sharedCirc_edgeR,c(group1,group2))
    dev.off()
    
    ########heatmap
    # png("heatmap1.png", width=1800, height=1800, res = 300)
    # heatmap_house(Circ_norm_edgeR[DE_list,],colData,title_hp = "Heatmap of Different expressed all DE circRNA")
    # dev.off()
    png("heatmap2.png", width=1200, height=1200, res = 300)
    heatmap_house(Circ_norm_edgeR[DE_list,],colData,
                  title_hp = "Heatmap of Different expressed all DE circRNA",cluster_rule = "row")
    dev.off()
    png("heatmap3.png", width=1200, height=1200, res = 300)
    heatmap_house(sharedCirc_DE_edgeR,colData,title_hp = "Heatmap of Different expressed shared_circRNA")
    dev.off()
    ########PCA
    # png("pca1.png", width=1800, height=1800, res = 300)
    # PCA_plot(DE_list,Circ_norm_edgeR,colData,withtext=FALSE)
    # dev.off()
    png("pca2.png", width=1200, height=1200, res = 300)
    PCA_plot(DE_list,Circ_norm_edgeR,colData,withtext=TRUE)
    dev.off()
    
    pdf("all_plot.pdf")
    volcano_plot(sharedCirc_edgeR,c(group1,group2))
    # heatmap_house(Circ_norm_edgeR[DE_list,],get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),title_hp = "Heatmap of Different expressed all DE circRNA")
    heatmap_house(Circ_norm_edgeR[DE_list,],get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),
                  title_hp = "Heatmap of Different expressed all DE circRNA",cluster_rule = "row")
    heatmap_house(sharedCirc_DE_edgeR,get_pheno(colnames(sharedCirc_DE_edgeR),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal"),title_hp = "Heatmap of Different expressed shared_circRNA")
    # PCA_plot(DE_list,Circ_norm_edgeR,colData,withtext=FALSE)
    PCA_plot(DE_list,Circ_norm_edgeR,colData,withtext=TRUE)
    dev.off()
    
    
    ###output table
    DE_sharedCirc_res$id=rownames(DE_sharedCirc_res)
    
    final_DE_res <- rbind(DE_sharedCirc_res)
    
    final_DE_res$logFC <- as.numeric(final_DE_res$logFC)
    final_DE_res = mutate(final_DE_res,
                          expres_regu = case_when(
                            logFC >= lfc ~ "up_regu",
                            logFC <= -lfc ~ "down_regu"
                          ))
    
    write.table(final_DE_res, file = "DE_result.xls",sep = "\t", quote = FALSE, row.names = F)
    
    ####################DE enrichment analysis
    #circ_feature= read_delim("/home/wqj/code/circPipe-master/uniq_circ_m6Astatus_anno.txt","\t")
    circ_name = paste(circ_feature$chr,circ_feature$start,circ_feature$end,circ_feature$strand,sep = "_")
    circ_feature$position = circ_name
    
    DE_circ_gene = filter(circ_feature,position%in%DE_list)
    gene_list=unique(DE_circ_gene$symbol)
    length(gene_list)
    
    
    GoKegg(gene_list,"./")
  }
  
}



