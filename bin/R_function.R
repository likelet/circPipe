

pie_plot <- function(dt){
#colnames must be Freq and type
#####Freq is Freq or percent , type is type or label
library(ggplot2)
library(ggsci)
#dt = dt[order(dt$Freq, decreasing = TRUE),]   ##  order the Freq 
  
myLabel = paste(dt$type, " ( ", round(dt$Freq/sum(dt$Freq)*100,2), "% )   ", sep = "")  
p = ggplot(dt, aes(x = "", y = Freq, fill = type)) + 
  geom_bar(stat = "identity", width = 1) + scale_fill_jama(name="",label=myLabel)+ 
  coord_polar(theta = "y") + theme_void()+
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank())#, legend.position = "bottom")
     
return(p)
}

########heatmap with default color
heatmap_house <- function(dt,colData="",title_hp="Heatmap of Different expressed",cluster_rule="no"){
  
  ######suit two class 
  library("pheatmap")
  anno_color <-c("#3A4ABE", "#F33D02")
  comb_inFunc <- combn(levels(as.factor(colData$Type)),2)
  names(anno_color) <- c(comb_inFunc[2,1],comb_inFunc[1,1]) 
  ann_colors = list(Type= anno_color)
  dt=dt[,as.character(rownames(colData))]
  if (cluster_rule=="row") {
    p <- pheatmap(log2(dt+1),cluster_rows=TRUE, cluster_cols=FALSE,scale = "row",show_rownames=F,show_colnames = F,
                  color = colorRampPalette(c(rep('#1C2B6F',1),'white', rep('#E31E26',1)))(50),
                  annotation_col=colData,annotation_colors = ann_colors,main = title_hp)
  } else if (cluster_rule=="col"){
    p <- pheatmap(log2(dt+1),cluster_rows=FALSE,cluster_cols=TRUE, scale = "row",show_rownames=F,show_colnames = F,
                  color = colorRampPalette(c(rep('#1C2B6F',1),'white', rep('#E31E26',1)))(50),
                  annotation_col=colData,annotation_colors = ann_colors,main = title_hp)
  } else if (cluster_rule=="all"){
    p <- pheatmap(log2(dt+1),cluster_rows=TRUE,cluster_cols=TRUE, scale = "row",show_rownames=F,show_colnames = F,
                  color = colorRampPalette(c(rep('#1C2B6F',1),'white', rep('#E31E26',1)))(50),
                  annotation_col=colData,annotation_colors = ann_colors,main = title_hp)
  } else if (cluster_rule=="no"){
    p <- pheatmap(log2(dt+1),cluster_rows=FALSE,cluster_cols=FALSE, scale = "row",show_rownames=F,show_colnames = F,
                  color = colorRampPalette(c(rep('#1C2B6F',1),'white', rep('#E31E26',1)))(50),
                  annotation_col=colData,annotation_colors = ann_colors,main = title_hp)
    
  }
  
  return(p)
}

########boxplot 
box_plot <- function(df,x,y,fill){
library(ggplot2)
library(ggpubr)
  p <- ggbarplot(df, x=x, y=y,
            fill = fill ,#"#3B4992"
            color = "white",
            sort.val = "desc",#rule of sort 
            sort.by.groups=FALSE,
            x.text.angle=60)+scale_y_continuous(expand = c(0,0))

  return(p)
}


####################
# my_comparisons <- list(c("m6A", "non-m6A"))
# p <- ggboxplot(circ_input_RPM_m6A, x="m6Astatus", y="Log2mean", fill = "m6Astatus", 
#                palette = c("#374E55FF", "#DF8F44FF"))+
#   labs(title="The expression level of m6A circRNA and non-m6A circRNA", x="", y="log2(SRPBM+1)")+
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif")+#不同组间的比较
#   stat_compare_means(label.y = 18)+theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.9,0.85))


###################
volcano_plot <- function(res_df,contrast_factor,ylab_variable="FDR",pval = 0.05,fc=1.5){
  ######res_df is DE result data.frame with all genes (with non-sig gene)
  #for example , contrast_factor=c("normal","tumor")
  library(scales)
  i=1
  lfc=signif(log2(fc),2)
  if (ylab_variable=="pvalue"){
    tab = data.frame(logFC = res_df$logFC, negLogPval = -log10(res_df$pvalue))
  }else {
    tab = data.frame(logFC = res_df$logFC, negLogPval = -log10(res_df$FDR))
  }
  
  nosigGene = (abs(tab$logFC) < lfc | tab$negLogPval < -log10(pval))
  signGenes_up = (tab$logFC > lfc & tab$negLogPval > -log10(pval))
  signGenes_down = (tab$logFC < -lfc & tab$negLogPval > -log10(pval))
  gap = max(tab$logFC)/50
  up_count = length(which(signGenes_up))
  down_count = length(which(signGenes_down))
  #plot
  par(mar = c(5, 6, 5, 5))
  plot(tab, pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), cex.lab = 1.5, col = alpha("black", 0))
  points(tab[nosigGene, ], pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), col = "black", bg = "grey")
  if (length(unique(signGenes_up)) > 1){
    points(tab[signGenes_up, ], pch = 21, col = "black", bg = "red")
  }
  if (length(unique(signGenes_down)) > 1){
    points(tab[signGenes_down, ], pch = 21, col = "black", bg = "cornflowerblue")
  }
  abline(h = -log10(pval), col = "green3", lty = 2)
  abline(v = c(-lfc, lfc), col = "orange", lty = 2)
  mtext(c(paste("-", lfc), paste("+", lfc)), side = 1, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
  mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
  mtext(c(as.character(contrast_factor[2]), as.character(contrast_factor[1])), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line = 2)
  mtext(c(paste(up_count,"genes",sep = " "), paste(down_count,"genes",sep = " ")), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=0.5)
  legend("top",legend = c("Up regulate","Down regulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
  mtext(c(as.character(contrast_factor[2]), as.character(contrast_factor[1])), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line = 2)
  mtext(c(paste(up_count,"genes",sep = " "), paste(down_count,"genes",sep = " ")), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=0.5)
  legend("top",legend = c("Up regulate","Down regulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
  
}


################################
plot_for_DE <- function(DE_res0,DE_peaks,colData,norm.mat,rPackage,outpath){
  ##peak is a charactor of peakID
  ###colname of result is not in common
  colnames(DE_res0) <- gsub("log2FoldChange","logFC",colnames(DE_res0))
  colnames(DE_res0) <- gsub("P.Value","pvalue",colnames(DE_res0))
  colnames(DE_res0) <- gsub("PValue","pvalue",colnames(DE_res0))
  colnames(DE_res0) <- gsub("adj.P.Val","FDR",colnames(DE_res0))
  colnames(DE_res0) <- gsub("padj","FDR",colnames(DE_res0))
  colData$Type <- as.factor(colData$Type)
  comb <- combn(levels(as.factor(colData$Type)),2)
  #order by foldchange
  DE_res <- DE_res0[as.character(DE_peaks),]
  DE_norm_mat=norm.mat[rownames(DE_res[order(DE_res$logFC, decreasing = TRUE),]),]
  plot_file <- paste0(outpath,"/",comb[2,1], "_vs_", comb[1,1], "_DE_",rPackage,".pdf")
  cat("plot ",rPackage ,"result in dir ",plot_file)
  pdf(file = plot_file)
  ########Vocano plot
  volcano_plot(DE_res0,c(comb[1,1],comb[2,1]),ylab_variable = "pvalue",pval = 1,fc = 2)
  ########heatmap
  heatmap_house(DE_norm_mat,colData,title_hp = "Heatmap of Different methylation all DE RNA")
  heatmap_house(DE_norm_mat,colData,title_hp = "Heatmap of Different methylation all DE RNA",cluster_rule = "row")
  ########PCA
  PCA_plot(as.character(DE_peaks),norm.mat,colData)
  dev.off()
}


#########
plot_for_DM <- function(DM_res0,DM_peaks,colData,norm.mat,rPackage,outpath){
  ##peak is a charactor of peakID
  ###colname of result is not in common
  colnames(DM_res0) <- gsub("log2FoldChange","logFC",colnames(DM_res0))
  colnames(DM_res0) <- gsub("P.Value","pvalue",colnames(DM_res0))
  colnames(DM_res0) <- gsub("PValue","pvalue",colnames(DM_res0))
  colnames(DM_res0) <- gsub("adj.P.Val","FDR",colnames(DM_res0))
  colnames(DM_res0) <- gsub("padj","FDR",colnames(DM_res0))
  colData$Type <- as.factor(colData$Type)
  comb <- combn(levels(as.factor(colData$Type)),2)
  #order by foldchange
  DM_res <- DM_res0[as.character(DM_peaks),]
  DM_norm_mat=norm.mat[rownames(DM_res[order(DM_res$logFC, decreasing = TRUE),]),]
  plot_file <- paste0(outpath,"/",comb[2,1], "_vs_", comb[1,1], "_DM_",rPackage,".pdf")
  cat("plot ",rPackage ,"result in dir ",plot_file)
  pdf(file = plot_file)
  ########Vocano plot
  volcano_plot(DM_res0,c(comb[1,1],comb[2,1]),ylab_variable = "pvalue")
  #volcano_plot(DM_res0,c(comb[1,1],comb[2,1]),ylab_variable = "FDR")
  ########heatmap
  heatmap_house(DM_norm_mat,colData,title_hp = "Heatmap of Different expression all DM RNA")
  heatmap_house(DM_norm_mat,colData,title_hp = "Heatmap of Different expression all DM RNA",cluster_rule = "row")
  ########PCA
  PCA_plot(as.character(DM_peaks),norm.mat,colData)
  dev.off()
}


################################
#clussterprofile

GoKegg=function(gene_list,outdir){
  library(clusterProfiler)
  library(DOSE)
  library(ReactomePA)
  library(pathview)
  outpath=paste0(outdir,'/GO_KEGG')
  if (!dir.exists(outpath)){
    dir.create(outpath)
  }
  #options(bitmapType = "cairo")
  gene_id_list=bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  gene_id_list<-as.character(gene_id_list[,2])
  #KEGG
  kk <- enrichKEGG(gene = gene_id_list,
                   organism ="human",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   minGSSize = 1,
                   use_internal_data =FALSE)
  write.table(as.data.frame(kk@result), file=paste0(outpath,"/kegg.txt"),sep = "\t",quote = F)
  #pdf(paste("GO_KEGG/"kegg.pdf",sep = ""),bg = "transparent")
  pdf(paste0(outpath,"/kegg.pdf"))
  print(dotplot(kk,showCategory = 10))
  dev.off()
  cat ("cluster: kegg is ok~\n")
  
  pdf(paste0(outpath,"/GO_enrichment_analysis.pdf"))
  for (myont in c("CC","BP","MF")){
    ego <- enrichGO(gene=gene_id_list,
                    OrgDb=org.Hs.eg.db,
                    ont = myont,
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.1,
                    readable = TRUE)
    write.table(as.data.frame(ego@result), file=paste0(outpath,"/clust_",myont,".txt",sep = ""),sep = "\t",quote = F)
    print(dotplot(ego,showCategory = 10)+labs(title=paste("enrichment analysis of GO",myont)))
  }
  
  dev.off()
  cat("cluster: GO is ok~\n")
}

####################correlation analysis by sample
#cor matrix : group1 & group2
cor_for_mat <- function(cor_df.mat,tag1,tag2,cormethod="spearman"){
  print("default cor test spearman")
  group1 <- grep(tag1,colnames(cor_df.mat))
  group2 <- grep(tag2,colnames(cor_df.mat))
  #cor per row
  cor_tmp <- function(x,group1,group2,cormethod="spearman"){
    cor_res_obj <- cor.test(as.numeric(x[group1]),as.numeric(x[group2]),method = cormethod)
    cor_res<-c(x[-c(group1,group2)],cor_value=as.numeric(cor_res_obj$estimate), p_value=cor_res_obj$p.value)
    return(cor_res)
  }
  #cor for matrix
  cor_df <- as.data.frame(t(apply(cor_df.mat,1,cor_tmp,group1=group1,group2=group2,cormethod=cormethod)))
  cor_df$cor_value <-  as.numeric(as.character(cor_df$cor_value))
  cor_df$p_value <- as.numeric(as.character(cor_df$p_value))
  return(cor_df)
}



####################correlation plot

Cor_plot <- function(CoxExpPlotData,gene1,gene2,cormethod="spearman"){
  library("ggplot2")
	# cor_df = data.frame(CoxExpPlotData, 
	# 	gene1 = CoxExpPlotData[,gene1],
	# 	gene2 = CoxExpPlotData[,gene2])
	pcutoff=0.05
	CoxTest = cor.test(CoxExpPlotData[, gene1], CoxExpPlotData[, gene2], method = cormethod)

  if(cormethod == "pearson"){
    plottitle <- paste0(cormethod,", R = ", signif(CoxTest$estimate[[1]], 3), "\nP value = ", signif(CoxTest$p.value, 4), sep = "")
  }else if(cormethod == "spearman"){
    plottitle <- paste0(cormethod,", rho = ", signif(CoxTest$estimate[[1]], 3), "\nP value = ", signif(CoxTest$p.value, 4), sep = "")
  }else if (cormethod == "kendall"){
    plottitle <- paste0(cormethod,", tau = ", signif(CoxTest$estimate[[1]], 3), "\nP value = ", signif(CoxTest$p.value, 4), sep = "")
  }
  
  print(paste("pvalue =", round(CoxTest$p.value, 4)))
  if(CoxTest$p.value <= 1){
    CoxExpPoint = ggplot(data = CoxExpPlotData, aes_string(x= gene1, y = gene2))+theme_classic()+
      geom_point(size = 2, color = "gray36")+
      annotate("text",x=-Inf,y=Inf,label=plottitle,hjust=-.2,vjust=2)+
      stat_smooth(method = lm, se = FALSE, colour = "#191970")+
      ylab(paste(gene2, "Exp.", sep = " "))+xlab(paste(gene1, "Exp.",  sep = " "))+
      scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
      theme(plot.title = element_text(size = 15, angle = 0, face = "plain", colour = "black", hjust = 0.5, vjust = -0.5),
        axis.title.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.title.y = element_text(size = 15, angle = 90, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.text.y = element_text(size = 15, angle = 0, face = "plain", colour = "black"))}
  return(CoxExpPoint)
}  

##############ECDF plot

ECDF_plot <- function(df,value_var,group_var,plot_title,test_result=""){
  library(ggplot2)
  
  if (test_result$p.value <=0){
    test_anno <- paste(test_result$method,"\n",names(test_result$statistic)," = ",signif(test_result$statistic, 3),
                        "\nP Value < 2.2e-16",sep = "")
  } else {
    test_anno <- paste(test_result$method,"\n",names(test_result$statistic)," = ",signif(test_result$statistic, 3),
                        "\nP Value= ",signif(test_result$p.value,3),sep = "")
  }
  
  p <-  ggplot(df,aes_string(x=value_var,group=group_var,color=group_var))+theme_test()+
    stat_ecdf( size = 1)+theme(legend.position=c(0.85,.15))+
    annotate("text",x=-Inf,y=Inf,vjust=1.5,hjust=-.12,label=test_anno)+ 
    scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))+
    labs(title= plot_title, y="Cumulative fraction")+
    theme(plot.title = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
          axis.title.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
          axis.title.y = element_text(size = 15, angle = 90, face = "plain", colour = "black"),
          axis.text.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
          axis.text.y = element_text(size = 15, angle = 0, face = "plain", colour = "black"))
  return(p)
}


###################edgeR for DE analysis with mode (paired or unpaired)

get_pheno <- function(x,label1,label2,group1,group2){
  #####x is colnames of expr_mat
  pheno_data = data.frame(Type = c(rep(group1,length(grep(label1,x))), rep(group2,length(grep(label2,x)))))
  rownames(pheno_data) = c(x[grep(label1,x)],x[grep(label2,x)])
  rownames(pheno_data)=sub("X","S",as.character(rownames(pheno_data)))
  return(pheno_data)
}

edgeR_test <- function(expre_mat,group_mat="",design_mode = "unpaired",test_method="QLFit",adjust_method = "BH" , pvalue = 0.05, lfc = 0.58){
  library(edgeR)
  #group_mat <- get_pheno(colnames(expre_mat))
  expre_mat <- expre_mat[,as.character(rownames(group_mat))]
  treatment_inFunc <- as.factor(group_mat$Type)
  patient <- as.factor(gsub("[T|N]","",rownames(group_mat)))
  deg_lst <- DGEList(counts = expre_mat,genes = as.character(rownames(expre_mat)))
  #fliter &TMM 
  keep <- rowSums(cpm(deg_lst)>0) >= 2 #a CPM>0 in at least 2 samples
  deg_lst <- deg_lst[keep,]
  deg_lst <- calcNormFactors(deg_lst)
  #design matrix  ###########treatment VS control factor should be place on last column
  if (design_mode=="paired"){
    design_mat <- model.matrix(~patient+treatment_inFunc)
  } else {
    design_mat <- model.matrix(~treatment_inFunc)
  }
  rownames(design_mat)<-colnames(deg_lst)
  
  #disper and test  
  deg_lst<-estimateDisp(deg_lst,design_mat)
  fit_in_func <- glmQLFit(deg_lst,design_mat)
  if ( test_method=="QLFit"){
    mode_res <- glmQLFTest(fit_in_func)
    topTags(mode_res)
  } else if (test_method=="LRT") {
    mode_res <- glmLRT(fit_in_func)
    topTags(mode_res)
  } else {
    cat("please provide test method : QLFit or LRT")
  }
  de_inFunc<-decideTestsDGE(mode_res,adjust.method = "none" , p.value = pvalue, lfc = lfc)
  print("non-adjust")
  print(summary(de_inFunc))
  de_inFunc<-decideTestsDGE(mode_res,adjust.method = adjust_method , p.value = pvalue, lfc = lfc)
  print(adjust_method)
  print(summary(de_inFunc))
  res_inFunction <- data.frame(mode_res$table, FDR = p.adjust(mode_res$table$PValue, method=adjust_method))
  return(res_inFunction)
}

##########DE limma

limmaDE <- function(expres_dfsiondata,colData,k=2,cutoff=0.05) {
  #####colname of group factor in colData should be Type
  library(limma)
  class=colData$Type
  design<-model.matrix(~0+factor(class))
  colnames(design)<-paste("K",1:k,sep="")
  K=colnames(design)
  ###voom normalization also for counts matrix
  #v <- voom(expres_dfsiondata, design, normalize.method = "quantile")
  fit<-lmFit(expres_dfsiondata,design) #fit linear model
  x <- c()
  for(i in 1:(k-1)){
    for(j in (i+1):k){
      x <- c(x, paste(K[j],K[i],sep="-"))
    }
  }
  contrast.matrix<-makeContrasts(contrasts=x,levels=design)#construct contrast matrix
  fit2<-contrasts.fit(fit,contrast.matrix) #get microarray linear model
  fit2<-eBayes(fit2) #empirical Bayes moderation,compute moderated t-statistics,modeated F-statistic
  output <- topTable(fit2,sort.by= "p",number=Inf)
  #genesignature<-topTable(fit2,number=20000,adjust.method='BH',p.value=cutoff)
  r <- list()
  #if(k>2){
  for(i in 1:(k*(k-1)/2)){
    geneset_i<-topTable(fit2,number=Inf,coef=i,adjust.method='BH',sort = "p",p.value=1)
    cat("sig peaks Pvalue(0.05)",sum(geneset_i$P.Value <= 0.05),"\n")
    cat("sig peaks adjust P(0.05)",sum(geneset_i$adj.P.Val <= 0.05))
    r[[length(r)+1]] <- geneset_i
  }
  # }
  return(r)
}


##########a function for wilcoxon sum rank test for expre_matrix
row_wilcox <- function(group_mat,x,test_mode=""){
  group_mat$Type <- as.factor(group_mat$Type)
  typelevel <- levels(as.factor(group_mat$Type))
  comb_2 <- combn(typelevel,2)
  group1 <- as.character(rownames(subset(group_mat,Type==comb_2[1])))
  group2 <- as.character(rownames(subset(group_mat,Type==comb_2[2])))
  if (test_mode=="paired"){
    res_wix0 <- wilcox.test(x[which(rownames(group_mat)%in%group1)],x[which(rownames(group_mat)%in%group2)],paired = T)
  } else {
    res_wix0 <- wilcox.test(x[which(rownames(group_mat)%in%group1)],x[which(rownames(group_mat)%in%group2)])
  }
  
  res_wix <- c(Pvalue=res_wix0$p.value,statistic=res_wix0$statistic)
  return(res_wix)
}

mat_wilcox <- function(expre_mat,group_mat,test_mode=""){
  cat("peak number ",dim(expre_mat)[1],"\n")
  expre_mat <- expre_mat[,as.character(rownames(group_mat))]
  res_wix_lst <- apply(expre_mat,1,function(x) row_wilcox(group_mat,x,test_mode))
  res_wix_lst = as.data.frame(t(res_wix_lst))
  res_wix_lst$FDR = p.adjust(res_wix_lst$Pvalue,method = "BH")
  res_wix_lst$BY = p.adjust(res_wix_lst$Pvalue,method = "bonferroni")
  cat("DM peaks Pvalue(0.05)",sum(res_wix_lst$Pvalue <=0.05),"\n")
  cat("DM peaks FDR(0.05)",sum(res_wix_lst$FDR <=0.05),"\n")
  return(res_wix_lst)
}


#############PCA
 PCA_plot<- function(DE_list,df,colData,withtext){
  colData$Type <- as.factor(colData$Type)
  comb <- combn(levels(as.factor(colData$Type)),2)
  colData <- subset(colData, Type == comb[1,1] | Type == comb[2,1])
  library(ggplot2)
  df <- df[DE_list, rownames(colData)] ########normalized matrix
  df = log2(df+1)
  pcaData <- as.data.frame(prcomp(df[DE_list,])$rotation)
  if(withtext==FALSE){
    pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=colData$Type)) +
    geom_point(size=3) +
    xlab("PC1") +
    ylab("PC2") +
    scale_colour_hue("Type") +
    #  coord_fixed() +
    theme_bw()
  print(pca_plot)
  }else{
    pca_plot_text <- ggplot(pcaData, aes(PC1, PC2, color=colData$Type)) +
    geom_text(aes(label = row.names(pcaData))) +
    xlab("PC1") +
    ylab("PC2") +
    scale_colour_hue("Type") +
    #  coord_fixed() +
    theme_bw()
  print(pca_plot_text)
  }
}