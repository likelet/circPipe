

#########load packages##############
library(data.table)
library(dplyr)
library(ggplot2)
options(stringsAsFactors=F)

########consent init######
p_value <- 0.05 # p value of FDR
lfc <- log2(1.5)  #logFC

#######arguments#########

args <-commandArgs(T)
function.script <- args[1] 
expres.mat.filePath <- args[2] 
design.filePath <- args[3] 
compare.filePath <- args[4] 
annotation.filePath <- args[5] 
outkegg<-args[6]
#source function script
source(function.script)

########load data########
design.mat = as.data.frame(fread(design.filePath))
compare.str = as.vector(read.table(compare.filePath))
expres.mat = as.data.frame(fread(expres.mat.filePath,sep = "\t"))
colnames(expres.mat)[1]="id"
annotation.df = as.data.frame(fread(annotation.filePath,sep = "\t"))%>%
  mutate(id = paste(chr,chromStart,chromEnd,strand,sep="_"))

expres.mat = select(expres.mat,id,design.mat$Sample_id)
countData = as.matrix(expres.mat[,-1])
rownames(countData) = expres.mat$id

type_level <- rev(strsplit((compare.str[1,1]),split = "_vs_")[[1]])#type_level <- levels(colData$Type)
colData = data.frame(Type= factor(design.mat$Type,levels = type_level))
comb <- combn(type_level,2)
rownames(colData) = design.mat$Sample_id

#########data clean#########
countData <- countData[,as.character(rownames(colData))]
countData <- countData[which(rowSums(countData > 0) >= 2),] #a Count>0 in at least 2 samples
dim(countData)

###########too few sample###########
if (length(design.mat$Sample_id) < 4 ){
  print("too few samples")
  quit()
}

##############edgeR############
length(unique(rownames(countData)))
#sharedCirc_edgeR_tmp <- edgeR_test(expre_mat = shared_circ,group_mat = colData,test_method = "LRT" )
Circ_edgeR.res <-edgeR_test(expre_mat = countData,group_mat = colData)
Circ_edgeR.res$id = rownames(Circ_edgeR.res)
Circ_edgeR.res.anno = mutate(Circ_edgeR.res,
                             expres_regu = case_when(logFC >= lfc ~ "up_regu",
                                                               logFC <= -lfc ~ "down_regu",
                                                               abs(logFC) < lfc ~ "not-sig"))%>%inner_join(annotation.df)

DE_Circ.res = filter(Circ_edgeR.res.anno,FDR <= p_value & abs(logFC) >= lfc)%>%arrange(desc(logFC))
DE_list = as.character(DE_Circ.res$id)
Circ_norm.mat=cpm(countData)
DE_Circ_norm.mat=Circ_norm.mat[DE_list,]


####output table##########
#revised by yy
write.table(Circ_edgeR.res.anno,file = paste0(outkegg,"/DE_analysis_result.xls"),sep = "\t",row.names = FALSE,quote = FALSE)

if (length(DE_list) >= 2) {

# ################plot#################
source(function.script)
#pdf(file = paste0("DE/",comb[2,1], "_vs_", comb[1,1], "_DE.edgeR.pdf"))
########Vocano plot####
volcano_plot(Circ_edgeR.res,type_level)
########heatmap####
heatmap_house(Circ_norm.mat[DE_list,],colData,title_hp = "Heatmap of Different expressed circRNA")
heatmap_house(Circ_norm.mat[DE_list,],colData,
              title_hp = "Heatmap of Different expressed circRNA",cluster_rule = "row")
heatmap_house(sharedCirc_DE_edgeR,colData,title_hp = "Heatmap of Different expressed circRNA")
########PCA#######
PCA_plot(DE_list,Circ_norm.mat[DE_list,],colData)

dev.off()

###########DE enrichment analysis###########

DE_circ_gene = filter(Circ_edgeR.res.anno,id%in%DE_list)
gene_list=unique(DE_circ_gene$symbol)
length(gene_list)

GoKegg(gene_list,outkegg,"DE_circ")


} else {

print("differnt expression circRNA < 2 ; cannot plot and enrichment analysis !")
}
