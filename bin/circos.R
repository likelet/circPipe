library(readr)
library(ggplot2)
library(circlize)
library(dplyr)
library(RColorBrewer)
options(stringsAsFactors=F)
args <-commandArgs(T)
setwd("./")
#args[1]<-"/home/wqj/code/circPipe-master/final.matrix"
finalmatrix <- as.data.frame(read_delim(args[1],delim = "\t"))
id<-as.character(finalmatrix$id)
chr<-sapply(strsplit(id,"_"),"[",1)
start<-sapply(strsplit(id,"_"),"[",2)
end<-sapply(strsplit(id,"_"),"[",3)
strand<-sapply(strsplit(id,"_"),"[",4)
a <- data.frame(chr = chr, start = as.numeric(start), end = as.numeric(end), value = rowSums(finalmatrix[,2:ncol(finalmatrix), drop = FALSE]))
#a <- data.frame(chr = len_m6Astatus2$chr, start = len_m6Astatus2$start, end = len_m6Astatus2$end, value = len_m6Astatus2$length)
a=filter(a, chr != "chrM")
b<-a$chr
b<-factor(b,ordered = TRUE,levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
a$chr<-b
a<-a[order(a[,1]),]
a$chr<-as.character(a$chr)
b=unique(a$chr)

png("circos.png",width=1800, height=1800, res = 300)
par(mar=c(1,1,1,1))
circos.initializeWithIdeogram(species = 'hg19',chromosome.index = b)
circos.genomicTrackPlotRegion(a,panel.fun = function(region, value, ...){
  circos.genomicPoints(region, value, cex = 0.3, pch = 16, col='red',...)})
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 12)
circos.trackHist(a$chr, a$start, bg.col = bg.col, col = "blue",bin.size = 3000000)
circos.clear()
dev.off()
