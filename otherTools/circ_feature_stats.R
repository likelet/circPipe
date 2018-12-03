
####distribution
library(readr)
library(ggplot2)
library(circlize)
library(dplyr)
args <-commandArgs(T)
source(args[1])
#circ_feature= read_delim("/home/wqj/code/circPipe-master/uniq_circ_m6Astatus_anno.txt","\t")
circ_feature= read_delim(args[2],"\t",col_names = FALSE)
colnames(circ_feature)<-c("chr","start","end","circ_type","exon_number","strand","strand2","ensemble_id","symbol","transcript","gene_feature","type1","type2","m6Astatus")
circ_name = paste(circ_feature$chr,circ_feature$start,circ_feature$end,circ_feature$strand,sep = "_")
circ_feature$position = circ_name

#circ_feature_m6A= subset(circ_feature,m6Astatus=="m6A")
#circ_feature_nonm6A= subset(circ_feature,m6Astatus=="non-m6A")
########feature: gene type 
plot_type_pie <- function(feature_mat) {
  df=as.data.frame(table(feature_mat$type1))
  names(df)[1] <- "type"
  df$Freq=as.numeric(as.character(df$Freq))
  p=pie_plot(df)
  return(p)
}

#pdf(file = "plot/11_7/genomic distribution of circRNAs.pdf", width=10.31,height=6.23)
#m6A_typepie=plot_type_pie(circ_feature_m6A)
#m6A_typepie+ggtitle("Genomic distribution of m6A circRNAs")
#nonm6A_typepie=plot_type_pie(circ_feature_nonm6A)
#nonm6A_typepie+ggtitle("Genomic distribution of nonm6A circRNAs")
all_circ_typepie=plot_type_pie(as.data.frame(circ_feature))

png(args[3],width=1200, height=1200, res = 300)
all_circ_typepie+ggtitle("Genomic distribution of all circRNAs")
dev.off()

# tmp <- xtabs(~m6Astatus+type1,circ_feature)
# tmp <- chisq.test(tmp)
# tmp$p.value
# ##################exon number################
# # circ_feature_m6A= m6A_circ
# # circ_feature_nonm6A= non_m6A_circ
# 
# within7_freq <- function(circ_feature_df){
#   circ_feature_stats <- as.data.frame(table(filter(circ_feature_df,circ_type != "ciRNA")$exon_number))
#   names(circ_feature_stats)[1] <- "exon_number"
#   #circ_feature_stats$exon_number <- as.character(circ_feature_stats$exon_number )
#   circ_feature_stats$Freq <- as.numeric(as.character(circ_feature_stats$Freq))
#   circ_feature_stats$exon_number <- as.character(circ_feature_stats$exon_number )
#   circ_stats_over7 <- c(paste0("over",circ_feature_stats[7,1]),sum(circ_feature_stats$Freq[8:dim(circ_feature_stats)[1]]))
#   circ_feature_stats <- rbind(circ_feature_stats[1:7,],circ_stats_over7)
#   circ_feature_stats$Freq <- as.numeric(as.character(circ_feature_stats$Freq))
#   return(circ_feature_stats)
# }
# 
# circ_feature_m6A_stats = within7_freq(circ_feature_m6A)
# circ_feature_nonm6A_stats = within7_freq(circ_feature_nonm6A)
# circ_feature_m6A_stats$type <- rep("m6A",dim(circ_feature_m6A_stats)[1])
# circ_feature_nonm6A_stats$type <- rep("nonm6A",dim(circ_feature_nonm6A_stats)[1])
# ###########
# all_circ_exon <- rbind(circ_feature_m6A_stats,circ_feature_nonm6A_stats)
# all_circ_exon$Freq <- as.numeric(as.character(all_circ_exon$Freq))
# #####percent
# all_circ_exon$Percentage <- c(circ_feature_m6A_stats$Freq/sum(circ_feature_m6A_stats$Freq),
#                               circ_feature_nonm6A_stats$Freq/sum(circ_feature_nonm6A_stats$Freq))
# 
# all_circ_exon2 = filter(all_circ_exon, exon_number != "over7")
# #pdf(file = "plot/11_7/Percentage of exons circRNAs spanning.pdf")
# ggplot(all_circ_exon2, aes(x = exon_number, y = Freq,fill=type)) + scale_fill_jama()+theme_classic()+
#    theme(legend.position=c(0.9,0.85))+scale_y_continuous(expand = c(0,0))+#scale_x_continuous(expand = c(0,0))+
#   geom_bar(stat = "identity",position = "dodge")+ggtitle("the number of exons circRNAs spanning ")
# 
# ggplot(all_circ_exon2, aes(x = exon_number, y = Percentage,fill=type)) + scale_fill_jama()+theme_classic()+
#   theme(legend.position=c(0.9,0.85))+scale_y_continuous(expand = c(0,0))+
#   geom_bar(stat = "identity",position = "dodge")+ggtitle("the number of exons circRNAs spanning ")
# #dev.off()


#### the length of circRNA

circ_feature$length = circ_feature$end-circ_feature$start
sum(circ_feature$length > 50000)
len_m6Astatus = filter(circ_feature, length>100 & length<10000)

summary(len_m6Astatus$length)

png(args[4],width=1200, height=1200, res = 300)
boxplot(len_m6Astatus$length)
dev.off()

len_df= as.data.frame(len_m6Astatus$length)
#plot(ecdf(len_m6Astatus2$length),do.point=F,verticals=T)

png(args[5],width=2400, height=1200, res = 300)
ggplot()+stat_ecdf(data = len_m6Astatus,aes(x = length), 
                   geom = "step", size = 1) +scale_color_lancet()+theme_bw()+
  ggtitle("the distrubution of circRNA spanning distance")+ylab("Frequency")+xlab("Length")
dev.off()


#####new code
#circ_feature$length = circ_feature$end-circ_feature$start

all=sum(circ_feature$length > 10000)
len_m6Astatus2 = filter(circ_feature, circ_type != "ciRNA", length < 10000)

len_df2= data.frame(length = len_m6Astatus2$length)

x_labs.txt = c(100,seq(1000,10000, 1000))
circ_length.hist <- ggplot(data = len_df2, aes(x = length))+
  geom_histogram(fill ="blue" ,color = "white",breaks = x_labs.txt)+ #color"#00006D"
  ylab("Frequency")+xlab("")+
  scale_x_continuous(limits  = c(0,10000),breaks =x_labs.txt)+
  scale_y_continuous(limits  = c(0,all),breaks = seq(0,all,all/5),expand = c(0,0))+theme(
    axis.title = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black",hjust = 1,vjust = 0.9,angle = 60),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.key.size = unit(0.3,"cm"),
    legend.title= element_blank(),
    legend.text = element_text(size = 10, face = NULL, color = "black"), 
    legend.background = element_rect(fill = NA),
    panel.background = element_blank())
# ggtitle("the distrubution of circRNA spanning distance")+
# ylab("Frequency")

png(args[6],width=2400, height=1200, res = 300)
circ_length.hist
dev.off()

######draw circle
a <- data.frame(chr = len_m6Astatus2$chr, start = len_m6Astatus2$start, end = len_m6Astatus2$end, value = len_m6Astatus2$length)
a=filter(a, chr != "chrM")

png(args[7],width=2400, height=2400, res = 300)
circos.genomicInitialize(a,plotType="NULL")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex = 1,
              facing = "clockwise", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)
circos.genomicTrackPlotRegion(a,panel.fun = function(region, value, ...){
  circos.genomicLines(region, value, type = "s",col='green',...)})
col <- rep(c("#FF0000", "#00FF00"), 12)
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 12)
circos.trackHist(a$chr, a$end, bg.col = bg.col, col = col)
circos.genomicTrackPlotRegion(a,panel.fun = function(region, value, ...){
  circos.genomicPoints(region, value, cex = 0.1, pch = 16, col='red',...)})
circos.clear()
dev.off()



pdf(args[8])
all_circ_typepie+ggtitle("Genomic distribution of all circRNAs")
boxplot(len_m6Astatus$length)
ggplot()+stat_ecdf(data = len_m6Astatus,aes(x = length), 
                   geom = "step", size = 1) +scale_color_lancet()+theme_bw()+
  ggtitle("the distrubution of circRNA spanning distance")+ylab("Frequency")+xlab("Length")
circ_length.hist
circos.genomicInitialize(a,plotType="NULL")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex = 1,
              facing = "clockwise", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)
circos.genomicTrackPlotRegion(a,panel.fun = function(region, value, ...){
  circos.genomicLines(region, value, type = "s",col='green',...)})
col <- rep(c("#FF0000", "#00FF00"), 12)
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 12)
circos.trackHist(a$chr, a$end, bg.col = bg.col, col = col)
circos.genomicTrackPlotRegion(a,panel.fun = function(region, value, ...){
  circos.genomicPoints(region, value, cex = 0.1, pch = 16, col='red',...)})
circos.clear()
dev.off()



# ###第一种画法
# circos.genomicInitialize(a,plotType="NULL")
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   chr = CELL_META$sector.index
#   xlim = CELL_META$xlim
#   ylim = CELL_META$ylim
#   circos.text(mean(xlim), mean(ylim), chr, cex = 1,
#               facing = "clockwise", niceFacing = TRUE)
# }, track.height = 0.1, bg.border = NA)
# circos.genomicTrackPlotRegion(a,panel.fun = function(region, value, ...){
#   circos.genomicPoints(region, value, cex = 0.01, pch = 16, col="red",...)})
# 
# circos.genomicTrackPlotRegion(a,panel.fun = function(region, value, ...){
#   circos.genomicLines(region, value, type = "s",col='green',...)})
# 
# circos.genomicTrackPlotRegion(a,panel.fun = function(region, value, ...) {
#   cex = (value[[2]] - min(value[[2]]))/(max(value[[2]]) - min(value[[2]]))#
#   i = getI(...)
#   circos.genomicPoints(region, value, cex = 0.5, pch = 1, col = "red", ...)})
# circos.clear()
# 
# #####第二种画法
# circos.initializeWithIdeogram(a,plotType = NULL)
# circos.trackPlotRegion(ylim = c(0,1) , bg.border= NA, panel.fun = function(x,y){
#   circos.text(CELL_META$xcenter,CELL_META$ycenter, CELL_META$sector.index, cex = 0.8,font = 2, facing = "clockwise", niceFacing = TRUE)
# })
# # par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)
# # circos.par(track.height = 0.1)
# # circos.initialize(factors = a$chr, x = a$start) #初始化，factors来控制track数目，初始化里只有x， 没有y。这一步相当于ggplot()
# circos.trackPlotRegion(factors = a$chr, x = a$start, y = a$end,
#                        panel.fun = function(x, y) { 
#                          circos.axis()})
# col <- rep(c("#FF0000", "#00FF00"), 12) #自定义一下颜色# 这里先解释一下，一个track有好几个cell，具体数目由factors决定的，向本数据集中factors有八个，因此绘制一个track，其包含八个cell。含有前缀circos.track的函数会在所有的cel里添加基本元素，而只有前缀circos.的函数可以在特定的track、cell里添加基本元素。具体看下演示。
# circos.trackPoints(a$chr, a$end, col = col, pch = 16, cex = 0.5) #所有的cell里都绘制点图
# #circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1) #在track 1中的标记为a的cell里添加
# #circos.text(1, 0.5, "right", sector.index = "a")
# bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 12)
# circos.trackHist(a$chr, a$end, bg.col = bg.col, col = col)
# circos.trackPlotRegion(factors = a$chr, x = a$start, y = a$end, 
#                        panel.fun = function(x, y) {
#                          grey = c("#FFFFFF", "#CCCCCC", "#999999") 
#                          sector.index = get.cell.meta.data("sector.index") #这个是第三个track，因为我们刚刚创建，这里这一步不用也可。
#                          xlim = get.cell.meta.data("xlim")
#                          ylim = get.cell.meta.data("ylim") 
#                          circos.points(x[1:mean(x)], y[1:mean(x)], col = "red", pch = 16, cex = 0.6) 
#                          circos.points(x[mean(x):max(x)], y[mean(x):max(x)], col = "blue", cex = 0.6)})
# # # update第2个track中标记为d的sector
# # circos.updatePlotRegion(sector.index = "d", track.index = 2)
# # circos.points(x = -2:2, y = rep(0, 5))
# # xlim <- get.cell.meta.data("xlim")
# # ylim <- get.cell.meta.data("ylim")
# # circos.text(mean(xlim), mean(ylim), "updated")
# circos.trackPlotRegion(factors = a$factors, y = a$y)
# circos.trackLines(a$factors[1:100], a$x[1:100], a$y[1:100], type = "h")
# # circos.link("a", 0, "b", 0, h = 0.3) #point to point
# # circos.link("c", c(-0.5, 0.5), "d", c(-0.5, 0.5), col = "red", border = NA, h = 0.2) #intreval to interval
# # circos.link("e", 0, "g", c(-1, 1), col = "green", border = "black", lwd = 2, lty = 2) #point to interval
# 
# circos.clear()

