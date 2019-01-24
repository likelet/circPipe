#####distribution
library(readr)
library(ggplot2)
library(circlize)
library(dplyr)
library(RColorBrewer)
options(stringsAsFactors=F)
setwd("./")

args <-commandArgs(T)

source(args[1])

circ_feature=read.delim(args[2],header = FALSE,sep="\t",fill = T,col.names = c("chr","start","end","circ_type","exon_number","strand","strand2","ensemble_id","symbol","transcript","gene_feature","type1","type2","m6Astatus"))

#source("/home/wqj/code/circPipe-master/R_function.R")
#circ_feature= read_delim("matrix/uniq_circ_m6Astatus_anno.txt","\t")
#circ_feature= read_delim("matrix/11_7/all_filter_circ_m6Astatus_anno.txt","\t")
circ_name = paste(circ_feature$chr,circ_feature$start,circ_feature$end,circ_feature$strand,sep = "_")
circ_feature$position = circ_name


##########plot init############
##font size
size_axisblacktext <- 18 #24
size_axisgreytext <- 12 #18
size_legendtext <- 9  #22
size_axistitle <- 15

##colors for ggplot 
ggcolors <- c("#BC3C29","#0072B5","#E18727","#20854E","#80397B")#pal_npg("nrc")(5)#gg_color_hue(5)
#show_col(ggcolors)

plotcolor = c("#BC3C29", "grey", "#0072B5")

## basic config for all ggplot
gg_theme <- theme(
  axis.title = element_text(size = size_axistitle, color = "black"),
  axis.text.x = element_text(size = size_axisgreytext, color = "black"),
  axis.text.y = element_text(size = size_axisblacktext, color = "black"),
  legend.key.size = unit(0.3,"cm"),
  legend.title= element_blank(),
  legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
  legend.background = element_rect(fill = NA))


################feature: gene type########################
plot_type_pie <- function(feature_mat) {
  df=as.data.frame(table(feature_mat$type1))
  names(df)[1] <- "type"
  df$Freq=as.numeric(as.character(df$Freq))
  p=pie_plot(df)
  return(p)
}

all_circ_typepie=plot_type_pie(as.data.frame(circ_feature))

png("genomic_distribution.png",width=1800, height=1800, res = 300)
all_circ_typepie+ggtitle("Genomic distribution of all circRNAs")
dev.off()

# tmp <- xtabs(~m6Astatus+type1,circ_feature)
# tmp <- chisq.test(tmp)
# tmp$p.value

################### the length of circRNA#################

circ_feature$length = circ_feature$end-circ_feature$start

all=sum(circ_feature$length < 10000)
len_m6Astatus2 = filter(circ_feature, circ_type != "ciRNA", length < 10000)

len_df= data.frame(length = len_m6Astatus2$length)

x_labs.txt = c(100,seq(1000,10000,1000))
circ_length.hist <- ggplot(data = len_df, aes(x = length))+
  geom_histogram(fill =ggcolors[5] ,color = "white",breaks = x_labs.txt)+ #color"#00006D"
  ylab("Frequency")+xlab("")+
  scale_x_continuous(limits  = c(0,10000),breaks =x_labs.txt)+
  scale_y_continuous(limits  = c(0,all),breaks = seq(0,all/2,all/5),expand = c(0,0))+theme(
    axis.title = element_text(size = size_axistitle, color = "black"),
    axis.text.x = element_text(size = size_axisgreytext, color = "black",hjust = 1,vjust = 0.9,angle = 60),
    axis.text.y = element_text(size = size_axisblacktext, color = "black"),
    legend.key.size = unit(0.3,"cm"),
    legend.title= element_blank(),
    legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
    legend.background = element_rect(fill = NA))
  # ggtitle("the distrubution of circRNA spanning distance")+
  # ylab("Frequency")

png("frequency.png",width=2400, height=1200, res = 300)
circ_length.hist
dev.off()

################

test=read.delim(args[3],header = FALSE,sep="\t",fill = T)
distance <- test$V14
feature <- test$V13


get_density_df <- function(raw_density_df,region_type){
  # region_type="circRNA"
  raw_density_df <- filter(raw_density_df,type==region_type)
  
  peak_freq <- c(filter(raw_density_df,feature == "5UTR")$distance , 
                 filter(raw_density_df,feature == "3UTR")$distance + 200,
                 filter(raw_density_df,feature == "CDS")$distance + 100)
  peak_freq2 <- data.frame(as.data.frame(table(peak_freq)/length(peak_freq)),
                           type= rep(region_type,nrow(as.data.frame(table(peak_freq)/length(peak_freq)))))
  peak_freq3 <- data.frame(distance=peak_freq,type=rep(region_type,length(peak_freq)))
  res_lst <- list(stats = peak_freq2,row =peak_freq3)
  return(res_lst)
}

circRNA_distance_df <- data.frame(distance=as.numeric(distance),feature=feature,type="circRNA")
circRNA_distance_df2 = get_density_df(circRNA_distance_df,"circRNA")$row

f2a_circRNA <- ggplot(circRNA_distance_df2,aes(x = distance))+geom_density(color=ggcolors[2])+
  scale_x_continuous(expand = c(0,0),breaks = c(100,200,300),labels = c("5'UTR","CDS","3'UTR"))+
  theme(
    axis.title = element_text(size = size_axistitle, color = "black"),
    axis.text.x = element_text(size = size_axisgreytext, color = "black",hjust = 1),
    axis.text.y = element_text(size = size_axisblacktext, color = "black"),
    legend.key.size = unit(0.3,"cm"),
    legend.position = "top",
    legend.title= element_blank(),
    legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
    legend.background = element_rect(fill = NA))+scale_y_continuous(expand = c(0,0))

png("distance_raw.png",width=2400, height=1200, res = 300)
f2a_circRNA  
dev.off()


##################density row line
density_plot_df <- data.frame()
for (i in unique(circRNA_distance_df$type) ){
  tmp <- get_density_df(circRNA_distance_df,i)
  tmp <- tmp$stats
  tmp$scale <- scale(tmp$Freq)+abs(min(scale(tmp$Freq)))
  density_plot_df <- rbind(density_plot_df,tmp)
}
density_plot_df$peak_freq <- as.numeric(as.character(density_plot_df$peak_freq))
density_plot_df$scale <- as.numeric(as.character(density_plot_df$scale))
density_plot_df$type <- as.factor(density_plot_df$type)
#
png("distance_rowline.png",width=2400, height=1200, res = 300)
ggplot(density_plot_df, aes(x = peak_freq,y=scale,group=type,color=type))+geom_line()
dev.off()

##################density smooth
density_plot_df2 <- data.frame()
for (i in unique(circRNA_distance_df$type) ){
  tmp <- get_density_df(circRNA_distance_df,i)
  #tmp$scale <- scale(tmp$Freq)+abs(min(scale(tmp$Freq)))
  density_plot_df2 <- rbind(density_plot_df2,tmp$row)
}
density_plot_df2$type <- as.factor(density_plot_df2$type)

png("distance_smooth.png",width=2400, height=1200, res = 300)
ggplot(density_plot_df2, aes(x = distance,group=type,color=type))+
  geom_density()+scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = ggcolors[1:2])+gg_theme
dev.off()
 
pdf("circRNA_feature.pdf")
  all_circ_typepie+ggtitle("Genomic distribution of all circRNAs")
  circ_length.hist
  f2a_circRNA 
  ggplot(density_plot_df, aes(x = peak_freq,y=scale,group=type,color=type))+geom_line()
  ggplot(density_plot_df2, aes(x = distance,group=type,color=type))+
    geom_density()+scale_y_continuous(expand = c(0,0))+
    scale_color_manual(values = ggcolors[1:2])+gg_theme
dev.off()


