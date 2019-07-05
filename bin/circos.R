library(readr)
library(ggplot2)
library(circlize)
library(dplyr)
library(RColorBrewer)
options(stringsAsFactors=F)
args <-commandArgs(T)
setwd("./")
######args########
#in-house script
function.script <- args[1] 
#chrName\tlength
##UCSC specific term : hg19/hg38/mm10
species <- args[2]
chrom.size.file <- args[3]
##chr start end name  value(total_reads) strand
final_circ.lst <- args[4] 

##############load data###############
#chromosomes
chrom.size <- read.delim(chrom.size.file,sep  = "\t",header = FALSE)
colnames(chrom.size) = c("chrName","length","offset","linebases","linewidth")
#circRNA final expression matrix
finalmatrix <- read_delim(final_circ.lst,delim = "\t")
colnames(finalmatrix)[1:5] = c("id","chr", "start", "end", "strand")
finalmatrix <- select(finalmatrix,-(id),-(strand))%>%filter( chr%in%chrom.size$chrName)
finalmatrix$value = rowMeans(finalmatrix[,4:ncol(finalmatrix)])
finalmatrix = filter(finalmatrix,value > 0)
finalmatrix$value = log2(finalmatrix$value) ###log2 normalization
circ_plot.df = select(finalmatrix,"chr", "start","end","value")

#####get cytoband file from UCSC or provide path to cytoband file########
####this step for config Ideogram
cytoband.lst = read.cytoband(species = species)
cytoband.df = cytoband.lst$df
#chromesome in cytoband in line with data
cytoband.df = filter(cytoband.df,V1%in%unique(circ_plot.df$chr))

###########plot function############

plot_circos <- function(circ_plot.df,cytoband.df,chrom.size){
  library(circlize)
  #color config
  ggcolors <- c("#BC3C29","#0072B5","#E18727","#20854E","#80397B")
  #angle of plot
  circos.par(start.degree = 90)
  #Ideogram
  circos.initializeWithIdeogram(cytoband = cytoband.df,sort.chr=TRUE,
                                plotType = c("ideogram", "labels"), ideogram.height = 0.05,track.height= 0.1)
  #point plot track for expression level
  count.mean = mean(circ_plot.df$value)
  count.1st = quantile(circ_plot.df$value,0.25)
  count.3st = quantile(circ_plot.df$value,0.75)
  circos.genomicTrackPlotRegion(circ_plot.df, panel.fun = function(region, value, ...) {
    col = ifelse( value > count.mean, ggcolors[1], ggcolors[4])
    circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 20)
    cell.xlim = get.cell.meta.data("cell.xlim")
    for(h in c(min(circ_plot.df$value), count.1st, count.mean, count.3st, max(circ_plot.df$value))) {
      circos.lines(cell.xlim, c(h, h), col = "#00000040")
    }
  }, track.height = 0.25)
  
  #hist plot track for circRNAs density in bins   
  bg.col <- rep(c("#EFEFEF", "#CCCCCC"), length(unique(circ_plot.df$chr)))[1:length(unique(circ_plot.df$chr))]
  circos.trackHist(factors = circ_plot.df$chr,x=(circ_plot.df$start+circ_plot.df$end)/2 ,force.ylim = FALSE,
                   track.height=0.2, bg.col = bg.col, col = ggcolors[2], border = ggcolors[2],bin.size = 3000000)
  circos.clear()
  
}

############plot###################
#png
png("circos.png",width=1800, height=1800, res = 300)

plot_circos(circ_plot.df = circ_plot.df,cytoband.df = cytoband.df , chrom.size = chrom.size)

dev.off()


#pdf
pdf("circos.pdf")

plot_circos(circ_plot.df = circ_plot.df,cytoband.df = cytoband.df , chrom.size = chrom.size)

dev.off()


print("circos plot done !")
