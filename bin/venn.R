
#!/usr/bin/env Rscript
# Script to plot venn diagram for used tools
options(stringsAsFactors=F)
#library(getopt)
library(readr)
library(dplyr)
library (VennDiagram)
library(grid)

args <-commandArgs(T)

#################

tools_matrix=args[1]
outprefix=args[2]
#arguments <- matrix(c(
#  'help', 'h', 0, "logical",
#  'tools_matrix' , 't', 2, "character",
#  'outprefix','o',1,'character'
#), ncol=4, byrow=T)



#opt <- getopt(arguments)

#if (!is.null(opt$help) || is.null(opt$tools_matrix) || is.null(opt$outprefix) ) {
#  cat(paste(getopt(opt, usage = T), "\n"))
#  q()
#}


checkReadable <- function(filename) {
  res <- file.access(names=filename, mode=4) == 0
  if (!res) {
    warning(paste(filename, "is not readable", sep=" "))
  }
  res
}



tools_matrix = tools_matrix 
outprefix = outprefix
checkReadable(tools_matrix)

print(tools_matrix)
print(outprefix)
#load data
venndata <- as.data.frame(read_delim(tools_matrix,delim = "\t"))

# ncol
if (ncol(venndata) <= 2 ){
  stop("there is just one tools, cannot plot Venn Diagram.")
  
}else if (ncol(venndata) > 6 ){
  stop("Too many tools , VennDiagram package don't support it.")
}

#convert df to list for Venn.plot
total_circ.lst = list()
for (i in 2:ncol(venndata)){
  tool_name = colnames(venndata)[i]
  circ.lst = as.character(venndata[which(venndata[,i] > 0),]$id)
  total_circ.lst[[(i-1)]] = circ.lst
} 

tools_number = ncol(venndata) - 1
names(total_circ.lst) = colnames(venndata)[2:ncol(venndata)]
#str(total_circ.lst)
#plot

ggcolors <- c("#BC3C29","#0072B5","#E18727","#20854E","#80397B")#pal_npg("nrc")(5)#gg_color_hue(5)

venn.plot <- venn.diagram(x=total_circ.lst[-2],filename = NULL,scaled = TRUE,
                          col = rep("white",tools_number),fill= ggcolors[1:tools_number],
                          cat.cex = 1.5,cex = 1.5,
                          alpha=rep(0.6,tools_number),resolution = 300)


#output venn plot
png(paste0(outprefix,".png"), width=1800, height=1800, res = 300)
grid.draw(venn.plot)
dev.off()

pdf(file = paste0(outprefix,".pdf"))
grid.draw(venn.plot)
dev.off()

