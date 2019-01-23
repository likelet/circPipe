library(dplyr)
library(readr)
options(stringsAsFactors=F)
setwd("./")
args <-commandArgs(T)
mergedata = read_delim(args[1],delim = "\t")
tools_num=ncol(mergedata)-1
kk=apply(mergedata[,2:ncol(mergedata)],1,sum)
mergedata$sum=kk
mergedata=mergedata[mergedata$sum==tools_num,]
write.table(mergedata, file = "all_tools_intersect.matrix",sep = "\t", quote = FALSE, row.names = F)
