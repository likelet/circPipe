library(readr)
library(dplyr)
options(stringsAsFactors=F)
setwd("./")

args <-commandArgs(T)

#testdata = read_delim("/data2/wqj/test/Result/Combination_Matrix/all_tools_merge.matrix",delim = "\t")
testdata = read_delim(args[1],delim = "\t")
testdata = testdata[-nrow(testdata),]
id<-as.character(testdata$id)
chr<-sapply(strsplit(id,"_"),"[",1)
start<-sapply(strsplit(id,"_"),"[",2)
end<-sapply(strsplit(id,"_"),"[",3)
strand<-sapply(strsplit(id,"_"),"[",4)

testdata1 = testdata[,2:ncol(testdata)]
for(i in 1:nrow(testdata1)){
  for (j in 1:ncol(testdata1)){
    if(testdata1[i,j]==1){
      testdata1[i,j]=colnames(testdata1[i,j])
    }else{
      testdata1[i,j]="no"
    }
  }
} 
write.table(testdata1, file = "whatisthis.txt",sep = ",", quote = FALSE, row.names = F)
testdata2=read_delim("whatisthis.txt",delim = "\t")

testdata3=data.frame(chr=chr,start=start,end=end,calculate=testdata2,nothing=".",strand=strand)
write.table(testdata3, file = "modify_merge.bed",sep = "\t", quote = FALSE, row.names = F,col.names = F)
