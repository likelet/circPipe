######changing the reads count to the order one
#library("readr")
library(dplyr)
args <-commandArgs(T)
#args[1] <- "/home/qijin/perl5/ciri_merge.matrix"
test1 <- as.data.frame(read.table(args[1],header=FALSE),stringsAsFactors = F)
rownames(test1)<-test1$V1
#split the id 
id<-as.character(test1$V1)
chr<-sapply(strsplit(id,"_"),"[",1)
start<-sapply(strsplit(id,"_"),"[",2)
end<-sapply(strsplit(id,"_"),"[",3)
strand<-sapply(strsplit(id,"_"),"[",4)

test1<-test1[,-1]

####change the matrix
df = test1
df.vector <- unlist(df)
uniq_count <- sort(unique(df.vector))
map.df <- data.frame(new_count = 1:length(uniq_count),raw_count = uniq_count)
tmp <- data.frame(raw_mat=df.vector)
tmp$new_df.vector <- map.df$new_count[match(df.vector,map.df$raw_count)]
new_count.df = matrix(tmp$new_df.vector,ncol = ncol(df))
rownames(new_count.df)= rownames(df)
colnames(new_count.df)= colnames(df)

test2<-data.frame(chr,start,end,strand,new_count.df)
write.table(test2, file = args[2],sep = "\t", quote = FALSE, row.names = F, col.names = F)


