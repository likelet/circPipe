library("readr")
library(dplyr)
library (VennDiagram)
library(grid)
options(stringsAsFactors=F)
args <-commandArgs(T)
#args[1]<-"/home/wqj/code/circPipe-master/all_tools_merge.matrix"
venndata <- as.data.frame(read_delim(args[1],delim = "\t"))

#hehe=nrow(venndata[(venndata$find_circ==1&venndata$segemehl==1&venndata$ciri==1&venndata$circexplorer2==1),])

venn.plot <- draw.quintuple.venn(
  area1 = venndata[nrow(venndata),2],
  area2 = venndata[nrow(venndata),3],
  area3 = venndata[nrow(venndata),4],
  area4 = venndata[nrow(venndata),5],
  area5 = venndata[nrow(venndata),6],
  n12 = nrow(venndata[(venndata$find_circ==1&venndata$circexplorer2==1),]),
  n13 = nrow(venndata[(venndata$find_circ==1&venndata$ciri==1),]),
  n14 = nrow(venndata[(venndata$find_circ==1&venndata$mapsplice==1),]),
  n15 = nrow(venndata[(venndata$find_circ==1&venndata$segemehl==1),]),
  n23 = nrow(venndata[(venndata$circexplorer2==1&venndata$ciri==1),]),
  n24 = nrow(venndata[(venndata$circexplorer2==1&venndata$mapsplice==1),]),
  n25 = nrow(venndata[(venndata$circexplorer2==1&venndata$segemehl==1),]),
  n34 = nrow(venndata[(venndata$ciri==1&venndata$mapsplice==1),]),
  n35 = nrow(venndata[(venndata$ciri==1&venndata$segemehl==1),]),
  n45 = nrow(venndata[(venndata$mapsplice==1&venndata$segemehl==1),]),
  n123 = nrow(venndata[(venndata$find_circ==1&venndata$circexplorer2==1&venndata$ciri==1),]),
  n124 = nrow(venndata[(venndata$find_circ==1&venndata$circexplorer2==1&venndata$mapsplice==1),]),
  n125 = nrow(venndata[(venndata$find_circ==1&venndata$circexplorer2==1&venndata$segemehl==1),]),
  n134 = nrow(venndata[(venndata$find_circ==1&venndata$mapsplice==1&venndata$ciri==1),]),
  n135 = nrow(venndata[(venndata$find_circ==1&venndata$segemehl==1&venndata$ciri==1),]),
  n145 = nrow(venndata[(venndata$find_circ==1&venndata$mapsplice==1&venndata$segemehl==1),]),
  n234 = nrow(venndata[(venndata$circexplorer2==1&venndata$ciri==1&venndata$mapsplice==1),]),
  n235 = nrow(venndata[(venndata$circexplorer2==1&venndata$ciri==1&venndata$segemehl==1),]),
  n245 = nrow(venndata[(venndata$circexplorer2==1&venndata$mapsplice==1&venndata$segemehl==1),]),
  n345 = nrow(venndata[(venndata$mapsplice==1&venndata$ciri==1&venndata$segemehl==1),]),
  n1234 = nrow(venndata[(venndata$find_circ==1&venndata$circexplorer2==1&venndata$ciri==1&venndata$mapsplice==1),]),
  n1235 = nrow(venndata[(venndata$find_circ==1&venndata$circexplorer2==1&venndata$ciri==1&venndata$segemehl==1),]),
  n1245 = nrow(venndata[(venndata$find_circ==1&venndata$circexplorer2==1&venndata$segemehl==1&venndata$mapsplice==1),]),
  n1345 = nrow(venndata[(venndata$find_circ==1&venndata$segemehl==1&venndata$ciri==1&venndata$mapsplice==1),]),
  n2345 = nrow(venndata[(venndata$segemehl==1&venndata$circexplorer2==1&venndata$ciri==1&venndata$mapsplice==1),]),
  n12345 = nrow(venndata[(venndata$find_circ==1&venndata$circexplorer2==1&venndata$ciri==1&venndata$mapsplice==1&venndata$segemehl==1),]),
  category = c(paste0("Find_circ", "\n(n=", venndata[nrow(venndata),2],")", sep=""), paste0("CIRCexplorer2", "\n(n=", venndata[nrow(venndata),3],")", sep=""), paste0("CIRI", "\n(n=", venndata[nrow(venndata),4],")", sep=""),paste0("Mapsplice", "\n(n=", venndata[nrow(venndata),5],")", sep="") , paste0("Segemehl", "\n(n=", venndata[nrow(venndata),6],")", sep="")),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1.4,
  lwd = c(1,1,1,1,1),
  cat.pos = c(0,-20,-160,180,0),
  cat.dist = c(0.23,0.23,0.2,0.22,0.23),
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
)
png(args[2], width=1800, height=1800, res = 300)
grid.draw(venn.plot)
dev.off()

