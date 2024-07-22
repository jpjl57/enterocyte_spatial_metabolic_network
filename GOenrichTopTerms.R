GOenrichTopTerms<-function(fisherTopGOresult,topK,max_gene,max_char,filename){
  maxNodes<-length(fisherTopGOresult$resultFisher@score)
  x<-GenTable(fisherTopGOresult$TopGOdata,fisherTopGOresult$resultFisher,topNodes=maxNodes,numChar=max_char)
  x<-x[x[,3]<max_gene,] #GO terms with at least 10 genes (implemented in topGO) and at most 500 genes
  x$result1[grep("<",x$result1)]<-1e-30
  x$result1<-as.numeric(x$result1)
  x$FDR<-p.adjust(x$result1,method="BH")
  top<-x[1:topK,]
  write.table(top,file=filename,sep="\t",quote=F,col.names=T,row.names=F)
}
