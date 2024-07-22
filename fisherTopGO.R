fisherTopGO<-function(genelist,background,category,method){
  #genelist:character vector of interesting genes
  #background:character vector of all genes
  #category :BP, MF, CC
  #method:elimCount or classic
  require(topGO)
  geneList<-factor(as.integer(background %in% genelist))
  names(geneList)<-background
  TopGOdata<-new("topGOdata",ontology=category,allGenes=geneList,nodeSize=5,
                  annot=annFUN.org,mapping="org.Mm.eg.db",ID="symbol")
  resultFisher <- runTest(TopGOdata, algorithm = method, statistic = "fisher")
  return(list(TopGOdata=TopGOdata,resultFisher=resultFisher))
}

