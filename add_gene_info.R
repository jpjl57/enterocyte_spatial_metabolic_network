add_gene_info<-function(res,ensg2symbol, df){
  ind<-match(res,ensg2symbol$Gene.stable.ID.version)
  df$gene.symbol<-ensg2symbol$Gene.name[ind]
  df$gene.names<-ensg2symbol$Gene.type[ind]
  #res$geneID = ensg2symbol$NCBI.gene.ID
  
  return(df)
}
