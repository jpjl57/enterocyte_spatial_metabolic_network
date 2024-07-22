setwd("~/Desktop/metabol_recon/mouse_intestine_bulk3/")

hepchol = read.table("Int_deleteModelGenes_optimizeCbModel_for_relative_section1and5.txt", header = F, sep = ",")
colnames(hepchol) = c("NCBI", "HepGRratio", "CholGRratio")

load("~/Downloads/recon1_genes.RData")

ncbi2symbol = function(data, recon1_genes){
  ind = match(data[,1], recon1_genes$entrezgene)
  data$hgnc = recon1_genes$hgnc_symbol[ind]
  data$biotype = recon1_genes$gene_biotype[ind]
  
  result = data
  result = result[complete.cases(result),]
  return(result)
}

hepchol = ncbi2symbol(hepchol, recon1_genes)

hepchol$newHepGR = ifelse(hepchol$HepGRratio<0.1, 0, hepchol$HepGRratio)
hepchol$newCholGR = ifelse(hepchol$CholGRratio<0.1, 0, hepchol$CholGRratio)

diffLHep = hepchol[hepchol$newHepGR<0.1 & hepchol$newCholGR>0.9,]
diffLChol = hepchol[hepchol$CholGRratio<0.1 & hepchol$HepGRratio>0.9,]

diffLHep2 = hepchol[hepchol$newHepGR<0.2 & hepchol$newCholGR>0.8,]
diffLChol2 = hepchol[hepchol$CholGRratio<0.2 & hepchol$HepGRratio>0.8,] #let's use this one

diffLHep3 = hepchol[hepchol$newHepGR<0.25 & hepchol$newCholGR>0.75,]
diffLChol3 = hepchol[hepchol$CholGRratio<0.25 & hepchol$HepGRratio>0.75,]

write.csv(diffLHep2, "int_OptimizeCbModel_diff_lethal_in_relative_section1.csv", row.names = F) #lethal in section1 and not lethal in section 5
write.csv(diffLChol2, "int_OptimizeCbModel_diff_lethal_in_relative_section5.csv", row.names = F) #vice versa
