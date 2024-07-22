setwd("~/Desktop/practice/Mouse_intestines_bulk/")
dat.files  <- list.files(path="~/Desktop/practice/Mouse_intestines_bulk/exprssion",
                         recursive=T,
                         pattern="mouse"
                         ,full.names=T)


lst = lapply(dat.files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
lst1 = lapply(lst, "[", c("target_id", "tpm"))


#order components in list
for (i in 1:length(lst1)){
  lst2 = order(lst1[[i]]$target_id)
  lst1[[i]] = lst1[[i]][lst2, ]
}

temp = do.call(cbind, lst1)
rn = temp$target_id
col_temp = which(colnames(temp) %in% "target_id")
df = temp[, -col_temp]
master_file = data.frame(rn, df)
rownames(master_file) = NULL
sum(duplicated(master_file$rn))

df = read.csv("~/Downloads/mart_export.txt")
ind = match(master_file$rn, df$Transcript.stable.ID)
res = df$Gene.name[ind]
master_file$rn = res
sum(duplicated(master_file$rn))
master_file = master_file[complete.cases(master_file),]

library(dplyr)
master_file = master_file %>% group_by(rn) %>% summarise(tpm = sum(tpm), tpm.1 = sum(tpm.1), tpm.2 = sum(tpm.2), tpm.3 = sum(tpm.3),
                                                         tpm.4 = sum(tpm.4),tpm.5 = sum(tpm.5),tpm.6 = sum(tpm.6),tpm.7 = sum(tpm.7),
                                                         tpm.8 = sum(tpm.8),tpm.9 = sum(tpm.9),tpm.10 = sum(tpm.10),tpm.11 = sum(tpm.11)
                                                         ,tpm.12 = sum(tpm.12),tpm.13 = sum(tpm.13),tpm.14 = sum(tpm.14))

colnames(master_file) = c("rn", "m1s1", "m1s2", "m1s3", "m1s4", "m1s5",
                               "m2s1", "m2s2", "m2s3", "m2s4", "m2s5",
                               "m3s1", "m3s2", "m3s3", "m3s4", "m3s5")
master_file = as.data.frame(master_file)



mf = master_file
rownames(mf) = mf$rn
background = rownames(mf)
mf = mf[,-1]


#filter
mf = mf[rowSums(mf)>ncol(mf),] #filter
mf = mf[apply(mf, 1, var) > 1,]

#rownames(master_file) = master_file$rn
#master_file = master_file[,-1]
#save(master_file, file = "unfiltered_master_file.RData")
#save(mf, file = "filtered_master_file.RData")

########################################################################################################################################
#Clustering of global gene expression (not using differentially expresseed) for general pattern
########################################################################################################################################
library(ComplexHeatmap)
library(circlize)
library(R.utils)

colnames(mf)= capitalize(tolower(colnames(mf)))

msf = as.matrix(t(scale(t(mf),center=T)))
msf = na.omit(msf)

#Adjust row_km for number of clusters
p = Heatmap(msf, col = colorRamp2(c(-4, 0, 4), c("deepskyblue1", "black", "firebrick1")), cluster_columns = T, name = "Z-Score", 
            show_row_names = F , show_column_names = T, show_row_dend = F, row_km = 2,
            row_title_gp = gpar(fontsize = 12, fontface = "bold"), cluster_row_slices = T,row_title = "Cluster %s", use_raster = F)

hm = ComplexHeatmap::draw(p, merge_legends = TRUE, heatmap_legend_side = "right")

#check gene categories
load("ensg2symbol.RData")

gene_category_info<-function(res,ensg2symbol){
  ind<-match(res,ensg2symbol$Gene.name)
  gene.symbol<-ensg2symbol$Gene.name[ind]
  gene.names<-ensg2symbol$Gene.type[ind]
  gene.description<-ensg2symbol$WikiGene.description[ind]
  cate = cbind(gene.symbol, gene.names, gene.description)
  #res$geneID = ensg2symbol$NCBI.gene.ID
  
  
  
  return(as.data.frame(cate))
}

categ = gene_category_info(rownames(mf), ensg2symbol) #output = genes with categories
categ = categ[complete.cases(categ),]
summarize_category = function(categ){
  lv = levels(factor(categ$gene.names))
  df = data.frame("Numbers")
  for (i in 1:length(lv)){
    num = sum(categ$gene.names==lv[i])
    df = cbind(df, num)}
  
  df = df[,-1]
  colnames(df) = lv
  return(df)
}




df = summarize_category(categ) #counts number of genes per category

#plot number of genes in each category
library(ggplot2)
df = t(df)
ndf = data.frame(rownames(df), df)
rownames(ndf) = NULL
colnames(ndf) = c("category", "number")
ndf$category = gsub("_", " ", ndf$category)

(bp =ggplot(data = ndf, aes(x=ndf[,1], y=ndf[,2])) + 
    geom_bar(stat='identity')+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
          axis.text.y=element_text(face = "bold",size=13),
          axis.title.x = element_text(face = "bold",size = 15), axis.title.y = element_text(face = "bold",size=15),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    xlab("Gene category") + ylab("Number of genes"))

tiff("number_of_genes_in_heatmap.tiff", width = 9, height = 6, units = "in", res = 300)
bp
dev.off()

#extract genes in gene cluster
set.seed(123456789)  #For reproduciblity


r.dend <- row_dend(hm)  #Extract row dendrogram
rcl.list <- row_order(hm)  #Extract clusters (output is a list)

lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

# loop to extract genes for each cluster.
for (i in 1:length(row_order(hm))){
  if (i == 1) {
    set.seed(123456789)
    clu <- row.names(msf[row_order(hm)[[i]],])
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Heatmap_Cluster")
  } else {
    set.seed(123456789)
    clu <- row.names(msf[row_order(hm)[[i]],])
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}


#check
library(openxlsx)
out = as.data.frame(out)
write.xlsx(out, "genes_in_cluster.xlsx")


#enrich genes in each cluster
library(dplyr)
library(R.utils)
library(stringr)
library(openxlsx)
library(ggplot2)
library(topGO)
#load the functions from working directory (I put them in the same directory); if you want to use the human annotation, open the function and change "org.Mm.eg.db" to "org.Hs.eg.db" after installing.
source("fisherTopGO.R")
source("GOenrichTopTerms.R")

#Gene Ontology enrichment
cluster1 = out[grep("cluster1", out$Heatmap_Cluster),]
cluster2 = out[grep("cluster2", out$Heatmap_Cluster),]
#cluster3 = out[grep("cluster3", out$Heatmap_Cluster),]
#cluster4 = out[grep("cluster4", out$Heatmap_Cluster),]

lst = list(cluster1, cluster2)#, cluster3)#, cluster4)
for (i in 1:length(lst)){
  if (length(unique(subset(lst[[i]], select="GeneID",drop=T))) > 0){
    bp_up<-fisherTopGO(unique(subset(lst[[i]], select="GeneID",drop=T)),
                       unique(rownames(mf)),"BP","classic") #BP > MF or CC
    filename<-paste("top_GO_enrich_in_cluster",i,".xlsx",sep="")
    GOenrichTopTerms(bp_up,topK=500,max_gene=500,max_char=100,filename=filename)
    go<-read.table(filename,header=T,sep="\t",stringsAsFactors = F,comment.char = "",quote="")
    go$Term<-capitalize(tolower(go$Term))
    go$Term<-str_wrap(go$Term,width=38)
    go$Term<-factor(go$Term,levels=rev(go$Term),ordered=T)
    p2 = ggplot(go[1:20,], aes(x=Term, y=-log10(FDR),fill=-log10(FDR))) + 
      geom_bar(stat='identity') + coord_flip() +
      xlab("GO Terms") + ylab("-log10(FDR)") +
      scale_fill_gradient(low="yellow",high="red",name="-log10(FDR)") + 
      theme(axis.text.x = element_text(size=10),axis.text.y=element_text(size=15),
            axis.title=element_text(size=15),legend.position=c(1,0),legend.justification=c(1,0),
            legend.title=element_blank())+
      geom_hline(yintercept = 1)
    png(paste("top_GO_enrich_in_cluster",i,".png",sep=""),height=900,width=600)
    print(p2)
    dev.off()}
  else{
    next
  }
}


#Take a look at the plots showing the top 20 enriched GO terms
#You will need to manually select representative GO terms and type in a few below
#as text1, text2, text3, etc. 
#the number of "text" variable that you will need to create will depend on how manu clusters you made in the above steps

text1 = c("Regulation of cytosolic calcium ion", "Sensory organ development",
          "Gamma-aminobutyric acid signaling pathway", "Taxis",
          "Extracellular matrix organization")
text2 = c("Nucleoside and nucleotide metabolism", 'Cofactor metabolism', "Translation",
          "Mitochondrion organization", "Cellular respiration", "ATP metabolism")

text_list = list(text1 = text1, text2 = text2)


# note how we set the width of this empty annotation
ra = rowAnnotation(foo = anno_empty(border = FALSE,
                                    width = max_text_width(unlist(text_list)) + unit(3, "mm")), name = "Enriched GO", 
                   show_annotation_name = T)

annot_all<- data.frame(Condition = c("Section 1", "Section 2", "Section 3", "Section 4", "Section 5", 
                                     "Section 1", "Section 2", "Section 3", "Section 4", "Section 5",
                                     "Section 1", "Section 2", "Section 3", "Section 4", "Section 5"))
row.names(annot_all) = NULL

#Here, set your group names and assign colors
head(msf)
colnames(msf) = c("Section 1", "Section 2", "Section 3", "Section 4", "Section 5", 
                  "Section 1", "Section 2", "Section 3", "Section 4", "Section 5",
                  "Section 1", "Section 2", "Section 3", "Section 4", "Section 5")
col <- list(Condition = c("Section 1" = "red", "Section 2" ="lightsalmon",
                         "Section 3" = "yellow", "Section 4" = "deepskyblue",
                         "Section 5" = "navy"))
ha<- HeatmapAnnotation(df = annot_all, col = col, show_annotation_name = F)

p = Heatmap(msf, col = colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")), top_annotation = ha, cluster_columns = T, name = "Z-Score", 
            show_row_names = F , show_column_names = F, show_row_dend = F, row_km = 2, right_annotation = ra,
            row_title_gp = gpar(fontsize = 12, fontface = "bold"), cluster_row_slices = T,row_title = "Cluster %s", use_raster = F)

hm = ComplexHeatmap::draw(p, merge_legends = TRUE, heatmap_legend_side = "right")
for(i in 1:2) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(0.2, "mm"), gp = gpar(fill = i, col = "black"), just = "left")
    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(1, "mm"),gp = gpar(fontsize = 11, fontface = "bold"), just = "left")
  })
}

#use dev.copy
dev.copy(tiff, "Intestine_everything_heatmap_using_devcopy.tiff", width = 9, height = 7,units = "in", res = 300)
dev.off()

########################################################################################################################################
#Do PCA
########################################################################################################################################
set.seed(123456)
msf = as.matrix(t(scale(t(mf), scale = T)))
pca<-princomp(na.omit(msf))
sd<-pca$sdev
loadings<-pca$loadings
var<-sd^2
var.percent<-var/sum(var)*100
scores <- data.frame(sample.groups=annot_all,
                     pca$loadings[,1:3])

colnames(scores)<-c("sample.groups", "PC1","PC2","PC3")

#scores$PC1<- -(scores$PC1)
#scores$PC2<- -(scores$PC2)
min_axis<-min(scores[,2:4])
max_axis<-max(scores[,2:4])
library(grid)
library(ggplot2)
cbPalette<-c("azure4","black","blue","brown","cadetblue","chartreuse","cyan",
             "darkorange","darkorchid","deeppink","gold","lightcoral",
             "lightseagreen","magenta","red","lightsalmon","yellow2","mediumorchid4",
             "deepskyblue","mediumvioletred","olivedrab","cornsilk","lavender","navajowhite4", "navy", "yellowgreen")


(p1=ggplot(data=scores, aes(x=PC1,y=PC2,color=sample.groups))+ 
    geom_point(size = 8, aes(shape = sample.groups)) +
    scale_shape_manual(values = c(16, 17,20, 15,18))+
    #geom_point(size = 0.1, stroke = 0, shape = 16) +
    theme(legend.key = element_rect(fill = "white"), legend.key.width=unit(0.5,"cm"),legend.key.height=unit(0.3,"cm"),
          legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
          legend.text=element_text(size=15, face = "bold"),axis.text.x = element_text(size=15),
          axis.text.y=element_text(size=15),axis.title=element_text(face = "bold",size=20),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    guides(color=guide_legend(ncol=5)) +
    xlim(min(scores$PC1),max(scores$PC1)) + 
    ylim(min(scores$PC2),max(scores$PC2)) + 
    scale_color_manual(values=c("dodgerblue4", "steelblue", "lightsteelblue4", "lightgoldenrod3", "yellow2")) + 
    xlab(paste("PC 1 (",round(var.percent[1],1),"%)",sep="")) + 
    ylab(paste("PC 2 (",round(var.percent[2],1),"%)",sep="")) )

tiff("Intestine_PCA_all_groups.tiff",units="in", width = 8, height = 5, res=300)
p1
dev.off()

########################################################################################################################################
#Get metabolic genes
########################################################################################################################################
s1 = rowMeans(mf[,c(1, 6, 11)])
s2 = rowMeans(mf[,c(2, 7, 12)])
s3 = rowMeans(mf[,c(3, 8, 13)])
s4 = rowMeans(mf[,c(4, 9, 14)])
s5 = rowMeans(mf[,c(5, 10, 15)])

s1 = data.frame(rn = rownames(mf), mean_expr = s1)
s2 = data.frame(rn = rownames(mf), mean_expr = s2)
s3 = data.frame(rn = rownames(mf), mean_expr = s3)
s4 = data.frame(rn = rownames(mf), mean_expr = s4)
s5 = data.frame(rn = rownames(mf), mean_expr = s5)

rownames(s1) = NULL
rownames(s2) = NULL
rownames(s3) = NULL
rownames(s4) = NULL
rownames(s5) = NULL

clust = list(s1, s2, s3, s4, s5)
write.xlsx(clust, "Mouse_intestine_bulk_avg_expr.xlsx")

setwd("~/Desktop/practice/Mouse_intestines_bulk/")
library(openxlsx)

clust = list()
clust[[1]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 1)
clust[[2]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 2)
clust[[3]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 3)
clust[[4]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 4)
clust[[5]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 5)

load("iMM1415.RData")

iMM1415 = iMM1415[!duplicated(iMM1415$NCBI.gene.ID),]

for (i in 1:5){
  ind = match(iMM1415$Gene.name,clust[[i]][,1])
  clust[[i]] = clust[[i]][ind,]
  clust[[i]] = clust[[i]][complete.cases(clust[[i]][,1]),]
  
  #clust[[i]] = clust[[i]][,colSums2(as.matrix(clust[[i]][,2:ncol(clust[[i]])]))>0.000001]
  #rowSums2(as.matrix(clust[[i]][,2:ncol(clust[[i]])]))
}

head(clust[[1]])
dim(clust[[1]])

symbol2ncbi = function(data, iMM1415){
  ind = match(data[,1], iMM1415$Gene.name)
  res = iMM1415$NCBI.gene.ID[ind]
  
  result = data.frame('genes' = res, 'state' = data[,2])
  result = result[complete.cases(result),]
  return(result)
}


clustv2 = list()
for (i in 1:5) {
  clustv2[[i]] = symbol2ncbi(clust[[i]], iMM1415)
}

dim(clustv2[[1]])
dim(clustv2[[2]])
dim(clustv2[[3]])
dim(clustv2[[4]])
dim(clustv2[[5]])

#add unmatched recon1 genes 
#expression frequency set to 0 for those

add = data.frame(genes = iMM1415$NCBI.gene.ID[!iMM1415$NCBI.gene.ID %in% clustv2[[1]]$genes],
                 state = rep(0,length(iMM1415$NCBI.gene.ID[!iMM1415$NCBI.gene.ID %in% clustv2[[1]]$genes])))

for (i in 1:5) {
  clustv2[[i]] = as.data.frame(dplyr::bind_rows(clustv2[[i]], add))
}

sum(clustv2[[1]]$genes%in%iMM1415$NCBI.gene.ID)
sum(clustv2[[2]]$genes%in%iMM1415$NCBI.gene.ID)
sum(clustv2[[3]]$genes%in%iMM1415$NCBI.gene.ID)

dim(clustv2[[1]][clustv2[[1]]$state>20,])
dim(clustv2[[2]][clustv2[[2]]$state>20,])
dim(clustv2[[3]][clustv2[[3]]$state>20,])
dim(clustv2[[4]][clustv2[[4]]$state>20,])
dim(clustv2[[5]][clustv2[[5]]$state>20,])

cluster.ids <- c("Section_1", "Section_2", "Section_3", "Section_4", "Section_5")

for (i in 1:5) {
  write.table(clustv2[[i]], file = paste("Intestine_bulk_expr_", cluster.ids[i], ".txt", sep = ""), 
              quote = F, sep = "\t", row.names = F, col.names = T)
}


########################################################################################################################################
#Get metabolic genes after converting to hgnc
########################################################################################################################################
s1 = rowMeans(mf[,c(1, 6, 11)])
s2 = rowMeans(mf[,c(2, 7, 12)])
s3 = rowMeans(mf[,c(3, 8, 13)])
s4 = rowMeans(mf[,c(4, 9, 14)])
s5 = rowMeans(mf[,c(5, 10, 15)])

s1 = data.frame(rn = rownames(mf), mean_expr = s1)
s2 = data.frame(rn = rownames(mf), mean_expr = s2)
s3 = data.frame(rn = rownames(mf), mean_expr = s3)
s4 = data.frame(rn = rownames(mf), mean_expr = s4)
s5 = data.frame(rn = rownames(mf), mean_expr = s5)

rownames(s1) = NULL
rownames(s2) = NULL
rownames(s3) = NULL
rownames(s4) = NULL
rownames(s5) = NULL

clust = list(s1, s2, s3, s4, s5)
#write.xlsx(clust, "Mouse_intestine_bulk_avg_expr.xlsx")

setwd("~/Desktop/practice/Mouse_intestines_bulk/")
library(openxlsx)

clust = list()
clust[[1]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 1)
clust[[2]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 2)
clust[[3]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 3)
clust[[4]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 4)
clust[[5]] = read.xlsx("Mouse_intestine_bulk_avg_expr.xlsx", sheet = 5)

load("human2mouse.RData")
load("~/Downloads/recon1_genes.RData")
recon1_genes = recon1_genes[!duplicated(recon1_genes$hgnc_symbol),]

for (i in 1:5){
  ind = match(clust[[i]][,1], human2mouse$mgi)
  res = human2mouse$hgnc[ind]
  clust[[i]]$rn = res
  clust[[i]] = clust[[i]][complete.cases(clust[[i]]$rn),]
  
  ind = match(recon1_genes$hgnc_symbol,clust[[i]]$rn)
  clust[[i]] = clust[[i]][ind,]
  clust[[i]] = clust[[i]][complete.cases(clust[[i]][,1]),]
  
  #clust[[i]] = clust[[i]][,colSums2(as.matrix(clust[[i]][,2:ncol(clust[[i]])]))>0.000001]
  #rowSums2(as.matrix(clust[[i]][,2:ncol(clust[[i]])]))
}

#compare recon1 with homolog result
load("iMM1415.RData")
mgenes = clust[[1]]
iMM1415 = iMM1415[!duplicated(iMM1415$Gene.name),]

ind = match(human2mouse$hgnc, recon1_genes$hgnc_symbol)
met = recon1_genes[ind,]
met = met[complete.cases(met$hgnc_symbol),]
met = met[order(met$entrezgene),]
met = met[!duplicated(met$entrezgene),]
dim(met)
dim(recon1_genes)

head(clust[[1]])
dim(clust[[1]])

setdiff(met$hgnc_symbol, clust[[i]][,1])


ind = match(recon1_genes$hgnc_symbol, human2mouse$hgnc)
res = human2mouse[ind,]
res = res[complete.cases(res),]

ind = match(iMM1415$Gene.name, mgenes$rn)
mgenes = mgenes[ind,]
mgenes = mgenes[complete.cases(mgenes$rn),]

ind = match(res$hgnc, mgenes$rn)
mmet = mgenes[!ind,]


symbol2ncbi = function(data, recon1_genes){
  ind = match(data[,1], recon1_genes$hgnc_symbol)
  res = recon1_genes$entrezgene[ind]
  
  result = data.frame('genes' = res, 'state' = data[,2])
  result = result[complete.cases(result),]
  return(result)
}


clustv2 = list()
for (i in 1:5) {
  clustv2[[i]] = symbol2ncbi(clust[[i]], recon1_genes)
}

dim(clustv2[[1]])
dim(clustv2[[2]])
dim(clustv2[[3]])
dim(clustv2[[4]])
dim(clustv2[[5]])

#add unmatched recon1 genes 
#expression frequency set to 0 for those

add = data.frame(genes = recon1_genes$entrezgene[!recon1_genes$entrezgene %in% clustv2[[1]]$genes],
                 state = rep(0,length(recon1_genes$entrezgene[!recon1_genes$entrezgene %in% clustv2[[1]]$genes])))

for (i in 1:5) {
  clustv2[[i]] = as.data.frame(dplyr::bind_rows(clustv2[[i]], add))
}

sum(clustv2[[1]]$genes%in%recon1_genes$entrezgene)
sum(clustv2[[2]]$genes%in%recon1_genes$entrezgene)
sum(clustv2[[3]]$genes%in%recon1_genes$entrezgene)

dim(clustv2[[1]][clustv2[[1]]$state>10,])
dim(clustv2[[2]][clustv2[[2]]$state>10,])
dim(clustv2[[3]][clustv2[[3]]$state>10,])
dim(clustv2[[4]][clustv2[[4]]$state>10,])
dim(clustv2[[5]][clustv2[[5]]$state>10,])


cluster.ids <- c("Section_1", "Section_2", "Section_3", "Section_4", "Section_5")

for (i in 1:5) {
  write.table(clustv2[[i]], file = paste("Intestine_bulk_expr_mouse2human", cluster.ids[i], ".txt", sep = ""), 
              quote = F, sep = "\t", row.names = F, col.names = T)
}

#absolute and relative changes
clustv3 = list()
state5 = ((clustv2[[5]]$state+1)/4)-((clustv2[[1]]$state+1)/(clustv2[[5]]$state+1))
state1 = ((clustv2[[1]]$state+1)/4)-((clustv2[[5]]$state+1)/(clustv2[[1]]$state+1))

clustv3[[1]] = data.frame(genes = clustv2[[1]]$genes, state = state1)
clustv3[[5]] = data.frame(genes = clustv2[[5]]$genes, state = state5)

dim(clustv3[[1]][clustv3[[1]]$state>4,])
dim(clustv3[[5]][clustv3[[5]]$state>4,])

clustv3[[1]]$state[clustv3[[1]]$state>4] == clustv3[[5]]$state[clustv3[[5]]$state>4]

write.table(clustv3[[1]], file = "Intestine_bulk_expr_mouse2human_relative_section1.txt", 
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(clustv3[[5]], file = "Intestine_bulk_expr_mouse2human_relative_section5.txt", 
            quote = F, sep = "\t", row.names = F, col.names = T)


########################################################################################################################################
#Run differential expression analysis
########################################################################################################################################
library(dplyr)
library(R.utils)
library(stringr)
library(openxlsx)
library(ggplot2)
library(topGO)
source("fisherTopGO.R")
source("GOenrichTopTerms.R")

#Differential expression, Volcano plot, Gene Ontology 
setwd("~/Desktop/practice/Mouse_intestines_bulk/")
dat.files  <- list.files(path="~/Desktop/practice/Mouse_intestines_bulk/exprssion/",
                         recursive=T,
                         pattern="mouse"
                         ,full.names=T)


lst = lapply(dat.files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
lst1 = lapply(lst, "[", c("target_id", "est_counts"))


#order components in list
for (i in 1:length(lst1)){
  lst2 = order(lst1[[i]]$target_id)
  lst1[[i]] = lst1[[i]][lst2, ]
}

temp = do.call(cbind, lst1)
rn = temp$target_id
col_temp = which(colnames(temp) %in% "target_id")
df = temp[, -col_temp]
master_file = data.frame(rn, df)
rownames(master_file) = NULL
sum(duplicated(master_file$rn))

df = read.csv("~/Downloads/mart_export.txt")
ind = match(master_file$rn, df$Transcript.stable.ID)
res = df$Gene.stable.ID.version[ind]
master_file$rn = res
sum(duplicated(master_file$rn))
master_file = master_file[complete.cases(master_file),]
background = df$Gene.name[ind]

library(dplyr)
master_file = master_file %>% group_by(rn) %>% summarise(est_counts = sum(est_counts), est_counts.1 = sum(est_counts.1), est_counts.2 = sum(est_counts.2), est_counts.3 = sum(est_counts.3),
                                                         est_counts.4 = sum(est_counts.4),est_counts.5 = sum(est_counts.5),est_counts.6 = sum(est_counts.6),est_counts.7 = sum(est_counts.7),
                                                         est_counts.8 = sum(est_counts.8),est_counts.9 = sum(est_counts.9),est_counts.10 = sum(est_counts.10),est_counts.11 = sum(est_counts.11)
                                                         ,est_counts.12 = sum(est_counts.12),est_counts.13 = sum(est_counts.13),est_counts.14 = sum(est_counts.14))

colnames(master_file) = c("rn", "m1s1", "m1s2", "m1s3", "m1s4", "m1s5",
                          "m2s1", "m2s2", "m2s3", "m2s4", "m2s5",
                          "m3s1", "m3s2", "m3s3", "m3s4", "m3s5")
meta = c("Section_1", "Section_2", "Section_3", "Section_4", "Section_5", 
         "Section_1", "Section_2", "Section_3", "Section_4", "Section_5",
         "Section_1", "Section_2", "Section_3", "Section_4", "Section_5")
meta = data.frame(cols = colnames(master_file)[-1], meta = meta)
#meta$meta = factor(meta$meta, levels = c("Section_1", "Section_2", "Section_3", "Section_4", "Section_5"), ordered = T)
master_file = as.data.frame(master_file)



mf = master_file
rownames(mf) = master_file$rn
mf = mf[,-1]

#filter
#mf = mf[rowSums(mf)>ncol(mf),] #filter
#mf = mf[apply(mf, 1, var) > 1,]

lines_i<-c("Section_1")
lines_j<-c("Section_2", "Section_3", "Section_4", "Section_5")

inflamm.l = list()
trp.l = list()
oxi.l = list()
epi.l = list()
ph1.l = list()
ph2.l = list()
nr.l = list()
ba.l = list()
tf.l = list()
int.l = list()
load("phaseI.RData")
load("phaseII.RData")
load("oxi_stress.RData")
load("epi_modification.RData")
load("inflamm.RData")
load("nuclear_receptors.RData")
load("bile_acid.RData")
load("transporters.RData")
load("TF_drug_metabolism.RData")
int = read.xlsx("~/Downloads/Mammalian_Metabolic_Final.xlsx")
int = int$Gene.Symbol

return_match = function(match_list,DE_genelist){
  ind = match(match_list,DE_genelist[,7])
  match_res = DE_genelist[ind, ]
  match_res[complete.cases(match_res), ]
}


mf = round(mf)
source("add_gene_info.R")
library(DESeq2)
for (i in seq_along(lines_j)){
  #DEseq
  countData = as.data.frame(c(mf[,meta$meta %in% c(lines_i)],
                              mf[,meta$meta %in% c(lines_j[i])]))
  rownames(countData) = rownames(mf)
  colData = bind_rows(meta[meta$meta %in% c(lines_i),],
                      meta[meta$meta %in% c(lines_j[i]),])
  dds<-DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~0+meta)
  dds$meta = factor(colData$meta, levels = c(paste(lines_i), paste(lines_j[i])))
  dds = DESeq(dds)
  res = results(dds)
  res = res[complete.cases(res),]
  result = as.data.frame(res[complete.cases(res),])
  result<-add_gene_info(result,df)
  write.xlsx(result, paste("DESeq2_all_genes",lines_i,"_vs_",lines_j[i],".xlsx",sep = ""))
  sig = result[result$padj<0.1,]
  write.xlsx(sig, paste("DESeq2_sig_all_genes",lines_i,"_vs_",lines_j[i],".xlsx",sep = ""))
  
  #match with genes of interest
  
  inflamm.l[[i]] = return_match(inflamm_gene_ls$X1, sig)
  trp.l[[i]] = return_match(tr.unique, sig)
  oxi.l[[i]] = return_match(oxi, sig)
  epi.l[[i]] = return_match(epi, sig)
  ph1.l[[i]] = return_match(ph1, sig)
  ph2.l[[i]] = return_match(ph2, sig)
  nr.l[[i]] = return_match(nr, sig)
  ba.l[[i]] = return_match(ba, sig)
  tf.l[[i]] = return_match(tf, sig)
  int.l[[i]] = return_match(int,sig)
  #do PCA
  #vsd = vst(dds, blind = T)
  #pcaData = plotPCA(vsd, intgroup = c("Treatment"), returnData = T)
  #percentVar = round(100*attr(pcaData, "percentVar"))
  #p=ggplot(pcaData, aes(x=PC1,y=PC2,color=group))+ 
  #  geom_point(size = 3) +
  #  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
  #        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
  #        legend.text=element_text(size=20),axis.text.x = element_text(size=12),
  #        axis.text.y=element_text(size=12),axis.title=element_text(size=20),legend.title=element_blank(), 
  #        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  #  guides(color=guide_legend(ncol=10)) +
  #  xlim(min(pcaData$PC1),max(pcaData$PC1)) + 
  #  ylim(min(pcaData$PC2),max(pcaData$PC2)) + 
  #  scale_color_manual(values= c("darkred", "lightblue")) + 
  #  xlab(paste0("PC1: ", percentVar[1],"% variance"))+
  #  ylab(paste0("PC2: ", percentVar[2],"% variance"))
  #geom_text(aes(label=group), size=3)
  #png(paste("PCA",lines_i[i],"_",lines_j,".png", sep = ""), width = 1200, height = 600)
  #print(p)
  #dev.off()
  #Volcano plot
  result$score<- -log10(result$padj)*result$log2FoldChange
  threshold<-rep(1,nrow(result))
  threshold[result$log2FoldChange>log2(1.5) & result$padj<0.1]=2
  threshold[result$log2FoldChange< -log2(1.5) & result$padj<0.1]=3
  result_temp<-result[order(-result$score),]
  gene2label<-na.omit(c(head(result_temp$gene.symbol),tail(result_temp$gene.symbol)))
  gene2label<-gene2label[gene2label %in% result$gene.symbol[result$padj<0.1]]
  points<-rep("",dim(result)[1])
  ind<-match(gene2label,result$gene.symbol)
  points[ind]<-gene2label
  p1=ggplot(data=result, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(colour=as.factor(threshold))) +
    scale_colour_manual(values=c("grey60","yellow2","dodgerblue4"),breaks=c(1,2,3),
                        labels=c("No sig. diff.",paste("Higher in Section",(i+1)),paste("Higher in Section 1")))+
    #geom_text(aes(label=points),size=2.8,color="black") +
    theme(legend.position="top",text = element_text(size=15),legend.title=element_blank(),
          legend.text=element_text(size=15, face = "bold"), 
          axis.text.x  = element_text(size=15),axis.text.y  = element_text(size=15),
          axis.title.x=element_text(size=20, face = "bold"),axis.title.y=element_text(size=20, face = "bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = "white")) +
    xlab("log2(fold change)") + ylab("-log10(adj.p-value)") + scale_size(guide="none") +
    geom_vline(xintercept = 1.5, linetype = "dashed") +
    geom_vline(xintercept = -1.5, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed")
  
  
  
  png(paste("Volcano_",lines_i,"_",lines_j[i],".png",sep=""),width=7,height=5, units = "in", res = 300)
  print(p1)
  dev.off()
  #Gene Ontology enrichment
  if (length(unique(subset(result,padj<0.1 & log2FoldChange>log2(1.5),select="gene.symbol",drop=T))) > 0){
    bp_up<-fisherTopGO(unique(subset(result,padj<0.1 & log2FoldChange>log2(1.5),select="gene.symbol",drop=T)),
                       unique(background),"BP","classic") #BP > MF or CC
    filename<-paste("top_GO_BP_up_in_",lines_i,"_vs_",lines_j[i],".xlsx",sep="")
    GOenrichTopTerms(bp_up,topK=500,max_gene=500,max_char=100,filename=filename)
    go<-read.table(filename,header=T,sep="\t",stringsAsFactors = F,comment.char = "",quote="")
    go$Term<-capitalize(tolower(go$Term))
    go$Term<-str_wrap(go$Term,width=25)
    go$Term<-factor(go$Term,levels=rev(go$Term),ordered=T)
    p2 = ggplot(go[1:5,], aes(x=Term, y=-log10(FDR),fill=-log10(FDR))) + 
      geom_bar(stat='identity') + coord_flip() +
      xlab("GO Terms") + ylab("-log10(FDR)") +
      scale_fill_gradient(low="grey90",high="yellow2",name="-log10(FDR)") + 
      theme(axis.text.x = element_text(size=15, face = "bold"),axis.text.y=element_text(size=23, face = "bold"),
            axis.title=element_text(size=20, face = "bold"),
            legend.position=c(1,0),legend.justification=c(1,0),legend.title=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      geom_hline(yintercept = 1, colour = "steelblue", size = 2, linetype = "dashed")
    png(paste("top_GO_BP_up_in_",lines_i,"_vs_",lines_j[i],".png",sep=""),height=12,width=8, units = "in", res = 300)
    print(p2)
    dev.off()}
  if (length(unique(subset(result,padj<0.1 & log2FoldChange< -log2(1.5),select="gene.symbol",drop=T))) > 0){
    bp_down<-fisherTopGO(unique(subset(result,padj<0.1 & log2FoldChange< -log2(1.5),select="gene.symbol",drop=T)),
                       unique(background),"BP","classic") #BP > MF or CC
    filename<-paste("top_GO_BP_down_in_",lines_i,"_vs_",lines_j[i],".xlsx",sep="")
    GOenrichTopTerms(bp_down,topK=500,max_gene=500,max_char=100,filename=filename)
    go.d<-read.table(filename,header=T,sep="\t",stringsAsFactors = F,comment.char = "",quote="")
    go.d$Term<-capitalize(tolower(go.d$Term))
    go.d$Term<-str_wrap(go.d$Term,width=25)
    go.d$Term<-factor(go.d$Term,levels=rev(go.d$Term),ordered=T)
    p3 = ggplot(go.d[1:5,], aes(x=Term, y=-log10(FDR),fill=-log10(FDR))) + 
      geom_bar(stat='identity') + 
      coord_flip() +xlab("GO Terms") + 
      ylab("-log10(FDR)") + 
      scale_fill_gradient(low="grey90",high="dodgerblue4",name="-log10(FDR)") + 
      theme(axis.text.x = element_text(size=15, face = "bold"),axis.text.y=element_text(size=23, face = "bold"),
            axis.title=element_text(size=20, face = "bold"),
            legend.position=c(1,0),legend.justification=c(1,0),legend.title=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      geom_hline(yintercept = 1, colour = "orangered", size = 2, linetype = "dashed")
    
    png(paste("top_GO_BP_down_in_",lines_i,"_vs_",lines_j[i],".png",sep=""),height=12,width=8, units = "in", res = 300)
    print(p3)
    dev.off()}
  else{
    next
  }
  
  
  
}


library(edgeR)
library(statmod)
y = DGEList(counts = mf,
            group = meta$meta)
treatment = factor(meta$meta)
design = model.matrix(~0 + treatment)
colnames(design) <- levels(treatment)
y = estimateDisp(y, design, robust = T)
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)

s5vs1 = makeContrasts(Section_5 - Section_1, levels=design)
s4vs1 = makeContrasts(Section_4 - Section_1, levels=design)
s3vs1 = makeContrasts(Section_3 - Section_1, levels=design)
s2vs1 = makeContrasts(Section_2 - Section_1, levels=design)

source("add_gene_info.R")
lines_c = list(s5vs1, s4vs1, s3vs1, s2vs1)
for (i in seq_along(lines_c)){
  #fit LRT
  res = glmLRT(fit, contrast = lines_c[[i]])
  topTags(res)
  res$table$FDR = p.adjust(res$table$PValue, method = "BH")
  res = res$table[res$table$FDR<0.1,]
  result<-add_gene_info(res,ensg2symbol)
  result = as.data.frame(result[complete.cases(result),])
  write.xlsx(result, paste("edgeR_all_genes_LRT",lines_j[i],"_vs_",lines_i,".xlsx",sep = ""))
  #fit QLF-test
  res <- glmQLFTest(fit, contrast= lines_c[[i]])
  topTags(res)
  res$table$FDR = p.adjust(res$table$PValue, method = "BH")
  res = res$table[res$table$FDR<0.1,]
  result1<-add_gene_info(res,ensg2symbol)
  result1 = as.data.frame(result[complete.cases(result),])
  write.xlsx(result1, paste("edgeR_all_genes_QLFT",lines_j[i],"_vs_",lines_i,".xlsx",sep = ""))
  
}
