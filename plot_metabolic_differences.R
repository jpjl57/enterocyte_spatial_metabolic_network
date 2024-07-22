#this script plots stacked bar plots and volcano plots comparing section 1 and 5 
#villus dataset

library(dplyr)
library(R.utils)
library(stringr)
library(openxlsx)
library(ggplot2)

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
ftf.l = list()
hist.l = list()

load("phaseI.RData")
load("phaseII.RData")
load("oxi_stress.RData")
load("epi_modification.RData")
load("inflamm.RData")
load("nuclear_receptors.RData")
load("bile_acid.RData")
load("transporters.RData")
load("TF_drug_metabolism.RData")
load("mouse_TF_full.RData")
ftf = ftf$Symbol
#load("histone_modification.RData")

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
  #res = res[res$padj<0.1,]
  result = as.data.frame(res[complete.cases(res),])
  result<-add_gene_info(result,df)
  
  #match with genes of interest
  sig = result
  inflamm.l[[i]] = return_match(inflamm_gene_ls$X1, sig)
  trp.l[[i]] = return_match(tr.unique, sig)
  oxi.l[[i]] = return_match(oxi, sig)
  epi.l[[i]] = return_match(epi, sig)
  ph1.l[[i]] = return_match(ph1, sig)
  ph2.l[[i]] = return_match(ph2, sig)
  nr.l[[i]] = return_match(nr, sig)
  ba.l[[i]] = return_match(ba, sig)
  tf.l[[i]] = return_match(tf, sig)
  ftf.l[[i]] = return_match(ftf, sig)
  #hist.l[[i]] = return_match(hist, sig)
  }


#plot only section 1 vs 5

result = list(inflamm.l[[4]], trp.l[[4]], oxi.l[[4]], epi.l[[4]], ph1.l[[4]], ph2.l[[4]], nr.l[[4]], ba.l[[4]], tf.l[[4]], ftf.l[[4]])#, hist.l[[4]])
name = c("inflammation", "transporters", "oxidative_stress", "epigenetic_modification", "phase1", "phase2", "nuclear_receptor", "bile_acid", "transcription_factor",
         "transcription_factor_full")#, "histone_modification")
for (i in seq_along(result)) {
  result[[i]]$score<- -log10(result[[i]]$padj)*result[[i]]$log2FoldChange
  threshold<-rep(1,nrow(result[[i]]))
  threshold[result[[i]]$log2FoldChange>log2(1.5) & result[[i]]$padj<0.1]=2
  threshold[result[[i]]$log2FoldChange< -log2(1.5) & result[[i]]$padj<0.1]=3
  temp<-result[[i]][order(-result[[i]]$score),]
  gene2label<-na.omit(c(head(temp$gene.symbol),tail(temp$gene.symbol)))
  gene2label<-gene2label[gene2label %in% result[[i]]$gene.symbol[result[[i]]$padj<0.1]]
  points<-rep("",dim(result[[i]])[1])
  ind<-match(gene2label,result[[i]]$gene.symbol)
  points[ind]<-gene2label
  p1=ggplot(data=result[[i]], aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(colour=as.factor(threshold))) +
    scale_colour_manual(values=c("grey60","yellow2","dodgerblue4"),breaks=c(1,2,3),
                        labels=c("No sig. diff.",paste("Higher in Section 5"),paste("Higher in Section 1")))+
    geom_text(aes(label=points),size=2.8,color="black") +
    theme(legend.position="top",text = element_text(size=10),legend.title=element_blank(),
          legend.text=element_text(size=10, face = "bold"), axis.text.x  = element_text(size=10),axis.text.y  = element_text(size=10),
          axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = "white")) +
    xlab("log2 fold change") + ylab("-log10 adj.p-value") + scale_size(guide="none") +
    geom_vline(xintercept = 1.5, linetype = "dashed") +
    geom_vline(xintercept = -1.5, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed")
  
  
  png(paste("Volcano_",name[i],"_section1_vs_5.png",sep=""),width=5.5,height=6, units = "in", res = 300)
  print(p1)
  dev.off()
  }
  




#another way is to plot results in stacked bar plots
#up = fold change >0 and down = fold change <0

