#alternative gene in different datasets
setwd("~/Desktop/practice/Mouse_intestines_bulk/alternative-genes-other-data/")
library(openxlsx)
library(ggplot2)
library(DESeq2)
library(dplyr)
######1) E-MTAB-7765 -- #Nt5e, Cd38, Slc43a2

#read data
df = read.xlsx("Counts_table_FLox-Lck.xlsx")
rn = df$`#gene`
rownames(df) = rn
df = df[,-c(colnames(df) == "#gene")]
colnames(df)

colnames(df) = gsub(pattern = "-[0-9]*", replacement = "", colnames(df))
colnames(df)

countData = as.data.frame(c(df[,colnames(df) == "FLox_FLox.counts"],
                            df[,colnames(df) == "Lck_Lck.counts"]))
countData = data.frame(apply(countData, 2, function(x) as.numeric(as.character(x))))
colnames(countData) = gsub(pattern = "[0-9]", replacement = "", colnames(countData))
colnames(countData) = gsub(pattern = "\\.", replacement = "", colnames(countData))

rownames(countData) = rownames(df)
colData = data.frame(colData = c(colnames(df)[colnames(df) == "FLox_FLox.counts"],
                                 colnames(df)[colnames(df) == "Lck_Lck.counts"]))
dds<-DESeqDataSetFromMatrix(countData = countData,
                            colData = colData,
                            design = ~0+colData)
dds$colData = factor(colData$colData, levels = c("FLox_FLox.counts", "Lck_Lck.counts"))
dds = DESeq(dds)
res = results(dds)
res = res[complete.cases(res),]
res = res[res$padj<0.1,]
result = as.data.frame(res[complete.cases(res),])
#result<-add_gene_info(result,df)
write.xlsx(result, paste("DESeq2_E-MTAB-7765.xlsx",sep = ""))


##E-MTAB-2350 -- Cuffdiff txt file -- Nt5e, Slc43a2

##E-GEOD-73832 -- Nt5e, Slc43a2; no raw expr data
#read metadata
meta = read.xlsx("small_intestinal_neuroendocrine_tumor_metadata.xlsx")
meta[1:10,1:10]
meta = meta[1:133,]
#select only small intestine
meta = meta[meta$`Characteristics.[sample.site]` == "small intestine",]
meta = meta[complete.cases(meta$Source.Name),]
sum(meta$`Comment.[Sample_source_name]` == "normal")
sum(meta$`Comment.[Sample_source_name]` == "primary")

##GSE9576 (small intestine neuroendocrine tumor = midgut carcinoid tumor) -- NT5E, CD38

