genes = c("Gmps", "Nadsyn1", "Slc36a1")
ind = match(genes, rownames(mf))
all.gl = mf[ind,]

rn = rownames(all.gl)
all.gl = as.data.frame(c(all.gl[,c(1,6,11)], 
                         all.gl[,c(5,10,15)]))

rownames(all.gl) = rn

#make box plots for all primary genes
library(stringr)
library(R.utils)
all.gl = as.data.frame(t(all.gl))
colnames(all.gl) = capitalize(tolower(colnames(all.gl)))
all.gl$x = c(rep("Sec1",3), rep("Sec5",3))
all.gl$x = factor(all.gl$x, levels = c("Sec1","Sec5"), ordered = T)

class(all.gl)  

#how about making data frames for each genes and then plot?
df = list(data.frame())
for (i in 1:(ncol(all.gl)-1)) {
  df[[i]] = data.frame(all.gl[,i],all.gl[,ncol(all.gl)])
  colnames(df[[i]]) = c(colnames(all.gl)[i], "x")
}

pl = list()
for (i in 1:3) {
  min_axis<-ifelse(min(df[[i]][,1])<5,0,round(min(df[[i]][,1])-5,0))
  max_axis<-round(ifelse((max(df[[i]][ ,1] )>1000),max(df[[i]][ ,1] )+40,
                         ifelse((max(df[[i]][ ,1] )>500&max(df[[i]][ ,1] )<1000),
                                max(df[[i]][ ,1] )+30,
                                ifelse((max(df[[i]][ ,1] )<500&max(df[[i]][ ,1] )>100),
                                       max(df[[i]][ ,1] )+20,
                                       ifelse((max(df[[i]][ ,1] )<100&max(df[[i]][ ,1] )>50),
                                              max(df[[i]][ ,1] )+10,
                                              ifelse((max(df[[i]][ ,1] )<50&max(df[[i]][ ,1] )>10),
                                                     max(df[[i]][ ,1] )+7,
                                                     ifelse((max(df[[i]][ ,1] )<10&max(df[[i]][ ,1] )>5),
                                                            max(df[[i]][ ,1] )+2,
                                                            ifelse((max(df[[i]][ ,1] )<5&max(df[[i]][ ,1] )>3),
                                                                   max(df[[i]][ ,1] )+1,
                                                                   ifelse((max(df[[i]][ ,1] )<0.1&max(df[[i]][ ,1] )>0.01),
                                                                          max(df[[i]][ ,1] )+0.05,
                                                                          max(df[[i]][ ,1] )+0.5)))))))),1)
  ##match each gene and use if statement to put asterisks if significant
  p1 <- ggplot(df[[i]], aes(x = factor(x), y = df[[i]][,1]), fill = df[[i]]$x) +
    geom_boxplot(alpha = 0.80, color = c("dodgerblue4", "yellow2"), size = 3) +
    geom_point(size = 3, shape = 21,alpha = 0) +
    geom_jitter(size = 2, alpha = 0.3, width = 0.2)+
    xlab("") + ylab("TPM") + 
    theme(axis.text.y=element_text(size=20), 
          legend.text = element_text(size = 0),
          axis.text.x=element_text(size=30, face = "bold"),
          axis.title.y=element_text(size=30, face = "bold"),
          legend.position=c(10,1),legend.justification=c(1,1),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 40, face = "bold"))+
    labs(title = colnames(df[[i]])[1])+
    ylim(min_axis,max_axis)
  
  
  tiff(paste("Sec1vs5_boxplot_for_",colnames(df[[i]])[1],".tiff"), 
       width = 6, height = 8, units = "in", res =300)
  print(p1)
  dev.off()
}
