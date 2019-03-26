library(ggplot2)
library(reshape2)

#Read in data
data=read.table("all_genes/600S5_uncorrected_all_genes_metaplot.txt",header=T,sep="\t",row.names=1)
data2=melt(data)
data2$window=c(1:60)

#Plot data
ggplot(data2,aes(x=window,y=value,color=variable))+geom_line(size=0.8) +
  theme(panel.background=element_blank(), panel.grid=element_blank(),
        axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
        axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
        legend.position="right", axis.line=element_line(color="black")) +
  ylab("Weighted methylation") + 
  xlab( "" ) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
  geom_vline(xintercept=20, linetype="longdash", color="grey55") +
  geom_vline(xintercept=40, linetype="longdash", color="grey55") +
  scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60)) #+
  #scale_color_manual("",values=c("","",""))
