#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

samples <- read.csv('../../../misc/samples.csv',header=T)

#CG & CHG metaplots
for(x in c('all_genes','R500_genes','R500_syntelogs','R500_non-syntelogs','TO1000_genes','TO1000_syntelogs','TO1000_non-syntelogs')){
  df <- matrix(nrow=0,ncol=0)
  for(i in names <- c('100S10','100S5','1100S10','1100S5','200S10','200S1','300S10','300S5','400S10','400S5','600S10','600S1','600S5','TO1000')){
    input=paste(x,paste(i,"uncorrected",x,"metaplot.txt",sep="_"),sep="/")
    df2 <- read.table(input,header=T,sep="\t")[2:4]
    colnames(df2) <- c('CG','CHG','CHH')
    df2$sample <- c(i)
    df <- rbind(df,melt(df2))
    rm(df2)
  }
  df$window <- c(1:60)
  
  for(y in c('CG','CHG')){
    df2 <- df[df$variable==y,]
    plot <- ggplot(df2, aes(x=window,y=value,color=sample)) + 
      geom_line(size=0.8) +
      theme(panel.background=element_blank(),panel.grid=element_blank(),axis.text.y=element_text(color="black"),
            axis.text.x=element_text(color="black"),axis.ticks=element_line(color="black"),
            axis.title=element_text(color="black"),legend.position="right", axis.line=element_line(color="black")) +
      ylab("Percent methylation") + 
      xlab("") +
      scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
      geom_vline(xintercept=20,linetype="longdash",color="grey55") +
      geom_vline(xintercept=40,linetype="longdash",color="grey55") +
      scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1, 20, 40, 60))
    
    ggsave(paste(x,y,"_metaplot.pdf",sep="_"),plot,path="plots/")
  }
}  

#CHH metaplots
for(x in c('all_genes','R500_genes','R500_syntelogs','R500_non-syntelogs','TO1000_genes','TO1000_syntelogs','TO1000_non-syntelogs')){
  df <- matrix(nrow=0,ncol=0)
  for(i in names <- c('100S10','100S5','1100S10','1100S5','200S10','200S1','300S10','300S5','400S10','400S5','600S10','600S1','600S5','TO1000')){
    input=paste(x,paste(i,"uncorrected",x,"metaplot.txt",sep="_"),sep="/")
    df2 <- read.table(input,header=T,sep="\t")[2:4]
    colnames(df2) <- c('CG','CHG','CHH')
    df2$sample <- c(i)
    df <- rbind(df,melt(df2))
    rm(df2)
  }
  df$window <- c(1:60)
  
  for(y in c('CHH')){
    df2 <- df[df$variable==y,]
    plot <- ggplot(df2, aes(x=window,y=value,color=sample)) + 
      geom_line(size=0.8) +
      theme(panel.background=element_blank(),panel.grid=element_blank(),axis.text.y=element_text(color="black"),
            axis.text.x=element_text(color="black"),axis.ticks=element_line(color="black"),
            axis.title=element_text(color="black"),legend.position="right", axis.line=element_line(color="black")) +
      ylab("Percent methylation") + 
      xlab("") +
      scale_y_continuous(limits=c(0,0.1),expand=c(0,0),breaks=c(0.025,0.05,0.075,0.1),labels=c("2.5%","5.0%","7.5%","10.0%")) +
      geom_vline(xintercept=20,linetype="longdash",color="grey55") +
      geom_vline(xintercept=40,linetype="longdash",color="grey55") +
      scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1, 20, 40, 60))
    
    ggsave(paste(x,y,"_metaplot.pdf",sep="_"),plot,path="plots/")
  }
} 
