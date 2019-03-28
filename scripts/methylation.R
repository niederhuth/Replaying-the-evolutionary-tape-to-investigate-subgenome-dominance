#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

samples <- read.csv('../misc/samples.csv',header=T)

#Metaplots
setwd("metaplots/")
#All samples CG & CHG metaplots
for(x in c('all_genes','R500_genes','R500_syntelogs','R500_non-syntelogs','TO1000_genes','TO1000_syntelogs','TO1000_non-syntelogs')){
  df <- matrix(nrow=0,ncol=0)
  for(i in names <- c('100S5','100S10','1100S1','1100S5','1100S10','200S1','200S5','200S10','300S5','300S10','400S1','400S5','400S10','600S1','600S5','600S10','TO1000')){
    input=paste(x,paste(i,"uncorrected",x,"metaplot.txt",sep="_"),sep="/")
    df2 <- read.table(input,header=T,sep="\t")[2:4]
    colnames(df2) <- c('CG','CHG','CHH')
    df2$sample <- c(i)
    df2 <- melt(df2)
    df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
    df <- rbind(df,df2)
    rm(df2)
  }
  df$window <- c(1:60)
  
  for(y in c('CG','CHG')){
    df2 <- df[df$variable==y,]
    plot <- ggplot(df2, aes(x=window,y=value,group=sample,color=generation)) + 
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
      scale_color_manual("",values=c("tomato3","orange3","magenta3","royalblue3"))
    
    ggsave(paste(x,y,"_metaplot.pdf",sep="_"),plot,path="plots/")
  }
}  

#All samples CHH metaplots
for(x in c('all_genes','R500_genes','R500_syntelogs','R500_non-syntelogs','TO1000_genes','TO1000_syntelogs','TO1000_non-syntelogs')){
  df <- matrix(nrow=0,ncol=0)
  for(i in names <- c('100S5','100S10','1100S1','1100S5','1100S10','200S1','200S5','200S10','300S5','300S10','400S1','400S5','400S10','600S1','600S5','600S10','TO1000')){
    input=paste(x,paste(i,"uncorrected",x,"metaplot.txt",sep="_"),sep="/")
    df2 <- read.table(input,header=T,sep="\t")[2:4]
    colnames(df2) <- c('CG','CHG','CHH')
    df2$sample <- c(i)
    df2 <- melt(df2)
    df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
    df <- rbind(df,df2)
    rm(df2)
  }
  df$window <- c(1:60)
  
  for(y in c('CHH')){
    df2 <- df[df$variable==y,]
    plot <- ggplot(df2, aes(x=window,y=value,group=sample,color=generation)) + 
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

#Line CG & CHG metaplots
for(x in c('all_genes','R500_genes','R500_syntelogs','R500_non-syntelogs','TO1000_genes','TO1000_syntelogs','TO1000_non-syntelogs')){
  for(j in line <- unique(samples$Line)){
    df <- matrix(nrow=0,ncol=0)
    for(i in names <- samples[samples$Line==j,]$Sample ){
      input=paste(x,paste(i,"uncorrected",x,"metaplot.txt",sep="_"),sep="/")
      df2 <- read.table(input,header=T,sep="\t")[2:4]
      colnames(df2) <- c('CG','CHG','CHH')
      df2$sample <- c(i)
      df2 <- melt(df2)
      df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
      df <- rbind(df,df2)
      rm(df2)
    }
    df$window <- c(1:60)
  
    for(y in c('CG','CHG')){
      df2 <- df[df$variable==y,]
      plot <- ggplot(df2, aes(x=window,y=value,group=sample,color=generation)) + 
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
        scale_color_manual("",values=c("tomato3","orange3","magenta3","royalblue3"))
    
      ggsave(paste(x,j,y,"_metaplot.pdf",sep="_"),plot,path="plots/")
    }
  }
}  

#Line CHH metaplots
for(x in c('all_genes','R500_genes','R500_syntelogs','R500_non-syntelogs','TO1000_genes','TO1000_syntelogs','TO1000_non-syntelogs')){
  for(j in line <- unique(samples$Line)){
    df <- matrix(nrow=0,ncol=0)
    for(i in names <- samples[samples$Line==j,]$Sample ){
      input=paste(x,paste(i,"uncorrected",x,"metaplot.txt",sep="_"),sep="/")
      df2 <- read.table(input,header=T,sep="\t")[2:4]
      colnames(df2) <- c('CG','CHG','CHH')
      df2$sample <- c(i)
      df2 <- melt(df2)
      df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
      df <- rbind(df,df2)
      rm(df2)
    }
    df$window <- c(1:60)
    
    for(y in c('CHH')){
      df2 <- df[df$variable==y,]
      plot <- ggplot(df2, aes(x=window,y=value,group=sample,color=generation)) + 
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
        scale_color_manual("",values=c("tomato3","orange3","magenta3","royalblue3"))
      
      ggsave(paste(x,j,y,"_metaplot.pdf",sep="_"),plot,path="plots/")
    }
  }
}  

#Total methylation levels
setwd("../total_methylation")
