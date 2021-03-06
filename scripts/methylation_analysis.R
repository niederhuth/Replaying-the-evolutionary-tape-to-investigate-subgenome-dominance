#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(scales)

samples <- read.csv('../misc/samples.csv',header=T)
plots <- c('all_genes','R500_genes','TO1000_genes','R500_syntelogs','R500_non-syntelogs','TO1000_syntelogs',
'TO1000_non-syntelogs','all_LTRs','R500_LTRs','TO1000_LTRs','R500_nonbiased','R500_BnAbiased','R500_BnCbiased',
'TO1000_nonbiased','TO1000_BnAbiased','TO1000_BnCbiased')
controls <- c("TO1000","IMB218","mock")

#Metaplots
setwd("metaplots/")
#All samples CG & CHG metaplots
for(x in plots){
  df1 <- matrix(nrow=0,ncol=0)
  for(i in samples$Sample){
    input=paste(x,paste(i,x,"metaplot.txt.gz",sep="_"),sep="/")
    if(file.exists(input)){
      df2 <- read.table(input,header=T,sep="\t")[2:4]
      colnames(df2) <- c('CG','CHG','CHH')
      df2$sample <- c(i)
      df2 <- melt(df2)
      df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
      df1 <- rbind(df1,df2)
      rm(df2)
    } else {
      next
    }
  }
  df1$window <- c(1:60)

  for(y in c('CG','CHG')){
    df2 <- df1[df1$variable==y,]
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

    ggsave(paste(x,y,"metaplot.pdf",sep="_"),plot,path="plots/")
  }

  for(y in c('CHH')){
    df2 <- df1[df1$variable==y,]
    plot <- ggplot(df2, aes(x=window,y=value,group=sample,color=generation)) +
      geom_line(size=0.8) +
      theme(panel.background=element_blank(),panel.grid=element_blank(),axis.text.y=element_text(color="black"),
            axis.text.x=element_text(color="black"),axis.ticks=element_line(color="black"),
            axis.title=element_text(color="black"),legend.position="right", axis.line=element_line(color="black")) +
      ylab("Percent methylation") +
      xlab("") +
      scale_y_continuous(limits=c(0,0.125),expand=c(0,0),breaks=c(0.025,0.05,0.075,0.1,0.125),labels=c("2.5%","5.0%","7.5%","10.0%","12.5%")) +
      geom_vline(xintercept=20,linetype="longdash",color="grey55") +
      geom_vline(xintercept=40,linetype="longdash",color="grey55") +
      scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1, 20, 40, 60))

    ggsave(paste(x,y,"metaplot.pdf",sep="_"),plot,path="plots/")
  }
}

#Line CG & CHG metaplots
for(x in plots){
  for(j in line <- unique(samples$Line)){
    df1 <- matrix(nrow=0,ncol=0)
    for(i in c(as.vector(samples[samples$Line==j,]$Sample),controls)){
      input=paste(x,paste(i,x,"metaplot.txt.gz",sep="_"),sep="/")
      if(file.exists(input)){
        df2 <- read.table(input,header=T,sep="\t")[2:4]
        colnames(df2) <- c('CG','CHG','CHH')
        df2$sample <- c(i)
        df2 <- melt(df2)
        df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
        df1 <- rbind(df1,df2)
        rm(df2)
      } else {
        next
      }
    }
  }
    df1$window <- c(1:60)

  for(y in c('CG','CHG')){
    df2 <- df1[df1$variable==y,]
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

    ggsave(paste(x,j,y,"metaplot.pdf",sep="_"),plot,path="plots/")
  }

  for(y in c('CHH')){
    df2 <- df1[df1$variable==y,]
    plot <- ggplot(df2, aes(x=window,y=value,group=sample,color=generation)) +
      geom_line(size=0.8) +
      theme(panel.background=element_blank(),panel.grid=element_blank(),axis.text.y=element_text(color="black"),
            axis.text.x=element_text(color="black"),axis.ticks=element_line(color="black"),
            axis.title=element_text(color="black"),legend.position="right", axis.line=element_line(color="black")) +
      ylab("Percent methylation") +
      xlab("") +
      scale_y_continuous(limits=c(0,0.125),expand=c(0,0),breaks=c(0.025,0.05,0.075,0.1,0.125),labels=c("2.5%","5.0%","7.5%","10.0%","12.5%")) +
      geom_vline(xintercept=20,linetype="longdash",color="grey55") +
      geom_vline(xintercept=40,linetype="longdash",color="grey55") +
      scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1, 20, 40, 60))
      scale_color_manual("",values=c("tomato3","orange3","magenta3","royalblue3"))

    ggsave(paste(x,j,y,"metaplot.pdf",sep="_"),plot,path="plots/")
  }
}

#Bias plots
for(a in c("R500","TO1000")){
  df1 <- matrix(nrow=0,ncol=0)
  for(b in c("BnA","BnC","Non")){
    for(i in samples$Sample){
      input=paste(paste(a,"_",b,"biased",sep=""),paste(i,a,paste(b,"biased",sep=""),"metaplot.txt.gz",sep="_"),sep="/")
      if(file.exists(input)){
        df2 <- read.table(input,header=T,sep="\t")[2:4]
        colnames(df2) <- c('CG','CHG','CHH')
        df2$sample <- c(i)
        df2$sample_bias <- c(paste(i,b,sep="_"))
        df2 <- melt(df2)
        df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
        df2$bias <- c(b)
        df1 <- rbind(df1,df2)
        rm(df2)
      } else {
        next
      }
    }
  }
  df1$window <- c(1:60)

  for(d in c('CG','CHG')){
    df2 <- df1[df1$variable==d,]
    plot <- ggplot(df2, aes(x=window,y=value,group=sample_bias,color=bias)) +
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
    #scale_color_manual("",values=c("tomato3","orange3","magenta3","royalblue3"))

    ggsave(paste(a,'bias',d,"metaplot.pdf",sep="_"),plot,path="plots/")
  }

  for(d in c('CHH')){
    df2 <- df1[df1$variable==d,]
    plot <- ggplot(df2, aes(x=window,y=value,group=sample_bias,color=bias)) +
      geom_line(size=0.8) +
      theme(panel.background=element_blank(),panel.grid=element_blank(),axis.text.y=element_text(color="black"),
            axis.text.x=element_text(color="black"),axis.ticks=element_line(color="black"),
            axis.title=element_text(color="black"),legend.position="right", axis.line=element_line(color="black")) +
      ylab("Percent methylation") +
      xlab("") +
      scale_y_continuous(limits=c(0,0.125),expand=c(0,0),breaks=c(0.025,0.05,0.075,0.1,0.125),labels=c("2.5%","5.0%","7.5%","10.0%","12.5%")) +
      geom_vline(xintercept=20,linetype="longdash",color="grey55") +
      geom_vline(xintercept=40,linetype="longdash",color="grey55") +
      scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1, 20, 40, 60))

    ggsave(paste(a,'bias',d,"metaplot.pdf",sep="_"),plot,path="plots/")
  }
}

#Total methylation levels
setwd("../total_methylation")

for(x in c('combined','R500','TO1000')){
  df1 <- matrix(nrow=0,ncol=0)
  for(i in samples$Sample){
    input=paste(x,paste(i,x,"weighted_methylation.txt.gz",sep="_"),sep="/")
    df2 <- t(read.table(input,header=T,sep="\t")[6])
    colnames(df2) <- c('CG','CHG','CHH')
    df2 <- melt(df2)
    df2$sample <- c(i)
    df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
    df2$generation2 <- c(as.numeric(gsub("G","",samples[samples$Sample==i,]$Generation)))
    df2$order <- c(as.numeric(samples[samples$Sample==i,]$Order))
    df1 <- rbind(df1,df2)
    rm(df2)
  }
  df1$generation2 <- ifelse(is.na(df1$generation2),0,df1$generation2)

  plot <- ggplot(df1,aes(x=reorder(sample,order),y=value,fill=Var2)) +
    geom_bar(stat="identity",position="dodge") +
    theme(panel.background=element_blank(),
          panel.grid=element_blank(),
          axis.text.y=element_text(color="black"),
          axis.text.x=element_text(color="black",angle=315,hjust=0),
          axis.ticks=element_line(color="black"),
          axis.title=element_text(color="black"),
          legend.position="right",
          axis.line=element_line(color="black")) +
    ylab("Percent methylation") +
    xlab("") +
    scale_y_continuous(expand=c(0,0), labels=percent)
    #scale_fill_manual("",values=c("dodgerblue2","darkolivegreen3","tomato2"))
  ggsave(paste(x,"total_methylation.pdf",sep="_"),plot,path="plots/")

  for(y in c('CG','CHG','CHH')){
    df2 <- df1[df1$Var2==y,]
    plot <- ggplot(df2,aes(x=reorder(sample,generation2),y=value,fill=generation)) +
              geom_bar(stat="identity",position="dodge") +
              theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.text.y=element_text(color="black"),
                    axis.text.x=element_text(color="black",angle=315,hjust=0),
                    axis.ticks=element_line(color="black"),
                    axis.title=element_text(color="black"),
                    legend.position="right",
                    axis.line=element_line(color="black")) +
              ylab("Percent methylation") +
              xlab("") +
              scale_y_continuous(expand=c(0,0), labels=percent)
              #scale_fill_manual("",values=c("tomato3","orange3","magenta3","royalblue3"))
    ggsave(paste(x,y,"_total_methylation.pdf",sep="_"),plot,path="plots/")
  }
}

#Line total methylation
for(j in line <- unique(samples$Line)){
  for(x in c('combined','R500','TO1000')){
    df1 <- matrix(nrow=0,ncol=0)
    for(i in c(as.vector(samples[samples$Line==j,]$Sample),controls)){
      input=paste(x,paste(i,x,"weighted_methylation.txt.gz",sep="_"),sep="/")
      df2 <- t(read.table(input,header=T,sep="\t")[6])
      colnames(df2) <- c('CG','CHG','CHH')
      df2 <- melt(df2)
      df2$sample <- c(i)
      df2$generation <- c(as.character(samples[samples$Sample==i,]$Generation))
      df2$generation2 <- c(as.numeric(gsub("G","",samples[samples$Sample==i,]$Generation)))
      df2$order <- c(as.numeric(samples[samples$Sample==i,]$Order))
      df1 <- rbind(df1,df2)
      rm(df2)
    }
    df1$generation2 <- ifelse(is.na(df1$generation2),0,df1$generation2)
    for(y in c('CG','CHG','CHH')){
      df2 <- df1[df1$Var2==y,]
      plot <- ggplot(df2,aes(x=reorder(sample,generation2),y=value,fill=generation)) +
                geom_bar(stat="identity",position="dodge") +
                theme(panel.background=element_blank(),
                      panel.grid=element_blank(),
                      axis.text.y=element_text(color="black"),
                      axis.text.x=element_text(color="black",angle=315,hjust=0),
                      axis.ticks=element_line(color="black"),
                      axis.title=element_text(color="black"),
                      legend.position="right",
                      axis.line=element_line(color="black")) +
                ylab("Percent methylation") +
                xlab("") +
                scale_y_continuous(expand=c(0,0), labels=percent)
                #scale_fill_manual("",values=c("tomato3","orange3","magenta3","royalblue3"))
      ggsave(paste(j,x,y,"_total_methylation.pdf",sep="_"),plot,path="plots/")
    }
  }
}

#DMR analysis
setwd("../DMRs")
DMRs=c('CNN','CGN','CHG','CHH')
for(i in DMRs){
  input=paste(i,"_filtered_rms_results_collapsed.tsv.gz",sep="")
  df1 <- read.table(input,header=T,sep="\t")
  df2 <- as.matrix(log2(df1[7:length(df1)]*100+0.000001))
  df2 <- na.omit(df2)
  colnames(df2) <- gsub("methylation_level_","",colnames(df1[,7:length(df1)]))
  df2 <- df2[,order(c(1,7,13,2,8,14,3,9,15,4,10,16,5,11,17,6,12,18))]
  #Heatmap(df2,cluster_columns=F,heatmap_legend_param=list(title="Log2(%mCHH)"))

  df3 <- df1[7:length(df1)]
  colnames(df3) <- gsub("methylation_level_","",colnames(df1[,7:length(df1)]))
  df3 <- df3[,order(c(1,7,13,2,8,14,3,9,15,4,10,16,5,11,17,6,12,18))]
  df4 <- melt(df3)
  df4$generation <- gsub(".*S","S",df4$variable)
  x <- c()
  for(y in 1:18){x <- c(x,rep(y,nrow(df3)))}
  df4$order <- x
  rm(x)

  plot <- ggplot(df4, aes(x=variable,y=value,fill=generation)) +
    geom_boxplot(notch=T) +
    theme(panel.background=element_blank(),panel.grid=element_blank(),axis.text.y=element_text(color="black"),
          axis.text.x=element_text(color="black",angle=90,hjust=0,vjust=0),axis.ticks=element_line(color="black"),
          axis.title=element_text(color="black"),legend.position="right", axis.line=element_line(color="black")) +
    xlab("Sample") +
    ylab("Methylation level") +
    scale_y_continuous(labels=percent,expand=c(0,0))

  ggsave(paste(i,"by_line_boxplot.pdf",sep="_"),plot,path="plots/")

  plot <- ggplot(df4, aes(reorder(x=variable,order),y=value,fill=generation)) +
    geom_boxplot(notch=T) +
    theme(panel.background=element_blank(),panel.grid=element_blank(),axis.text.y=element_text(color="black"),
          axis.text.x=element_text(color="black",angle=90,hjust=0,vjust=0),axis.ticks=element_line(color="black"),
          axis.title=element_text(color="black"),legend.position="right", axis.line=element_line(color="black")) +
    xlab("Sample") +
    ylab("Methylation level") +
    scale_y_continuous(labels=percent,expand=c(0,0))

  ggsave(paste(i,"by_generation_boxplot.pdf",sep="_"),plot,path="plots/")
}

setwd("../")
