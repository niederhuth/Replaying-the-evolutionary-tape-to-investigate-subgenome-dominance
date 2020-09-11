#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(scales)

samples <- read.csv('../misc/samples.csv',header=T)
setwd("metaplots/")

controls <- c("TO1000","R500","mock")
colors=c("TO1000"="blue","R500"="red","mock"="black","G1"="orchid","G5"="orchid","G10"="orchid")

plots <- c('all_genes','R500_genes','TO1000_genes','R500_syntelogs','R500_non-syntelogs','TO1000_syntelogs',
'TO1000_non-syntelogs','R500_nonbiased','R500_BnAbiased','R500_BnCbiased','TO1000_nonbiased','TO1000_BnAbiased',
'TO1000_BnCbiased')
df1 <- matrix(nrow=0,ncol=0)
for(x in plots){
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

	for(y in c("G1","G5","G10")){
		sample_list <- c(y,controls)
		plot <- ggplot(df1[df1$variable=="CG" & df1$generation %in% sample_list,]) + 
			geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
			theme(panel.background=element_blank(),panel.grid=element_blank(),
				axis.text=element_text(color="black",size=12),
				axis.ticks=element_line(color="black"),
				axis.title=element_text(color="black",size=12),
				legend.position="right",axis.line=element_line(color="black"),
				text=element_text(size=12)) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
				labels=percent) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
			scale_color_manual("",values=colors)
		ggsave(paste(x,y,"CG_metaplot.pdf",sep="_"),plot,path="new_plots/")
		plot <- plot + theme(legend.position="None")
		ggsave(paste(x,y,"CG_metaplot.pdf",sep="_"),plot,path="new_plots_no_legend/")

		plot <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
			geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
			theme(panel.background=element_blank(),panel.grid=element_blank(),
				axis.text=element_text(color="black",size=12),
				axis.ticks=element_line(color="black"),
				axis.title=element_text(color="black",size=12),
				legend.position="right",axis.line=element_line(color="black"),
				text=element_text(size=12)) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,0.3),expand=c(0,0),breaks=c(0.1,0.2,0.3),
				labels=percent) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
			scale_color_manual("",values=colors)
		ggsave(paste(x,y,"CHG_metaplot.pdf",sep="_"),plot,path="new_plots/")
		plot <- plot + theme(legend.position="None")
		ggsave(paste(x,y,"CHG_metaplot.pdf",sep="_"),plot,path="new_plots_no_legend/")

		plot <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
			geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
			theme(panel.background=element_blank(),panel.grid=element_blank(),
				axis.text=element_text(color="black",size=12),
				axis.ticks=element_line(color="black"),
				axis.title=element_text(color="black",size=12),
				legend.position="right",axis.line=element_line(color="black"),
				text=element_text(size=12)) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,0.075),expand=c(0,0),breaks=c(0.025,0.05,0.075),
				labels=percent) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
			scale_color_manual("",values=colors)
		ggsave(paste(x,y,"CHH_metaplot.pdf",sep="_"),plot,path="new_plots/")
		plot <- plot + theme(legend.position="None")
		ggsave(paste(x,y,"CHH_metaplot.pdf",sep="_"),plot,path="new_plots_no_legend/")
	}
	df1 <- matrix(nrow=0,ncol=0)
}

plots <- c('all_LTRs','R500_LTRs','TO1000_LTRs')

df1 <- matrix(nrow=0,ncol=0)
for(x in plots){
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

	for(y in c("G1","G5","G10")){
		sample_list <- c(y,controls)
		plot <- ggplot(df1[df1$variable=="CG" & df1$generation %in% sample_list,]) + 
			geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
			theme(panel.background=element_blank(),panel.grid=element_blank(),
				axis.text=element_text(color="black",size=12),
				axis.ticks=element_line(color="black"),
				axis.title=element_text(color="black",size=12),
				legend.position="right",axis.line=element_line(color="black"),
				text=element_text(size=12)) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0.25,0.5,0.75,1),
				labels=percent) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
			scale_color_manual("",values=colors)
		ggsave(paste(x,y,"CG_metaplot.pdf",sep="_"),plot,path="new_plots/")
		plot <- plot + theme(legend.position="None")
		ggsave(paste(x,y,"CG_metaplot.pdf",sep="_"),plot,path="new_plots_no_legend/")

		plot <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
			geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
			theme(panel.background=element_blank(),panel.grid=element_blank(),
				axis.text=element_text(color="black",size=12),
				axis.ticks=element_line(color="black"),
				axis.title=element_text(color="black",size=12),
				legend.position="right",axis.line=element_line(color="black"),
				text=element_text(size=12)) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
				labels=percent) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
			scale_color_manual("",values=colors)
		ggsave(paste(x,y,"CHG_metaplot.pdf",sep="_"),plot,path="new_plots/")
		plot <- plot + theme(legend.position="None")
		ggsave(paste(x,y,"CHG_metaplot.pdf",sep="_"),plot,path="new_plots_no_legend/")

		plot <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
			geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
			theme(panel.background=element_blank(),panel.grid=element_blank(),
				axis.text=element_text(color="black",size=12),
				axis.ticks=element_line(color="black"),
				axis.title=element_text(color="black",size=12),
				legend.position="right",axis.line=element_line(color="black"),
				text=element_text(size=12)) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,0.125),expand=c(0,0),breaks=c(0.025,0.075,0.125),
				labels=percent) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
			scale_color_manual("",values=colors)
		ggsave(paste(x,y,"CHH_metaplot.pdf",sep="_"),plot,path="new_plots/")
		plot <- plot + theme(legend.position="None")
		ggsave(paste(x,y,"CHH_metaplot.pdf",sep="_"),plot,path="new_plots_no_legend/")
	}
	df1 <- matrix(nrow=0,ncol=0)
}
