#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(scales)

samples <- read.csv('../misc/samples.csv',header=T)
plots <- c('all_genes')
controls <- c("TO1000","R500","mock")

setwd("metaplots/")

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
				axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),
				axis.ticks=element_line(color="black"),axis.title=element_text(color="black"),
				legend.position="right",axis.line=element_line(color="black")) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
				labels=c("20%","40%","60%")) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1, 20, 40, 60)) +
			scale_color_manual("",values=c("violetred2","black","dodgerblue3","cyan3"))
		ggsave(paste(x,y,"CG_metaplot.pdf",sep="_"),plot,path="new_plots/")

		plot <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
			geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
			theme(panel.background=element_blank(),panel.grid=element_blank(),
				axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),
				axis.ticks=element_line(color="black"),axis.title=element_text(color="black"),
				legend.position="right",axis.line=element_line(color="black")) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,0.3),expand=c(0,0),breaks=c(0.1,0.2,0.3),
				labels=c("10%","20%","30%")) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1, 20, 40, 60)) +
			scale_color_manual("",values=c("violetred2","black","dodgerblue3","cyan3"))
		ggsave(paste(x,y,"CHG_metaplot.pdf",sep="_"),plot,path="new_plots/")

		plot <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
			geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
			theme(panel.background=element_blank(),panel.grid=element_blank(),
				axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),
				axis.ticks=element_line(color="black"),axis.title=element_text(color="black"),
				legend.position="right",axis.line=element_line(color="black")) +
			ylab("Percent methylation") + xlab("") +
			scale_y_continuous(limits=c(0,0.075),expand=c(0,0),breaks=c(0.025,0.05,0.075),
				labels=c("2.5%","5.0%","7.5%")) +
			geom_vline(xintercept=20,linetype="longdash",color="grey55") +
			geom_vline(xintercept=40,linetype="longdash",color="grey55") +
			scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1, 20, 40, 60)) +
			scale_color_manual("",values=c("violetred2","black","dodgerblue3","cyan3"))
		ggsave(paste(x,y,"CHH_metaplot.pdf",sep="_"),plot,path="new_plots/")
	}
}
