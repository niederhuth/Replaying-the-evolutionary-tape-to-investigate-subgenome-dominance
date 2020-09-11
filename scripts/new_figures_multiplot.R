#!/usr/bin/env Rscript

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(ggplot2)
library(reshape2)
library(scales)

samples <- read.csv('../misc/samples.csv',header=T)
setwd("metaplots/")
controls <- c("TO1000","R500","mock")
colors=c("TO1000"="blue","R500"="red","mock"="black","G1"="orchid","G5"="orchid","G10"="orchid")

#Genes
plots <- c('all_genes','R500_genes','TO1000_genes','R500_syntelogs','R500_non-syntelogs','TO1000_syntelogs',
'TO1000_non-syntelogs')
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

	sample_list <- c("G1",controls)
	plot1 <- ggplot(df1[df1$variable=="CG" & df1$generation %in% sample_list,]) + 
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black"),) + 
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot2 <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.3),expand=c(0,0),breaks=c(0.1,0.2,0.3),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot3 <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.075),expand=c(0,0),breaks=c(0.025,0.05,0.075),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	sample_list <- c("G5",controls)
	plot4 <- ggplot(df1[df1$variable=="CG" & df1$generation %in% sample_list,]) + 
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black"),) + 
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot5 <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.3),expand=c(0,0),breaks=c(0.1,0.2,0.3),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot6 <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.075),expand=c(0,0),breaks=c(0.025,0.05,0.075),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	sample_list <- c("G10",controls)
	plot7 <- ggplot(df1[df1$variable=="CG" & df1$generation %in% sample_list,]) + 
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black"),) + 
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot8 <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.3),expand=c(0,0),breaks=c(0.1,0.2,0.3),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot9 <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.075),expand=c(0,0),breaks=c(0.025,0.05,0.075),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	pdf(paste("new_plots_multiplots/",x,"_metaplot.pdf",sep=""))
		multiplot(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,cols=3)
	dev.off()
	df1 <- matrix(nrow=0,ncol=0)
}

#LTR
plots=c('all_LTRs','R500_LTRs','TO1000_LTRs')
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

	sample_list <- c("G1",controls)
	plot1 <- ggplot(df1[df1$variable=="CG" & df1$generation %in% sample_list,]) + 
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black"),) + 
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0.25,0.5,0.75,1),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot2 <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot3 <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.125),expand=c(0,0),breaks=c(0.025,0.075,0.125),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	sample_list <- c("G5",controls)
	plot4 <- ggplot(df1[df1$variable=="CG" & df1$generation %in% sample_list,]) + 
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black"),) + 
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0.25,0.5,0.75,1),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot5 <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot6 <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.125),expand=c(0,0),breaks=c(0.025,0.075,0.125),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	sample_list <- c("G10",controls)
	plot7 <- ggplot(df1[df1$variable=="CG" & df1$generation %in% sample_list,]) + 
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black"),) + 
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0.25,0.5,0.75,1),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot8 <- ggplot(df1[df1$variable=="CHG" & df1$generation %in% sample_list,]) +		
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.6),expand=c(0,0),breaks=c(0.2,0.4,0.6),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	plot9 <- ggplot(df1[df1$variable=="CHH" & df1$generation %in% sample_list,]) +
		geom_line(aes(x=window,y=value,color=generation,group=sample),size=1) + 
		theme(panel.background=element_blank(),panel.grid=element_blank(),
			axis.text=element_text(color="black",size=8),
			axis.ticks=element_line(color="black"),
			legend.position="None",axis.line=element_line(color="black")) +
		ylab("") + xlab("") +
		scale_y_continuous(limits=c(0,0.125),expand=c(0,0),breaks=c(0.025,0.075,0.125),
			labels=percent) +
		geom_vline(xintercept=20,linetype="longdash",color="grey55") +
		geom_vline(xintercept=40,linetype="longdash",color="grey55") +
		scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"),breaks=c(1,20,40,60)) +
		scale_color_manual("",values=colors)

	pdf(paste("new_plots_multiplots/",x,"_metaplot.pdf",sep=""))
		multiplot(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,cols=3)
	dev.off()
	df1 <- matrix(nrow=0,ncol=0)
}
