library(fitdistrplus)

setwd("../gene_methylation/")
rm(df1)
for(x in c('upstream','CDS')){
  for(y in c('CG','CHG','CHH')){
    for(i in samples$Sample){
      input=paste(x,paste(i,"all_genes",x,"methylation.txt",sep="_"),sep="/")
      if(exists("df1")){
        df2 <- read.table(input,header=T,sep="\t",row.names=1)[15]
        colnames(df2) <- c(i)
        df1 <- cbind(df1,df2)
        rm(df2)
        df <- read.table(input,header=T,sep="\t",row.names=1)[15]
        colnames(df) <- c(i)
      } else {
        df1 <- read.table(input,header=T,sep="\t",row.names=1)[15]
        colnames(df1) <- c(i)
      }
    }
  }
}


df2 <- melt(df1)
df2$gene <- row.names(df1)
df2 <- reorder(df2,df2$gene)
gen <- c()
for(i in samples$Generation){
  tmp <-rep(i,nrow(df1))
  gen <- c(gen,tmp)
  rm(tmp)
}
df2$generation <- gen
rm(gen)

G1vG5 <- matrix(nrow=0,ncol=7)
G1vG10 <- matrix(nrow=0,ncol=7)
G5vG10 <- matrix(nrow=0,ncol=7)
for(i in names <- unique(df2$gene)){
  tmp <- df2[df2$gene==i,]
  if( nrow(na.omit(tmp[tmp$generation=="G1",])) > 1 & nrow(na.omit(tmp[tmp$generation=="G5",])) > 1){
    tmp2 <- t.test(tmp[tmp$generation=="G1",]$value,tmp[tmp$generation=="G5",]$value)
    tmp3 <- data.frame(tmp2$estimate[1],
                       tmp2$estimate[2],
                       tmp2$estimate[2]/tmp2$estimate[1],
                       tmp2$stderr,
                       tmp2$conf.int[1],
                       tmp2$conf.int[2],
                       tmp2$p.value)
    row.names(tmp3) <- i
    G1vG5 <- rbind(G1vG5,tmp3)
  } else {
    tmp3 <- data.frame(NA,NA,NA,NA,NA,NA,NA)
    row.names(tmp3) <- i
    colnames(tmp3) <- colnames(G1vG5)
    G1vG5 <- rbind(G1vG5,tmp3)
  }
  if( nrow(na.omit(tmp[tmp$generation=="G1",])) > 1 & nrow(na.omit(tmp[tmp$generation=="G10",])) > 1){
    tmp2 <- t.test(tmp[tmp$generation=="G1",]$value,tmp[tmp$generation=="G10",]$value)
    tmp3 <- data.frame(tmp2$estimate[1],
                       tmp2$estimate[2],
                       tmp2$estimate[2]/tmp2$estimate[1],
                       tmp2$stderr,
                       tmp2$conf.int[1],
                       tmp2$conf.int[2],
                       tmp2$p.value)
    row.names(tmp3) <- i
    G1vG10 <- rbind(G1vG10,tmp3)
  } else {
    tmp3 <- data.frame(NA,NA,NA,NA,NA,NA,NA)
    row.names(tmp3) <- i
    colnames(tmp3) <- colnames(G1vG10)
    G1vG10 <- rbind(G1vG10,tmp3)
  }
  if( nrow(na.omit(tmp[tmp$generation=="G5",])) > 1 & nrow(na.omit(tmp[tmp$generation=="G10",])) > 1){
    tmp2 <- t.test(tmp[tmp$generation=="G5",]$value,tmp[tmp$generation=="G10",]$value)
    tmp3 <- data.frame(tmp2$estimate[1],
                       tmp2$estimate[2],
                       tmp2$estimate[2]/tmp2$estimate[1],
                       tmp2$stderr,
                       tmp2$conf.int[1],
                       tmp2$conf.int[2],
                       tmp2$p.value)
    row.names(tmp3) <- i
    G5vG10 <- rbind(G5vG10,tmp3)
  } else {
    tmp3 <- data.frame(NA,NA,NA,NA,NA,NA,NA)
    row.names(tmp3) <- i
    colnames(tmp3) <- colnames(G5vG10)
    G5vG10 <- rbind(G5vG10,tmp3)
  }
}

colnames(G1vG5) <- c("mean1","mean2","foldchange","stderr","conf1","conf2","p-value")
colnames(G1vG10) <- c("mean1","mean2","foldchange","stderr","conf1","conf2","p-value")
colnames(G5vG10) <- c("mean1","mean2","foldchange","stderr","conf1","conf2","p-value")



df1.log <- log(df1)
df1.samples <- colnames(df1)
df1.pca <- prcomp(df1.log,
                  center = TRUE,
                  scale. = TRUE,
                  na.action = na.omit)

ggbiplot(df1.pca, obs.scale = 1, var.scale = 1, 
         groups = ir.species, ellipse = TRUE, 
         circle = TRUE)


