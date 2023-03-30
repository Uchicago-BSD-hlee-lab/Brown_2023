#!$HOME/bin/Rscript

HOME <- getwd()

args <- commandArgs(TRUE)
intersect<-read.table(args[1],sep="\t")
out     <-args[2]
id.class<-args[3]

gene.id <- read.table(paste0(HOME,"/geneIDs.WS230"), sep="," );
colnames(gene.id)<-c("WBGene","name","cosmid")

intersect.sense<-intersect[intersect[,6]==intersect[,12],]
sense.ppm<-as.data.frame(tapply(intersect.sense[,5],intersect.sense[,10],sum))
colnames(sense.ppm)<-"ppm"
sense.ppm<-merge(sense.ppm,gene.id,by.x=0,by.y=id.class,all.x=T)
colnames(sense.ppm)[0]<-id.class

x <- sum(sense.ppm$ppm)
write(x, stdout())