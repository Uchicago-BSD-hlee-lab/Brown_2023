#!$HOME/bin/Rscript

HOME <- getwd()

library(ggplot2)
library(ggpubr)
library(scales)
library(ggsci)
library(reshape2)
library(rstatix)
library(gridExtra)

f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# To make small RNA box plots, load in small RNA master file, then compute relevant fold change and plot
# Below is an example of this code in action used to create the plot shown in Figure 2B

data_piRNA <- read.table(paste0(HOME,"/results/master.v0.21u.sense.txt"),sep="\t",header=T)
data_siRNA <- read.table(paste0(HOME,"/results/master.v0.22G.mrna.anti.txt"),sep="\t",header=T)
data_miRNA <- read.table(paste0(HOME,"/results/master.v0.mirna.sense.txt"),sep="\t",header=T)

# arb <- sort(unique(data_siRNA$JB.20180416_RRS.L4440),decreasing = F)[2]
arb <- 0

fold_npp7_piRNA <- log2((data_piRNA$JB.20180416_RRS.npp7+arb)/(data_piRNA$JB.20180416_RRS.L4440+arb))
fold_prp17_piRNA <- log2((data_piRNA$JB.20180416_RRS.prp17+arb)/(data_piRNA$JB.20180416_RRS.L4440+arb))
fold_ints1_piRNA <- log2((data_piRNA$JB.20181018_RRS.ints1+arb)/(data_piRNA$JB.20181018_RRS.L4440_ints+arb))
fold_dic1_piRNA <- log2((data_piRNA$JB.20181018_RRS.dic1+arb)/(data_piRNA$JB.20181018_RRS.L4440_ints+arb))

piRNA <- data.frame(npp7=fold_npp7_piRNA,
                    prp17=fold_prp17_piRNA,
                    ints1=fold_ints1_piRNA,
                    dic1=fold_dic1_piRNA)

fold_npp7_siRNA <- log2((data_siRNA$JB.20180416_RRS.npp7+arb)/(data_siRNA$JB.20180416_RRS.L4440+arb))
fold_prp17_siRNA <- log2((data_siRNA$JB.20180416_RRS.prp17+arb)/(data_siRNA$JB.20180416_RRS.L4440+arb))
fold_ints1_siRNA <- log2((data_siRNA$JB.20181018_RRS.ints1+arb)/(data_siRNA$JB.20181018_RRS.L4440_ints+arb))
fold_dic1_siRNA <- log2((data_siRNA$JB.20181018_RRS.dic1+arb)/(data_siRNA$JB.20181018_RRS.L4440_ints+arb))

siRNA <- data.frame(npp7=fold_npp7_siRNA[which(data_siRNA$WAGO1.target==T | data_siRNA$WAGO9.target==T)],
                    prp17=fold_prp17_siRNA[which(data_siRNA$WAGO1.target==T | data_siRNA$WAGO9.target==T)],
                    ints1=fold_ints1_siRNA[which(data_siRNA$WAGO1.target==T | data_siRNA$WAGO9.target==T)],
                    dic1=fold_dic1_siRNA[which(data_siRNA$WAGO1.target==T | data_siRNA$WAGO9.target==T)])

fold_npp7_miRNA <- log2((data_miRNA$JB.20180416_RRS.npp7+arb)/(data_miRNA$JB.20180416_RRS.L4440+arb))
fold_prp17_miRNA <- log2((data_miRNA$JB.20180416_RRS.prp17+arb)/(data_miRNA$JB.20180416_RRS.L4440+arb))
fold_ints1_miRNA <- log2((data_miRNA$JB.20181018_RRS.ints1+arb)/(data_miRNA$JB.20181018_RRS.L4440_ints+arb))
fold_dic1_miRNA <- log2((data_miRNA$JB.20181018_RRS.dic1+arb)/(data_miRNA$JB.20181018_RRS.L4440_ints+arb))

miRNA <- data.frame(npp7=fold_npp7_miRNA,
                    prp17=fold_prp17_miRNA,
                    ints1=fold_ints1_miRNA,
                    dic1=fold_dic1_miRNA)

df_piRNA <- melt(piRNA)
df_piRNA$type <- "piRNAs"

df_siRNA <- melt(siRNA)
df_siRNA$type <- "22G RNAs"

df_miRNA <- melt(miRNA)
df_miRNA$type <- "miRNAs"

df <- rbind(df_piRNA,df_siRNA)
df <- rbind(df,df_miRNA)

df$variable <- factor(df$variable,
                      levels=c("ints1","dic1","npp7","prp17"),
                      ordered=T)

df$type <- factor(df$type,
                  levels=c("piRNAs","22G RNAs","miRNAs"),
                  ordered=T)

df <- df[which(is.finite(df$value)),]

stat.test <- df %>%
  df_group_by(vars=c("type","variable")) %>%
  t_test(value~0,alternative = "two.sided",mu=0)

stat.test <- stat.test %>% mutate(y.position = c(1,2,3,4,3,4,2,3,3,3,3,4))

stat.test <- stat.test %>% mutate(group1=stat.test$variable)
stat.test <- stat.test %>% mutate(group2=stat.test$variable)


gg_all <- ggplot(data=df,aes(x=variable,y=value))+
  facet_wrap(~type)+
  stat_summary(fun.data=f,geom="boxplot",width=0.5)+
  scale_fill_uchicago(palette="light")+
  geom_hline(yintercept=0,linetype=2)+
  # coord_cartesian(ylim=c(-3,5))+
  theme_bw(base_size = 15)+theme(aspect.ratio = 1,
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 axis.title.x = element_blank(),
                                 strip.background = element_blank(),
                                 text = element_text(family = "sans",color="black"),
                                 legend.position = "o")+
  labs(y="log2(RNAi reads / control reads)")+
  stat_pvalue_manual(stat.test,label="p")+
  scale_x_discrete(label=c("ints-1","dic-1","npp-7","prp-17"))

ggsave(plot=gg_all,
       filename="gg_all.png",
       path=paste0(HOME,"/results/"),
       device="png",
       dpi=300,height=4,width=10)

ggsave(plot=gg_all,
       filename="gg_all.pdf",
       path=paste0(HOME,"/results/"),
       device="pdf",
       dpi=300,height=4,width=10)




# piRNA precursor length histogram shown in Figure S2C

pirna_bed <- read.table(paste0(HOME,"/reference/WS230/ce_WS230.pirna.bed"),sep="\t",header=F)

expand <- function(precursor.df,overhang.5,pirna_bed){
  if (!is.na(overhang.5)){
    precursor.sub <- precursor.df[which(precursor.df$overhang.5==overhang.5),]
  }  else {precursor.sub <- precursor.df}
  precursor.expand <- data.frame(NULL)
  for (i in 1:nrow(precursor.sub)){
    times=as.numeric(as.character(precursor.sub$n.reads[i]))
    lengths <- rep(x=as.numeric(as.character(precursor.sub$precursor.length[i])),times=times)
    if (max(grep(precursor.sub$piRNA[i],pattern="type2"),0)==1){
      type <- "typeII"
    } else {
      type <- "typeI"
    }
    
    types <- rep(x=type,times=times)
    
    chr <- as.character(pirna_bed[which(pirna_bed[,1]==as.character(precursor.sub$piRNA[i])),3])
    chrs <- rep(x=chr,times=times)
    
    df <- data.frame(lengths,types,chrs)
    
    precursor.expand <- rbind(precursor.expand,df)
  }
  print("expanded")
  precursor.expand
}

expand_uniq <- function(precursor.df,overhang.5,pirna_bed){
  if (!is.na(overhang.5)){
    precursor.sub <- precursor.df[which(precursor.df$overhang.5==overhang.5),]
  }  else {precursor.sub <- precursor.df}
  precursor.expand <- data.frame(NULL)
  for (i in 1:nrow(precursor.sub)){
    times=1
    lengths <- rep(x=as.numeric(as.character(precursor.sub$precursor.length[i])),times=times)
    if (max(grep(precursor.sub$piRNA[i],pattern="type2"),0)==1){
      type <- "typeII"
    } else {
      type <- "typeI"
    }
    
    types <- rep(x=type,times=times)
    
    chr <- as.character(pirna_bed[which(pirna_bed[,1]==as.character(precursor.sub$piRNA[i])),3])
    chrs <- rep(x=chr,times=times)
    
    df <- data.frame(lengths,types,chrs)
    
    precursor.expand <- rbind(precursor.expand,df)
  }
  print("expanded")
  precursor.expand
}


wt.df <- readRDS(paste0(HOME,"/termination/JB.20201013_ev/precursor.df"))
ints1.df <- readRDS(paste0(HOME,"/termination/JB.20201013_ints1/precursor.df"))

# wt.expand <- expand(wt.df,2,pirna_bed)
# ints1.expand <- expand(ints1.df,2,pirna_bed)

wt.expand.uniq <- expand_uniq(wt.df,2,pirna_bed)
ints1.expand.uniq <- expand_uniq(ints1.df,2,pirna_bed)

# wt_no.overhang.rq.expand <- expand(wt.df,NA,pirna_bed)
# ints1_no.overhang.rq.expand <- expand(ints1.df,NA,pirna_bed)



expand.df <- data.frame(precursor_length=c(wt.expand.uniq$lengths,
                                           ints1.expand.uniq$lengths),
                        precursor_chr=c(as.character(wt.expand.uniq$chrs),
                                        as.character(ints1.expand.uniq$chrs)),
                        precursor_type=c(as.character(wt.expand.uniq$types),
                                         as.character(ints1.expand.uniq$types)),
                        library=factor(x=c(rep("WT",times=nrow(wt.expand.uniq)),
                                           rep("ints-1",times=nrow(ints1.expand.uniq))),
                                       ordered=T,
                                       levels=c("WT","ints-1")))
ggdens <- function(df){
  ggplot(data=df,aes(x=precursor_length,fill=library,color=library))+
    theme_bw(base_size = 20)+theme(aspect.ratio = 0.5,
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   legend.title = element_blank())+
    geom_density(aes(y=..density..),
                 alpha=0.2)+
    scale_color_manual(values=c("red","blue"))+
    scale_fill_manual(values=c("red","blue"))+
    scale_x_continuous(limits=c(23,100))+
    labs(y="Fraction of sequences",x="Precursor length (nt)")
}

plot_density <- function(expand.df,lib1,lib2,sub){
  df <- expand.df[which(expand.df$library==lib1 | expand.df$library==lib2),]
  gg <- list()
  x <- 1
  if (sub==0){
    ggdens(df)
    
    
  } else {
    groups <- unique(df[,sub])
    for (i in groups){
      df_group <- df[which(df[,sub]==i),]
      gg[[x]] <- ggdens(df_group)
      x <- x+1
      
    }
    ggarrange(plotlist=gg,labels=groups,common.legend = T)
    
  }
  
  
  
  
}

dens <- plot_density(expand.df,"WT","ints-1",sub=0)

ggsave(
  plot=dens,
  filename="dens.pdf",
  path=paste0(HOME,"/results/"),
  device="pdf",
  dpi=300,height=5,width=5
)

ggsave(
  plot=dens,
  filename="dens.png",
  path=paste0(HOME,"/results/"),
  device="png",
  dpi=300,height=5,width=5
)


# define most snpc-4 dependent loci for figure 3A snpc-4 dependent integrator plots

prob.mat <- function(control,exp){
  p1 <- control
  p2 <- exp
  t1 <- sum(p1)
  t2 <- sum(p2)
  
  pc <-(p1+p2)/(t1+t2) 
  pr1 <- (pc^p1)*(pc^p2)*(1-pc^(t1-p1))*(1-pc^(t2-p2))
  
  pc1 <- p1/t1
  pc2 <- p2/t2
  pr2 <- (pc1^p1)*(pc2^p2)*(1-pc1^(t1-p1))*(1-pc2^(t2-p2))
  
  pr <- pr1/pr2
  pr[is.na(pr)] <- 0
  pr
}

arb <- 3

snpc4_log2 <- log2((data_piRNA$JB.20181018_RRS.snpc4+arb)/(data_piRNA$JB.20181018_RRS.L4440+arb))
snpc4_p <- prob.mat((data_piRNA$JB.20181018_RRS.L4440+arb),(data_piRNA$JB.20181018_RRS.snpc4+arb))
snpc4_nlog10 <- -log10(snpc4_p)
alpha <- (0.05/length(snpc4_p))

df <- data.frame(snpc4_log2,
                 snpc4_nlog10)

snpc4_coor <- which(df$snpc4_log2<(-1) & df$snpc4_nlog10>(-log10(alpha)))
df_snpc4 <- df[snpc4_coor,]
snpc4_names <- data_piRNA$Gene.ID[snpc4_coor]

# construct combined piRNA precursor dataframe and compare snpc-4 dependent loci (Figure 3A)

arb <- unique(sort(data_piRNA$JB.20181018_RRS.L4440))[2]


wt_snpc4.df <- readRDS(paste0(HOME,"/termination/JB.20181018_RRS.L4440/precursor.df"))
wt_ints.df <- readRDS(paste0(HOME,"/termination/JB.20181018_RRS.L4440_ints/precursor.df"))
ints1.df <- readRDS(paste0(HOME,"/termination/JB.20181018_RRS.ints1/precursor.df"))
dic1.df <- readRDS(paste0(HOME,"/termination/JB.20181018_RRS.dic1/precursor.df"))
snpc4.df <- readRDS(paste0(HOME,"/termination/JB.20181018_RRS.snpc4/precursor.df"))

precursor.df.conservative.construct <- function(precursor.df,req){
  if (req=="overhang"){
    which.arg <- as.numeric(as.character(precursor.df$overhang.5))==2
  } else if (req=="length"){
    which.arg <- as.numeric(as.character(precursor.df$precursor.length))>21
  } else {
    print("check requirement")
    break
  }
  
  uniq_names <- as.character(unique(precursor.df[which(which.arg),1]))
  precursor.df.conservative <- data.frame(NULL)
  for (i in uniq_names){
    df_tmp <- precursor.df[which(precursor.df$piRNA==i & which.arg),c(1,7)]
    row <- data.frame(piRNA=i,
                      n.reads=sum(as.numeric(as.character(df_tmp$n.reads))))
    precursor.df.conservative <- rbind(precursor.df.conservative,
                                       row)
  }
  
  precursor.df.conservative
}

rq <- "overhang"

wt_snpc4_conservative <- precursor.df.conservative.construct(wt_snpc4.df,req=rq)
wt_ints_conservative <- precursor.df.conservative.construct(wt_ints.df,req=rq)
ints1_conservative <- precursor.df.conservative.construct(ints1.df,req=rq)
dic1_conservative <- precursor.df.conservative.construct(dic1.df,req=rq)
snpc4_conservative <- precursor.df.conservative.construct(snpc4.df,req=rq)

precursor.df <- merge(wt_snpc4_conservative,wt_ints_conservative,by="piRNA",all=T)
precursor.df <- merge(precursor.df,ints1_conservative,by="piRNA",all=T)
precursor.df <- merge(precursor.df,dic1_conservative,by="piRNA",all=T)
precursor.df <- merge(precursor.df,snpc4_conservative,by="piRNA",all=T)

colnames(precursor.df) <- c("piRNA",
                            "reads_wt_snpc4",
                            "reads_wt_ints",
                            "reads_ints1",
                            "reads_dic1",
                            "reads_snpc4")

precursor.df[is.na(precursor.df)] <- 0



mature.df.conservative.construct <- function(mature.df){
  which.arg <- as.numeric(as.character(mature.df$overhang.5))==0 & as.numeric(as.character(mature.df$overhang.3))==0
  uniq_names <- as.character(unique(mature.df[which(which.arg),1]))
  mature.df.conservative <- data.frame(NULL)
  for (i in uniq_names){
    df_tmp <- mature.df[which(mature.df$piRNA==i & which.arg),c(1,7)]
    row <- data.frame(piRNA=i,
                      n.reads=sum(as.numeric(as.character(df_tmp$n.reads))))
    mature.df.conservative <- rbind(mature.df.conservative,
                                    row)
  }
  
  mature.df.conservative
}

wt_snpc4_conservative <- mature.df.conservative.construct(wt_snpc4.df)
wt_ints_conservative <- mature.df.conservative.construct(wt_ints.df)
ints1_conservative <- mature.df.conservative.construct(ints1.df)
dic1_conservative <- mature.df.conservative.construct(dic1.df)
snpc4_conservative <- mature.df.conservative.construct(snpc4.df)

mature.df <- merge(wt_snpc4_conservative,wt_ints_conservative,by="piRNA",all=T)
mature.df <- merge(mature.df,ints1_conservative,by="piRNA",all=T)
mature.df <- merge(mature.df,dic1_conservative,by="piRNA",all=T)
mature.df <- merge(mature.df,snpc4_conservative,by="piRNA",all=T)

colnames(mature.df) <- c("piRNA",
                         "reads_wt_snpc4",
                         "reads_wt_ints",
                         "reads_ints1",
                         "reads_dic1",
                         "reads_snpc4")

mature.df[is.na(mature.df)] <- 0

# determine normalization strategy used in bowtie_alignment.sh script and apply to piRNA read counts here

mature_call <- mature.df[which(mature.df$piRNA=="21ur-1000"),]
master_call <- data_piRNA[which(data_piRNA$Gene.ID=="21ur-1000"),]

wt_snpc4_total <- mature_call$reads_wt_snpc4/master_call$JB.20181018_RRS.L4440*10^6
wt_ints_total <- mature_call$reads_wt_ints/master_call$JB.20181018_RRS.L4440_ints*10^6
ints1_total <- mature_call$reads_ints1/master_call$JB.20181018_RRS.ints1*10^6
dic1_total <- mature_call$reads_dic1/master_call$JB.20181018_RRS.dic1*10^6
snpc4_total <- mature_call$reads_snpc4/master_call$JB.20181018_RRS.snpc4*10^6

precursor.df$rpm_wt_snpc4 <- precursor.df$reads_wt_snpc4*10^6/wt_snpc4_total
precursor.df$rpm_wt_ints <- precursor.df$reads_wt_ints*10^6/wt_ints_total
precursor.df$rpm_ints1 <- precursor.df$reads_ints1*10^6/ints1_total
precursor.df$rpm_dic1 <- precursor.df$reads_dic1*10^6/dic1_total
precursor.df$rpm_snpc4 <- precursor.df$reads_snpc4*10^6/snpc4_total

mature.df$rpm_wt_snpc4 <- mature.df$reads_wt_snpc4*10^6/wt_snpc4_total
mature.df$rpm_wt_ints <- mature.df$reads_wt_ints*10^6/wt_ints_total
mature.df$rpm_ints1 <- mature.df$reads_ints1*10^6/ints1_total
mature.df$rpm_dic1 <- mature.df$reads_dic1*10^6/dic1_total
mature.df$rpm_snpc4 <- mature.df$reads_snpc4*10^6/snpc4_total

precursor <- precursor.df
precursor_snpc4 <- precursor.df[which(precursor.df$piRNA%in%snpc4_names),]


lib_names <- c(rep(paste("All piRNA precursors\nsnpc-4\nn=",length(precursor$rpm_snpc4),sep=""),times=length(precursor$rpm_snpc4)),
               rep(paste("piRNA precursors down in snpc-4\nsnpc-4\nn=",length(precursor_snpc4$rpm_snpc4),sep=""),times=length(precursor_snpc4$rpm_snpc4)),
               rep(paste("All piRNA precursors\ninst-1\nn=",length(precursor$rpm_ints1),sep=""),times=length(precursor$rpm_ints1)),
               rep(paste("piRNA precursors down in snpc-4\ninst-1\nn=",length(precursor_snpc4$rpm_ints1),sep=""),times=length(precursor_snpc4$rpm_ints1)),
               rep(paste("All piRNA precursors\ndic-1\nn=",length(precursor$rpm_dic1),sep=""),times=length(precursor$rpm_dic1)),
               rep(paste("piRNA precursors down in snpc-4\ndic-1\nn=",length(precursor_snpc4$rpm_dic1),sep=""),times=length(precursor_snpc4$rpm_dic1)))

df_box <- data.frame(fold=c(log2((precursor$rpm_snpc4+arb)/(precursor$rpm_wt_snpc4+arb)),
                            log2((precursor_snpc4$rpm_snpc4+arb)/(precursor_snpc4$rpm_wt_snpc4+arb)),
                            log2((precursor$rpm_ints1+arb)/(precursor$rpm_wt_ints+arb)),
                            log2((precursor_snpc4$rpm_ints1+arb)/(precursor_snpc4$rpm_wt_ints+arb)),
                            log2((precursor$rpm_dic1+arb)/(precursor$rpm_wt_ints+arb)),
                            log2((precursor_snpc4$rpm_dic1+arb)/(precursor_snpc4$rpm_wt_ints+arb))),
                     group=factor(lib_names,
                                  levels=unique(lib_names),
                                  ordered=T),
                     lib=factor(c(rep("snpc-4",times=length(lib_names[which(lib_names==unique(lib_names)[1] | 
                                                                              lib_names==unique(lib_names)[2])])),
                                  rep("inst-1",times=length(lib_names[which(lib_names==unique(lib_names)[3] | 
                                                                              lib_names==unique(lib_names)[4])])),
                                  rep("dic-1",times=length(lib_names[which(lib_names==unique(lib_names)[5] | 
                                                                             lib_names==unique(lib_names)[6])]))),
                                levels=c("snpc-4","inst-1","dic-1"),
                                ordered=T),
                     type=factor(c(rep("Precursor",times=length(lib_names[which(lib_names==unique(lib_names)[1] | 
                                                                                  lib_names==unique(lib_names)[2])])),
                                   rep("Precursor",times=length(lib_names[which(lib_names==unique(lib_names)[3] | 
                                                                                  lib_names==unique(lib_names)[4])])),
                                   rep("Precursor",times=length(lib_names[which(lib_names==unique(lib_names)[5] | 
                                                                                  lib_names==unique(lib_names)[6])]))),
                                 levels=c("Precursor"),
                                 ordered=T))

df_box$group_collapse <- ""
df_box$group_collapse[grep(pattern="All piRNA ",x=df_box$group)] <- "All\nprecursors"
df_box$group_collapse[grep(pattern="precursors down ",x=df_box$group)] <- "snpc-4\ndependent\nprecursors"



stat.test <- df_box %>%
  group_by(lib) %>%
  wilcox_test(fold~group_collapse)

stat.test <- stat.test %>% mutate(y.position = c(7,7,7))

box_precursor <- ggplot(data=df_box,aes(x=group_collapse,y=fold,color=group))+
  facet_wrap(~lib)+
  geom_boxplot(outlier.alpha = 0,width=0.5)+
  geom_point(alpha=0.25,size=1,position = position_jitterdodge(jitter.width=1.5))+
  theme_bw(base_size = 20)+theme(aspect.ratio=2,
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 axis.title.x = element_blank(),
                                 legend.position = "o",
                                 strip.background = element_blank())+
  scale_color_manual(values=c("black","red","black","red","black","red","black","red"))+
  coord_cartesian(ylim=c(-8,8))+
  labs(y="Log2(knockdown/control)",
       x="")+
  stat_pvalue_manual(stat.test,label="p",tip.length = 0)

ggsave(box_precursor,
       path=paste0(HOME,"/results/"),
       filename="box_precursor.pdf",
       device="pdf",dpi=300,width=10,height=5)

ggsave(box_precursor,
       path=paste0(HOME,"/results/"),
       filename="box_precursor.png",
       device="png",dpi=300,width=10,height=5)



# type I vs type II piRNA accumulation shown in Figure 3B

int_typeI <- read.table(paste0(HOME,"/results/master.v0.21u.sense.txt"),sep="\t",header=T)
int_typeII <- read.table(paste0(HOME,"/results/master.v0.21u.type2.sense.txt"),sep="\t",header=T)

# arb <- sort(unique(int_typeII$JB.20181018_RRS.L4440.gen2),decreasing=F)[2]
arb <- 0

int_typeI$ints1_fold <- log2((int_typeI$JB.20181018_RRS.ints1+arb)/(int_typeI$JB.20181018_RRS.L4440_ints+arb))
int_typeI$dic1_fold <- log2((int_typeI$JB.20181018_RRS.dic1+arb)/(int_typeI$JB.20181018_RRS.L4440_ints+arb))

int_typeII$ints1_fold <- log2((int_typeII$JB.20181018_RRS.ints1+arb)/(int_typeII$JB.20181018_RRS.L4440_ints+arb))
int_typeII$dic1_fold <- log2((int_typeII$JB.20181018_RRS.dic1+arb)/(int_typeII$JB.20181018_RRS.L4440_ints+arb))

int_typeI_melt <- melt(int_typeI[,c(1,16:17)])
int_typeII_melt <- melt(int_typeII[,c(1,16:17)])

int_typeI_melt$type <- paste("TypeI piRNA\nN=",nrow(int_typeI),sep="")
int_typeII_melt$type <- paste("TypeII piRNA\nN=",nrow(int_typeII),sep="")

int <- rbind(int_typeI_melt,int_typeII_melt)
int$type <- factor(int$type,
                   levels=c(paste("TypeI piRNA\nN=",nrow(int_typeI),sep=""),
                            paste("TypeII piRNA\nN=",nrow(int_typeII),sep="")),
                   ordered=T)

x_ints1_typeI <- int$value[which(int$variable=="ints1_fold" & int$type==paste("TypeI piRNA\nN=",nrow(int_typeI),sep=""))]
x_ints1_typeI <- x_ints1_typeI[-which(x_ints1_typeI=="Inf")]
x_ints1_typeI <- x_ints1_typeI[-which(x_ints1_typeI=="-Inf")]

t_ints1_typeI <- t.test(x_ints1_typeI,mu=0,alternative="less")

x_ints1_typeII <- int$value[which(int$variable=="ints1_fold" & int$type==paste("TypeII piRNA\nN=",nrow(int_typeII),sep=""))]
x_ints1_typeII <- x_ints1_typeII[-which(x_ints1_typeII=="Inf")]
x_ints1_typeII <- x_ints1_typeII[-which(x_ints1_typeII=="-Inf")]

t_ints1_typeII <- t.test(x_ints1_typeII,mu=0,alternative="less")

x_dic1_typeI <- int$value[which(int$variable=="dic1_fold" & int$type==paste("TypeI piRNA\nN=",nrow(int_typeI),sep=""))]
x_dic1_typeI <- x_dic1_typeI[-which(x_dic1_typeI=="Inf")]
x_dic1_typeI <- x_dic1_typeI[-which(x_dic1_typeI=="-Inf")]

t_dic1_typeI <- t.test(x_dic1_typeI,mu=0,alternative="less")

x_dic1_typeII <- int$value[which(int$variable=="dic1_fold" & int$type==paste("TypeII piRNA\nN=",nrow(int_typeII),sep=""))]
x_dic1_typeII <- x_dic1_typeII[-which(x_dic1_typeII=="Inf")]
x_dic1_typeII <- x_dic1_typeII[-which(x_dic1_typeII=="-Inf")]

t_dic1_typeII <- t.test(x_dic1_typeII,mu=0,alternative="less")

gg_int <- ggplot(data=int,aes(x=variable,y=value,fill=variable))+
  facet_wrap(~type)+
  stat_summary(fun.data=f,geom="boxplot",width=0.5)+
  scale_fill_uchicago(palette="light")+
  geom_hline(yintercept=0,linetype=2)+
  # coord_cartesian(ylim=c(-3,5))+
  theme_bw(base_size = 15)+theme(aspect.ratio = 2,
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 axis.title.x = element_blank(),
                                 legend.position = "o")+
  labs(y="log2(integrator RNAi reads / control reads)")+
  scale_x_discrete(label=c("ints-1","dic-1"))

ggsave(plot=gg_int,
       filename="gg_int.png",
       path=paste0(HOME,"/results/"),
       device="png",
       dpi=300,height=5,width=5)

ggsave(plot=gg_int,
       filename="gg_int.pdf",
       path=paste0(HOME,"/results/"),
       device="pdf",
       dpi=300,height=5,width=5)

