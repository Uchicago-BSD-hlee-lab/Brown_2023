#!$HOME/bin/Rscript

HOME <- getwd()

library(Biostrings)
library(ggplot2)
library(parallel)

args <- commandArgs(T)
library <- args[1]
cores <- args[2]

genome <- readDNAStringSet(paste0(HOME,"/reference/WS230/c_elegans.WS230.genomic.fa"))

mature.bed <- read.table(paste0(HOME,"/reference/WS230/ce_WS230.pirna.bed"),sep="\t",header=F)
file_path <- paste0(HOME,"/termination/",library,"/",library,".precursor.v0.bowtie")
precursor.bed <- read.table(file=file_path,sep="\t",header=F,skipNul = T)
uniq.precursor.bed <- unique(precursor.bed)
print("data loaded")
precursor.df <- NULL
x <- seq(from=1, to=nrow(uniq.precursor.bed))
construct.df <- function(i){
  split <- unlist(strsplit(as.character(uniq.precursor.bed[i,1]),split=","))
  pirna <- split[1]
  mature <- split[2]
  precursor <- split[3]
  masked <- as.numeric(regexpr(pattern=mature,text=precursor))-1
  if (masked>0){
    overhang.5 <- length(unlist(strsplit(precursor,split=""))[1:masked])
    overhang.3 <- nchar(precursor)-(overhang.5+nchar(mature))
  } else if (masked==0&nchar(precursor)>21){
    overhang.5 <- 0
    overhang.3 <- nchar(precursor)-nchar(mature)} else {
      overhang.5 <- 0
      overhang.3 <- 0}
  
  chr <- as.character(uniq.precursor.bed[i,3])
  start <- uniq.precursor.bed[i,4]
  precursor.length <- nchar(as.character(uniq.precursor.bed[i,5]))

  n.reads <- length(which(precursor.bed[,5]==uniq.precursor.bed[i,5]))/(uniq.precursor.bed[i,7]+1)
  row.tmp <- c(pirna,mature,precursor,precursor.length,overhang.5,overhang.3,n.reads)
  prog <- round(i/nrow(uniq.precursor.bed)*100,2)

  row.tmp
}

precursor.ls <- mclapply(FUN=construct.df,X=x,mc.cores = cores)

saveRDS(precursor.ls,file=paste0(HOME,"/termination/",library,"/","precursor.ls"))
precursor.df <- as.data.frame(do.call(rbind,precursor.ls))
colnames(precursor.df) <- c("piRNA","mature.seq","precursor.seq","precursor.length","overhang.5","overhang.3","n.reads")

saveRDS(precursor.df,file=paste0(HOME,"/termination/",library,"/","precursor.df"))

write.table(precursor.df,file=paste0(HOME,"/termination/",library,"/","precursor.txt"),quote=F,row.names=F,sep="\t")