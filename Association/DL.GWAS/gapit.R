args = commandArgs(trailingOnly=TRUE)
prefix <- args[1]
pheno.file <- args[2]

library(multtest)
library(gplots)
#library(genetics)
library(compiler) #this library is already installed in R library("scatterplot3d")
library(scatterplot3d)
source("~/opt/gapit_functions.txt")
source("~/opt/emma.txt")

CV <-read.table("plink.eigenvec")
#The first 4 PCA
CV.all <- CV[,c(1,3:6)]


##pheno
path<-getwd()
setwd("../")
Pheno <- read.delim(pheno.file,check.names = F)
setwd(path)

CV.in <-CV.all[CV.all[,1] %in% Pheno[,1], ]

myY <- Pheno[ Pheno[,1] %in% CV.all[,1] ,]

myG <- read.table(paste0(prefix,".hmp.txt"), head = FALSE,comment.char="&")

cat(paste(prefix, "Samples:",nrow(myY),"SNPs:",nrow(myG)-1),paste0( prefix,".log"),append = T)

pdf(paste0(prefix,"_hist.pdf"))
hist(myY[,2],cex.lab=1.4,main="",xlab="Samples No.")
dev.off()

#CMLM
myGAPIT <- GAPIT(
  G=myG,
  Y=myY,
  CV=CV.in,
  model="CMLM",
  file.output=F
)

write.csv(myGAPIT$GWAS,   paste0(prefix ,"_CMLM.GWAS.result.csv"),row.names = F)

myGAPIT<-NULL

#FarmCPU
myGAPIT <- GAPIT(
  G=myG,
  Y=myY,
  CV=CV.in,
  model="FarmCPU",
  file.output=F
)
myG<-NULL
write.csv(myGAPIT$GWAS,   paste0(prefix ,"_FarmCPU.GWAS.result.csv"),row.names = F)

