args = commandArgs(trailingOnly=TRUE)

Prefix  <- args[1]
pheno.file <- args[2]


library(multtest)
library(gplots)
#library(genetics)
library(compiler) #this library is already installed in R library("scatterplot3d")
library(scatterplot3d)
source("~/opt/gapit_functions.txt")
source("~/opt/emma.txt")


myGM<-read.delim("~/RNASeq_delin/data//GmV2.Map.txt")

myGD.leaf<-read.delim( "~/RNASeq_delin/data/Leaf.Quan5.Geno.txt")
myGD.seed<-read.delim( "~/RNASeq_delin/data/Seed.Quan5.Geno.txt")

Phenotype<-read.delim(paste0("../",pheno.file), check.names = F)

PCA <- read.table("~/RNASeq_delin/SNPs/final/leaf/plink.eigenvec",check.names = F,header = F)

for( i in 2:ncol(Phenotype)){
  
  trait <- colnames(Phenotype)[i]
  
  myY<-Phenotype[  !is.na(Phenotype[,i]) ,c(1,i)]
  myY<-myY[order(as.character(myY[,1])),]

  
  ##leaf
  myGD.tem<-myGD.leaf[myGD.leaf$taxa %in% myY[,1],]
  myGD.tem<-myGD.tem[order(as.character(myGD.tem$taxa)),]
  myGD.tem <- myGD.tem[,apply(myGD.tem,2,function(x){return(sum(is.na(x))==0)})]
  myGM.tem<-myGM[myGM$Name %in% colnames(myGD.tem),]
  myY.tem<- myY[ myY[,1]%in%  myGD.tem$taxa  , ]
  myCV.tem <- PCA[PCA[,1] %in% myY.tem[,1] ,c(1,3,4,5)] 

  myGAPIT <- GAPIT(Y=myY.tem,
                   GD=myGD.tem, GM=myGM.tem,
                   CV=myCV.tem,
                   model= "CMLM",
                   SNP.MAF=0,
                   file.output=F
                   
  )
  write.csv(myGAPIT$GWAS,paste0(Prefix,"_",trait, "_Leaf.CMLM.SNPsPCA.GWAS.Results.csv" ),row.names = F)
  
  
  myFarmCPU  <- GAPIT(Y=myY.tem,
                      GD=myGD.tem, GM=myGM.tem,
                      CV=myCV.tem,
                      model= "FarmCPU",
                      SNP.MAF=0,
                      file.output=F
                      
  )
  write.csv(myFarmCPU$GWAS,paste0(Prefix,"_",trait, "_Leaf.FarmCPU.SNPsPCA.GWAS.Results.csv" ),row.names = F)

  myGAPIT <- GAPIT(Y=myY.tem,
                   GD=myGD.tem, GM=myGM.tem,
                   #CV=myCV.tem,
                   PCA.total=3,
                   model= "CMLM",
                   SNP.MAF=0,
                   file.output=F
                   
  )
  write.csv(myGAPIT$GWAS,paste0(Prefix,"_",trait, "_Leaf.CMLM.GWAS.Results.csv" ),row.names = F)
  
  
  myFarmCPU  <- GAPIT(Y=myY.tem,
                      GD=myGD.tem, GM=myGM.tem,
                      #CV=myCV.tem,
                      PCA.total=3,
                      model= "FarmCPU",
                      SNP.MAF=0,
                      file.output=F
                      
  )
  write.csv(myFarmCPU$GWAS,paste0(Prefix,"_",trait, "_Leaf.FarmCPU.GWAS.Results.csv" ),row.names = F)
  
  ##seed
  myGD.tem<-myGD.seed[myGD.seed$taxa %in% myY[,1],]
  myGD.tem<-myGD.tem[order(as.character(myGD.tem$taxa)),]
  myGD.tem <- myGD.tem[,apply(myGD.tem,2,function(x){return(sum(is.na(x))==0)})]
  myGM.tem<-myGM[myGM$Name %in% colnames(myGD.tem),]
  myY.tem<- myY[ myY[,1]%in%  myGD.tem$taxa  , ]
  myCV.tem <- PCA[PCA[,1] %in% myY.tem[,1] ,c(1,3,4,5)]
  
  myGAPIT <- GAPIT(Y=myY.tem,
                   GD=myGD.tem, GM=myGM.tem,
                   CV=myCV.tem,
                   PCA.total=0,
                   model= "CMLM",
                   SNP.MAF=0,
                   file.output=F
                   
  )
  write.csv(myGAPIT$GWAS,paste0(Prefix,"_",trait, "_Seed.CMLM.SNPsPCA.GWAS.Results.csv" ),row.names = F)
  
  
  myFarmCPU  <- GAPIT(Y=myY.tem,
                      GD=myGD.tem, GM=myGM.tem,
                      CV=myCV.tem,
                      PCA.total=0,
                      model= "FarmCPU",
                      SNP.MAF=0,
                      file.output=F
                      
  )
  write.csv(myFarmCPU$GWAS,paste0(Prefix,"_",trait, "_Seed.FarmCPU.SNPsPCA.GWAS.Results.csv" ),row.names = F)
  
  myGAPIT <- GAPIT(Y=myY.tem,
                   GD=myGD.tem, GM=myGM.tem,
                   #CV=myCV.tem,
                   PCA.total=3,
                   model= "CMLM",
                   SNP.MAF=0,
                   file.output=F
                   
  )
  write.csv(myGAPIT$GWAS,paste0(Prefix,"_",trait, "_Seed.CMLM.GWAS.Results.csv" ),row.names = F)
  
  
  myFarmCPU  <- GAPIT(Y=myY.tem,
                      GD=myGD.tem, GM=myGM.tem,
                      #CV=myCV.tem,
                      PCA.total=3,
                      model= "FarmCPU",
                      SNP.MAF=0,
                      file.output=F
                      
  )
  write.csv(myFarmCPU$GWAS,paste0(Prefix,"_",trait, "_Seed.FarmCPU.GWAS.Results.csv" ),row.names = F)
}
