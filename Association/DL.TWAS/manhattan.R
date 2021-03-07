args = commandArgs(trailingOnly=TRUE)
prefix <- args[1]
source("~/opt/DL.GWAS/R.functions.R")

gemma.files <- dir("./")
gemma.files <- gemma.files[grep("assoc.txt$",gemma.files)]

for(file in gemma.files){
GWAS.tem<-read.delim(file)

name.tem <- sub(  ".assoc.txt","", file)
cutoff.snp=0.05/nrow(GWAS.tem)
#GWAS.tem<-GWAS.tem[!is.na(GWAS.tem$p_wald),]
GWAS.tem$Chr.Num <- as.numeric(as.character(sub("Chr","",sub("Chr0","",GWAS.tem$chr))))

png(paste0(prefix,".",name.tem,".GWAS.GEMMA.png"), width = 20,height=8,units = "in",res = 200)
#layout(matrix(1:2, 2, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    ##
par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box

GAPIT.Manhattan.only(GWAS.tem[,c(13,3,12)] ,name.of.trait= "",plot.type = "Genomewise",cutOff= cutoff.snp,
                     highliht.sig=F,color1="#FFA500", color2="#2F4F4F",#color1="#C0C0C0", color2="#000000",
                     cex.lab=2.4,color1.sig="#FFA500", color2.sig="#2F4F4F")

dev.off()
}



gapit.files <- dir("./")
gapit.files <- gapit.files[grep("csv$",gapit.files)]

for(file in gapit.files){

GWAS.tem<-read.csv(file)
name.tem <- sub(  ".Results.csv","", file)

cutoff.snp=0.05/nrow(GWAS.tem)
GWAS.tem<-GWAS.tem[!is.na(GWAS.tem$P.value),]

png(paste0(prefix,".",name.tem,".GAPIT.png"), width = 20,height=8,units = "in",res = 200)
#layout(matrix(1:2, 2, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    ##
par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box


GAPIT.Manhattan.only(GWAS.tem[,-1] ,name.of.trait= "",plot.type = "Genomewise",cutOff= cutoff.snp,
                     highliht.sig=F,color1="#FFA500", color2="#2F4F4F",#color1="#C0C0C0", color2="#000000",
                     cex.lab=2.4,color1.sig="#FFA500", color2.sig="#2F4F4F")

dev.off()
}


