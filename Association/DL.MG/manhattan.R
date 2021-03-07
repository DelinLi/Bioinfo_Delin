args = commandArgs(trailingOnly=TRUE)
prefix <- args[1]
source("~/opt/DL.GWAS/R.functions.R")

Files <- dir("output")[grep("assoc.txt",dir("output"))]
file <- Files[grep(args[1],Files)]
GWAS.tem<-read.delim(paste0("output/",file))
cutoff.snp=0.0001 #0.1/nrow(GWAS.tem)
GWAS.tem<-GWAS.tem[!is.na(GWAS.tem$p_wald),]


png(paste0(args[1],".GWAS.GEMMA.png"), width = 20,height=8,units = "in",res = 200)
#layout(matrix(1:2, 2, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    ##
par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box

GAPIT.Manhattan.only(GWAS.tem[,c(1,3,12)] ,name.of.trait= args[1],plot.type = "Genomewise",cutOff= cutoff.snp,
                     highliht.sig=F,color1="#FFA500", color2="#2F4F4F",#color1="#C0C0C0", color2="#000000", 
                     cex.lab=2.4,color1.sig="#FFA500", color2.sig="#2F4F4F")

dev.off()
