args = commandArgs(trailingOnly=TRUE)
prefix <- args[1]
source("~/opt/DL.GWAS/R.functions.R")


GWAS.tem<-read.delim(paste0("output/",args[1],".assoc.txt"))
cutoff.snp=0.05/nrow(GWAS.tem)
GWAS.tem<-GWAS.tem[!is.na(GWAS.tem$p_wald),]


png(paste0(args[1],".GWAS.GEMMA.png"), width = 20,height=8,units = "in",res = 200)
#layout(matrix(1:2, 2, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    ##
par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box

GAPIT.Manhattan.only(GWAS.tem[,c(1,3,12)] ,name.of.trait= "",plot.type = "Genomewise",cutOff= cutoff.snp,
                     highliht.sig=F,color1="#FFA500", color2="#2F4F4F",#color1="#C0C0C0", color2="#000000", 
                     cex.lab=2.4,color1.sig="#FFA500", color2.sig="#2F4F4F")

dev.off()


GWAS.tem<-read.csv(paste0(prefix ,"_CMLM.GWAS.result.csv"))
cutoff.snp=0.05/nrow(GWAS.tem)
GWAS.tem<-GWAS.tem[!is.na(GWAS.tem$P.value),]


png(paste0(prefix,".GWAS.CMLM.png"), width = 20,height=8,units = "in",res = 200)
#layout(matrix(1:2, 2, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    ##
par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box


GAPIT.Manhattan.only(GWAS.tem[,-1] ,name.of.trait= "",plot.type = "Genomewise",cutOff= cutoff.snp,
                     highliht.sig=F,color1="#FFA500", color2="#2F4F4F",#color1="#C0C0C0", color2="#000000", 
                     cex.lab=2.4,color1.sig="#FFA500", color2.sig="#2F4F4F")

dev.off()


GWAS.tem<-read.csv(paste0(prefix ,"_FarmCPU.GWAS.result.csv"))
cutoff.snp=0.05/nrow(GWAS.tem)
GWAS.tem<-GWAS.tem[!is.na(GWAS.tem$P.value),]


png(paste0(prefix,".FarmCPU.CMLM.png"), width = 20,height=8,units = "in",res = 200)
#layout(matrix(1:2, 2, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    ##
par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box

GAPIT.Manhattan.only(GWAS.tem[,-1] ,name.of.trait= "",plot.type = "Genomewise",cutOff= cutoff.snp,
                     highliht.sig=F,color1="#FFA500", color2="#2F4F4F",#color1="#C0C0C0", color2="#000000", 
                     cex.lab=2.4,color1.sig="#FFA500", color2.sig="#2F4F4F")

dev.off()

GAPIT.Manhattan.Highligh <-  function(GI.MP = NULL, name.of.trait = "Trait",plot.type = "Genomewise", plot.type2 = "Chromosomewise",highliht.sig.loci=F,highligh.SNPs=NA,
         DPP=50000,cutOff=0.01,band=5,seqQTN=NULL, highliht.sig=F,  color1="#C0C0C0", color2="#000000", #color1="orangered", color2="navyblue",
         color1.sig="orangered", color2.sig="navyblue",cex.axis=1.5, cex.lab=2.2,cex=1.4,cex.sig=1.8,pvalue.max=NA){
  print(paste("Start to plot figure for trait: ", name.of.trait ,sep = ""))
  if(is.null(GI.MP)) return
  P.values <- as.numeric(GI.MP[,3])
  borrowSlot=4
  GI.MP[,borrowSlot]=0 #Inicial as 0
  if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
  
  index=which(GI.MP[,borrowSlot]==1  & is.na(GI.MP[,3]))
  GI.MP[index,3]=1
  GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))
  
  #Remove all SNPs that do not have a choromosome, bp position and p value(NA)
  GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
  GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
  GI.MP <- GI.MP[!is.na(GI.MP[,3]),]
  GI.MP <- GI.MP[!is.na(GI.MP[,4]),] # snps ID
  #Remove all SNPs that have P values between 0 and 1 (not na etc)
  GI.MP <- GI.MP[GI.MP[,3]>0,]
  GI.MP <- GI.MP[GI.MP[,3]<=1,]
  
  GI.MP <- GI.MP[GI.MP[,1]!=0,]
  
  numMarker=nrow(GI.MP)
  bonferroniCutOff=-log10(cutOff)####change
  
  #Replace P the -log10 of the P-values
  GI.MP[,3] <-  -log10(GI.MP[,3])
  
  y.lim <-ifelse(is.na(pvalue.max), as.integer(ceiling(max(GI.MP[,3])))*1.2 ,pvalue.max)  
  
  #y.lim = y.lim+1
  print("The max -logP vlaue is")
  print(y.lim)
  
  chm.to.analyze <- unique(GI.MP[,1])
  chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
  numCHR= length(chm.to.analyze)
  
  if(plot.type == "Genomewise")
  {
    nchr=max(chm.to.analyze)
    ncycle=ceiling(nchr/band)
    ncolor=band*ncycle
    palette(rainbow(ncolor+1))
    cycle1=seq(1,nchr,by= ncycle)
    thecolor=cycle1
    
    for(i in 2:ncycle){thecolor=c(thecolor,cycle1+(i-1))}
    GI.MP <- GI.MP[order(GI.MP[,2]),]
    GI.MP <- GI.MP[order(GI.MP[,1]),]
    color.vector <- rep(c("orangered","navyblue"),numCHR)
    ticks=NULL
    lastbase=0
    
    #change base position to accumulatives (ticks)
    for (i in chm.to.analyze)
    {
      index=(GI.MP[,1]==i)
      ticks <- c(ticks, lastbase+mean(GI.MP[index,2]))
      GI.MP[index,2]=GI.MP[index,2]+lastbase
      lastbase=max(GI.MP[index,2])
    }
    
    x0 <- as.numeric(GI.MP[,2])
    y0 <- as.numeric(GI.MP[,3])
    z0 <- as.numeric(GI.MP[,1])
    snps.plot  <- as.character(GI.MP[,3])
    position=order(y0,decreasing = TRUE)
    #index0=GAPIT.Pruning(y0[position],DPP=DPP)
    index=position#[index0]
    x=x0[index]
    y=y0[index]
    z=z0[index]
    snps.plot <- snps.plot[index]
    
    #Extract QTN
    QTN=GI.MP[which(GI.MP[,borrowSlot]==1),]
    
    #Draw circles with same size and different thikness
    size=1
    ratio=5
    base=1
    themax=max(y)
    themin=min(y)
    wd=((y-themin+base)/(themax-themin+base))*size*ratio
    s=size-wd/ratio/2
    y.lim <-ifelse(is.na(pvalue.max), as.integer(ceiling(max(GI.MP[,3])))*1.2 ,pvalue.max)   
    #y.lim=27
    #pdf(paste("GWAS.", name.of.trait,".Manhattan-QQ Plot.pdf" ,sep = ""), width = 18,height=4.5)
    #layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(3,1), heights=c(1,1))  ##  
    #par(mar = (c(3,4,2,2)+ 0.5), mgp=c(1.6,1,0))   ##  
    #par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box
    
    mycols=rep(c(color1,color2),max(z))  ## doulb dolor loop by chromosome 
    plot(y~x,ylab=expression(-log[10](italic(p))) , ylim=c(0,y.lim), xaxs="i", yaxs="i" ,xlim=c(0,max(x)*1.01),
         cex.axis=cex.axis, cex.lab=cex.lab ,col=mycols[z], axes=FALSE, type = "p",
         pch=20,lwd=wd,cex=cex, xlab="Chromosome",main=list(name.of.trait,cex = 2.5),mgp=c(3,1.8,0))
    
    if(highliht.sig){
      par(new=T)
      x.sig<-x[y>=bonferroniCutOff ]
      y.sig<-y[y>=bonferroniCutOff ]
      z.sig<-z[y>=bonferroniCutOff ]
      mycols.sig=rep(c(color1.sig,color2.sig),max(z))  ## doulb dolor loop by chromosome 
      
      points( x.sig,y.sig ,type = "p",pch=20,col=mycols.sig[z.sig],cex=cex.sig)
      
    }
    
    if(highliht.sig.loci){
      par(new=T)
      x.sig<-x[y>=bonferroniCutOff & snps.plot %in% highligh.SNPs]
      y.sig<-y[y>=bonferroniCutOff & snps.plot %in% highligh.SNPs]
      z.sig<-z[y>=bonferroniCutOff & snps.plot %in% highligh.SNPs]
      mycols.sig=rep(c(color1.sig,color2.sig),max(z))  ## doulb dolor loop by chromosome 
      
      points( x.sig,y.sig ,type = "p",pch=20,col=mycols.sig[z.sig],cex=cex.sig)
      
    }
    #mtext(name.of.trait, cex=2.5, font.main =1,  side=3, outer=TRUE,line=-1.5) 
    if(!is.null(dim(QTN)))abline(v=QTN[,2], lty = 2, lwd=1.5, col = "grey")
    abline(h=bonferroniCutOff,col="dimgray",lty=2, lwd=1)
    #title(xlab="Chromosome")
    axis(1, at=ticks,tck=-0.01, cex.axis=cex.axis,labels=chm.to.analyze,tick=T, lwd=1.5, padj=-1)    ## lwd: line width tick=T, 
    axis(2, tck=-0.01, cex.axis=cex.axis,lwd=1.5, padj=1)     ## tck=-0.01 let the tck shor            
    box()
    print("Manhattan-Plot.Genomewise finished!")
    
  }   
}

