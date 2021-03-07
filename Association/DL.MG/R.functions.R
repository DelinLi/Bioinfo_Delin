###from server 230
`GAPIT.Pruning` <-
  function(values,DPP=5000){
    #Object: Manhattan & QQ Plot for genome-wide association study
    #Output: Single figure with PDF format
    #Authors: Zhiwu Zhang
    #Modified by Zhengkui Zhou (zkzhou@126.com) 
    #Last update: Feb 8, 2014
    ##############################################################################################
    
    #No change if below the requirement
    if(length(values)<=DPP)return(c(1:length(values)))
    
    #values= log.P.values
    
    values=sqrt(values)  #This shift the weight a little bit to the low building.
    theMin=min(values)
    theMax=max(values)
    range=theMax-theMin
    interval=range/DPP
    
    ladder=round(values/interval)
    ladder2=c(ladder[-1],0)
    keep=ladder-ladder2
    index=which(keep>0)
    
    return(index)
  }#end of GAPIT.Pruning

`GAPIT.Manhattan` <-
  function(GI.MP = NULL, name.of.trait = "Trait",plot.type = "Genomewise", plot.type2 = "Chromosomewise",
           DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,color1="orangered", color2="navyblue", highliht.sig=F,
           color1.sig="orangered", color2.sig="navyblue",cex.axis=1.5, cex.lab=2.2,cex=1.4,cex.sig=1.8,pvalue.max=NA
  ){
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
      position=order(y0,decreasing = TRUE)
      #index0=GAPIT.Pruning(y0[position],DPP=DPP)
      index=position#[index0]
      x=x0[index]
      y=y0[index]
      z=z0[index]
      
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
        x.sig<-x[y>=bonferroniCutOff]
        y.sig<-y[y>=bonferroniCutOff]
        z.sig<-z[y>=bonferroniCutOff]
        mycols.sig=rep(c(color1.sig,color2.sig),max(z.sig))  ## doulb dolor loop by chromosome 
        
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
      
      ############# QQ-plot#################
      
      P.values=P.values[!is.na(P.values)]
      P.values=P.values[P.values>0]
      P.values=P.values[P.values<=1]
      
      if(length(P.values[P.values>0])<1) return(NULL)
      N=length(P.values)
      DPP=round(DPP/4) #Reduce to 1/4 for QQ plot
      P.values <- P.values[order(P.values)]
      
      #Set up the p-value quantiles
      #print("Setting p_value_quantiles...")
      p_value_quantiles <- (1:length(P.values))/(length(P.values)+1)
      
      log.P.values <- -log10(P.values)
      log.Quantiles <- -log10(p_value_quantiles)
      
      index=GAPIT.Pruning(log.P.values,DPP=DPP)
      log.P.values=log.P.values[index ]
      log.Quantiles=log.Quantiles[index]
      
      # par(mar = c(5,6,5,3))
      
      #Add conficence interval
      N1=length(log.Quantiles)
      ## create the confidence intervals
      c95 <- rep(NA,N1)
      c05 <- rep(NA,N1)
      for(j in 1:N1){
        i=ceiling((10^-log.Quantiles[j])*N)
        if(i==0)i=1
        c95[j] <- qbeta(0.95,i,N-i+1)
        c05[j] <- qbeta(0.05,i,N-i+1)
        #print(c(j,i,c95[j],c05[j]))
      }
      
      #CI Lines
      #plot(log.Quantiles, -log10(c05), xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, axes=FALSE, xlab="", ylab="",col="black")
      #par(new=T)
      #plot(log.Quantiles, -log10(c95), xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, axes=FALSE, xlab="", ylab="",col="black")
      
      #CI shade
      plot(NULL, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",cex.axis=cex.axis,
           mgp=c(3,1.8,0),cex.lab=cex.lab, lty = 1,  lwd = 2,  axes=FALSE, xlab="", ylab="",col="gray")
      index=length(c95):1
      polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col='gray',border=NA)
      
      #Diagonal line
      abline(a = 0, b = 1, col = "red",lwd=2)
      
      #data
      par(new=T)
      #if(plot.style=="FarmCPU"){
      #  plot(log.Quantiles, log.P.values, cex.axis=1.1, cex.lab=1.3, lty = 1,  lwd = 2, col = "Black" ,bty='l', xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste(name.of.trait,sep=""),pch=20)
      #}
      #if(plot.style=="rainbow"){
      plot(log.Quantiles, log.P.values, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), cex.axis=cex.axis, cex.lab=cex.lab, lty = 1,  lwd = 2, col = "Blue" ,
           xlab =expression(Expected~~-log[10](italic(p))),ylab = expression(Observed~~-log[10](italic(p))), 
           main ="",# paste(name.of.trait,sep=""),)
           pch=20,mgp=c(3.2,1,0))
      #}
      
      
    }
    
    #dev.off()
    print("QQ-Plot finished!")
    print(paste("Manhattan & QQ Plot for trait: ", name.of.trait," accomplished successfully!" ,sep = ""))
  }   


###

eRD.Manhattan <- function(GI.MP = NULL, name.of.trait = "Trait",plot.type = "Genomewise", plot.type2 = "Chromosomewise",
                          DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,color1="orangered", color2="navyblue", highliht.sig=F,
                          color1.sig="orangered", color2.sig="navyblue",cex.axis=1.5, cex.lab=2.2,cex=1.4,cex.sig=1.8,pvalue.max=NA
){
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
  
  #Remove all SNPs that have P values between 0 and 1 (not na etc)
  GI.MP <- GI.MP[GI.MP[,3]>0,]
  GI.MP <- GI.MP[GI.MP[,3]<=1,]
  
  GI.MP <- GI.MP[GI.MP[,1]!=0,]
  
  numMarker=nrow(GI.MP)
  bonferroniCutOff=-log10(cutOff)####change
  
  #The model Freq
  GI.MP[,3] <-   GI.MP[,3] 
  
  y.lim <-ifelse(is.na(pvalue.max), as.integer(ceiling(max(GI.MP[,3])))*1.2 ,pvalue.max)  
  
  #y.lim = y.lim+1
  print("The max model freq is")
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
    position=order(y0,decreasing = TRUE)
    #index0=GAPIT.Pruning(y0[position],DPP=DPP)
    index=position#[index0]
    x=x0[index]
    y=y0[index]
    z=z0[index]
    
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
    #plot(y~x,ylab="Model Frequency" , ylim=c(0,y.lim), xaxs="i", yaxs="i" ,xlim=c(0,max(x)*1.01),
    plot(y~x,ylab="Model Frequency" , ylim=c(0,y.lim), xaxs="i", yaxs="i" ,xlim=c(0,max(x)*1.01),
         cex.axis=cex.axis, cex.lab=cex.lab ,col=mycols[z], axes=FALSE, type = "p",
         pch=20,lwd=wd,cex=cex, xlab="Chromosome",main=list(name.of.trait,cex = 2.5),mgp=c(3,1.8,0))
    
    if(highliht.sig){
      par(new=T)
      x.sig<-x[y>=bonferroniCutOff]
      y.sig<-y[y>=bonferroniCutOff]
      z.sig<-z[y>=bonferroniCutOff]
      mycols.sig=rep(c(color1.sig,color2.sig),max(z.sig))  ## doulb dolor loop by chromosome 
      
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
    
    print(paste("Manhattan & QQ Plot for trait: ", name.of.trait," accomplished successfully!" ,sep = ""))
  }
  
}



GAPIT.Manhattan.only <-function(GI.MP = NULL, name.of.trait = "Trait",plot.type = "Genomewise", plot.type2 = "Chromosomewise",
                                DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,color1="orangered", color2="navyblue", highliht.sig=F,
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
    position=order(y0,decreasing = TRUE)
    #index0=GAPIT.Pruning(y0[position],DPP=DPP)
    index=position#[index0]
    x=x0[index]
    y=y0[index]
    z=z0[index]
    
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
      x.sig<-x[y>=bonferroniCutOff]
      y.sig<-y[y>=bonferroniCutOff]
      z.sig<-z[y>=bonferroniCutOff]
      mycols.sig=rep(c(color1.sig,color2.sig),max(z.sig))  ## doulb dolor loop by chromosome 
      
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


`GAPIT.Manhattan` <-
  function(GI.MP = NULL, name.of.trait = "Trait",plot.type = "Genomewise", plot.type2 = "Chromosomewise",
           DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,color1="orangered", color2="navyblue", highliht.sig=F,
           color1.sig="orangered", color2.sig="navyblue",cex.axis=1.5, cex.lab=2.2,cex=1.4,cex.sig=1.8,pvalue.max=NA
  ){
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
      position=order(y0,decreasing = TRUE)
      #index0=GAPIT.Pruning(y0[position],DPP=DPP)
      index=position#[index0]
      x=x0[index]
      y=y0[index]
      z=z0[index]
      
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
        x.sig<-x[y>=bonferroniCutOff]
        y.sig<-y[y>=bonferroniCutOff]
        z.sig<-z[y>=bonferroniCutOff]
        mycols.sig=rep(c(color1.sig,color2.sig),max(z.sig))  ## doulb dolor loop by chromosome 
        
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
      
      ############# QQ-plot#################
      
      P.values=P.values[!is.na(P.values)]
      P.values=P.values[P.values>0]
      P.values=P.values[P.values<=1]
      
      if(length(P.values[P.values>0])<1) return(NULL)
      N=length(P.values)
      DPP=round(DPP/4) #Reduce to 1/4 for QQ plot
      P.values <- P.values[order(P.values)]
      
      #Set up the p-value quantiles
      #print("Setting p_value_quantiles...")
      p_value_quantiles <- (1:length(P.values))/(length(P.values)+1)
      
      log.P.values <- -log10(P.values)
      log.Quantiles <- -log10(p_value_quantiles)
      
      index=GAPIT.Pruning(log.P.values,DPP=DPP)
      log.P.values=log.P.values[index ]
      log.Quantiles=log.Quantiles[index]
      
      # par(mar = c(5,6,5,3))
      
      #Add conficence interval
      N1=length(log.Quantiles)
      ## create the confidence intervals
      c95 <- rep(NA,N1)
      c05 <- rep(NA,N1)
      for(j in 1:N1){
        i=ceiling((10^-log.Quantiles[j])*N)
        if(i==0)i=1
        c95[j] <- qbeta(0.95,i,N-i+1)
        c05[j] <- qbeta(0.05,i,N-i+1)
        #print(c(j,i,c95[j],c05[j]))
      }
      
      #CI Lines
      #plot(log.Quantiles, -log10(c05), xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, axes=FALSE, xlab="", ylab="",col="black")
      #par(new=T)
      #plot(log.Quantiles, -log10(c95), xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, axes=FALSE, xlab="", ylab="",col="black")
      
      #CI shade
      plot(NULL, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",cex.axis=cex.axis,
           mgp=c(3,1.8,0),cex.lab=cex.lab, lty = 1,  lwd = 2,  axes=FALSE, xlab="", ylab="",col="gray")
      index=length(c95):1
      polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col='gray',border=NA)
      
      #Diagonal line
      abline(a = 0, b = 1, col = "red",lwd=2)
      
      #data
      par(new=T)
      #if(plot.style=="FarmCPU"){
      #  plot(log.Quantiles, log.P.values, cex.axis=1.1, cex.lab=1.3, lty = 1,  lwd = 2, col = "Black" ,bty='l', xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste(name.of.trait,sep=""),pch=20)
      #}
      #if(plot.style=="rainbow"){
      plot(log.Quantiles, log.P.values, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), cex.axis=cex.axis, cex.lab=cex.lab, lty = 1,  lwd = 2, col = "Blue" ,
           xlab =expression(Expected~~-log[10](italic(p))),ylab = expression(Observed~~-log[10](italic(p))), 
           main ="",# paste(name.of.trait,sep=""),)
           pch=20,mgp=c(3.2,1,0))
      #}
      
      
    }
    
    #dev.off()
    print("QQ-Plot finished!")
    print(paste("Manhattan & QQ Plot for trait: ", name.of.trait," accomplished successfully!" ,sep = ""))
  }   


###

eRD.Manhattan <- function(GI.MP = NULL, name.of.trait = "Trait",plot.type = "Genomewise", plot.type2 = "Chromosomewise",
                          DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,color1="orangered", color2="navyblue", highliht.sig=F,
                          color1.sig="orangered", color2.sig="navyblue",cex.axis=1.5, cex.lab=2.2,cex=1.4,cex.sig=1.8,pvalue.max=NA
){
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
  
  #Remove all SNPs that have P values between 0 and 1 (not na etc)
  GI.MP <- GI.MP[GI.MP[,3]>0,]
  GI.MP <- GI.MP[GI.MP[,3]<=1,]
  
  GI.MP <- GI.MP[GI.MP[,1]!=0,]
  
  numMarker=nrow(GI.MP)
  bonferroniCutOff=-log10(cutOff)####change
  
  #The model Freq
  GI.MP[,3] <-   GI.MP[,3] 
  
  y.lim <-ifelse(is.na(pvalue.max), as.integer(ceiling(max(GI.MP[,3])))*1.2 ,pvalue.max)  
  
  #y.lim = y.lim+1
  print("The max model freq is")
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
    position=order(y0,decreasing = TRUE)
    #index0=GAPIT.Pruning(y0[position],DPP=DPP)
    index=position#[index0]
    x=x0[index]
    y=y0[index]
    z=z0[index]
    
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
    #plot(y~x,ylab="Model Frequency" , ylim=c(0,y.lim), xaxs="i", yaxs="i" ,xlim=c(0,max(x)*1.01),
    plot(y~x,ylab="Model Frequency" , ylim=c(0,y.lim), xaxs="i", yaxs="i" ,xlim=c(0,max(x)*1.01),
         cex.axis=cex.axis, cex.lab=cex.lab ,col=mycols[z], axes=FALSE, type = "p",
         pch=20,lwd=wd,cex=cex, xlab="Chromosome",main=list(name.of.trait,cex = 2.5),mgp=c(3,1.8,0))
    
    if(highliht.sig){
      par(new=T)
      x.sig<-x[y>=bonferroniCutOff]
      y.sig<-y[y>=bonferroniCutOff]
      z.sig<-z[y>=bonferroniCutOff]
      mycols.sig=rep(c(color1.sig,color2.sig),max(z.sig))  ## doulb dolor loop by chromosome 
      
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
    
    print(paste("Manhattan & QQ Plot for trait: ", name.of.trait," accomplished successfully!" ,sep = ""))
  }
  
}



GAPIT.Manhattan.FDR <-function(GI.MP = NULL, name.of.trait = "Trait",plot.type = "Genomewise", plot.type2 = "Chromosomewise",
                                DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,color1="orangered", color2="navyblue", highliht.sig=F,
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
    position=order(y0,decreasing = TRUE)
    #index0=GAPIT.Pruning(y0[position],DPP=DPP)
    index=position#[index0]
    x=x0[index]
    y=y0[index]
    z=z0[index]
    
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
    plot(y~x,ylab=expression(-log[10](italic(FDR))) , ylim=c(0,y.lim), xaxs="i", yaxs="i" ,xlim=c(0,max(x)*1.01),
         cex.axis=cex.axis, cex.lab=cex.lab ,col=mycols[z], axes=FALSE, type = "p",
         pch=20,lwd=wd,cex=cex, xlab="Chromosome",main=list(name.of.trait,cex = 2.5),mgp=c(3,1.8,0))
    
    if(highliht.sig){
      par(new=T)
      x.sig<-x[y>=bonferroniCutOff]
      y.sig<-y[y>=bonferroniCutOff]
      z.sig<-z[y>=bonferroniCutOff]
      mycols.sig=rep(c(color1.sig,color2.sig),max(z.sig))  ## doulb dolor loop by chromosome 
      
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
