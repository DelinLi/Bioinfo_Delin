args = commandArgs(trailingOnly=TRUE)
path<- getwd()
setwd("../")
Pheno <-read.delim(args[1],check.names = F)
setwd(path)
print(paste(args[1],args[2]))
Tem <-read.table(paste0(args[2],".fam"))
OUT<-NULL
Pheno<-Pheno[!is.na(Pheno[,2]),]
Pheno<-as.matrix(Pheno)
for(i in 1:nrow(Tem)){
  if(sum(Pheno[,1] %in% Tem[i,1])==1 ){
    OUT<- rbind(OUT,Pheno[Pheno[,1] %in% Tem[i,1],])
  }else{
    OUT<- rbind(OUT,c(Tem[i,1],rep(-9,ncol(Pheno)-1)))
  }
  
}
Tem <- cbind(Tem[,1:5],OUT[,-1])

write.table(Tem,paste0(args[2],".fam"),sep=" ",quote = F,row.names = F,col.names = F,na ="-9")
write.table(colnames(Pheno)[-1],"Triats.txt",sep="\t",col.names = F,quote = F )
