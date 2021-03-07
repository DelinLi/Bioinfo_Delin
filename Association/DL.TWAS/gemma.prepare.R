args = commandArgs(trailingOnly=TRUE)

Prefix  <- args[1]

pheno.file <- args[2]

Path <-"~/RNASeq_delin/data/"
Map <-read.delim("~/RNASeq_delin/data//GmV2.Map.txt")
 
Phenotype <-read.delim(paste0("../",pheno.file), check.names = F)
 
Files <- dir(Path)[grep("Quan5.Geno.txt",dir(Path))]

for( file in Files){
  Data <-read.delim(paste0(Path,file))
  Data<-Data[Data$taxa %in% Phenotype[,1] ,]
  Data <- Data[,apply(Data,2,function(x){return(sum(is.na(x))==0)})]
  row.names(Data) <- Data$taxa
  Data.in<-apply(Data[,-1],2,function(x){ return(round(x,2)) })
  
  Geno <- as.data.frame(t(Data.in))
  
  Geno$Major <-"A"
  Geno$Minor <-"T"
  
  write.table(Geno[,c(ncol(Geno)-1,ncol(Geno),1:(ncol(Geno)-2))],paste0("tables/",sub("txt",paste0(Prefix,".txt"),file)),
              row.names=T,col.names = F,sep=",",quote = F)
  
  write.table(Map[Map$Name %in% colnames(Data),c(1,3,2)],paste0("tables/",
                                                                sub("Geno","Map",sub("txt",paste0(Prefix,".txt"),file))) ,
              row.names=F,col.names = F,sep=",",quote = F )
  
  OUT.Pheno <-  merge(Data[,1:2],Phenotype ,by.x="taxa",by.y=colnames(Phenotype)[1],all.x=T)
  row.names(OUT.Pheno) <- OUT.Pheno$taxa
  
  OUT.Pheno.out <- OUT.Pheno[colnames(Geno)[1:(ncol(Geno)-2)],-c(1:2)]
  
  write.table(OUT.Pheno.out,paste0("tables/",sub("Geno","Pheno",
                                                       sub("txt",paste0(Prefix,".txt"),file))) ,
              row.names=F,col.names = F,sep="\t",quote = F )
  
}

write.table(colnames(OUT.Pheno.out),"tables/Triats.txt",sep="\t",col.names = F,quote = F )


