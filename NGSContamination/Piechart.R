setwd("Path for downloaded Files")

Blast.out<-read.delim("Id.txt",header=F) ##Read in the files with beat hits' GI ID (only one column)

Data.blast<-read.delim("IDs.taxi.txt",header=F) ##The file with two column, column1 GI & column2 Taxi_ID, which are unique lines
Data.blast<-merge(Blast.out,Data.blast,by.x="V1",by.y="V1",all.x=T)

Taxi.annotation<-read.table("tax_report.txt",sep="\t",header=T) 
##Taxi ID and it's indentity From  http://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi  
#and processed with "1,$ s/\t|//g" in vim
Taxi.annotation$Species<-gsub("[[:space:]].*","",Taxi.annotation$taxname)
#Genus ID to reduce organisms

Summary<-merge(Data.blast,Taxi.annotation[,c(2,4,5)],by.x="V2",by.y="taxid",all.x=T)
nrow(Summary)
table(Summary[,4])[which.max(table(Summary[,4]))]
Table.out<-table(Summary[,4])[(table(Summary[,4]))>500]##This cutoff (500) is really depends on situation
Table.out<-c(Table.out,sum(table(Summary[,4])[(table(Summary[,4]))<=500]))
names(Table.out)[4]<-"Others"

##paste Genus with it's percentage
lbls <- paste(names(Table.out), "\n", paste(round(Table.out/sum(Table.out),2)*100,"%",sep=""), sep="")

pdf("piechar.pdf")
pie(Table.out, labels = lbls, main="",cex=1.5)
dev.off()

