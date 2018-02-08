
####exmaple codes for Venn Diagram
library(VennDiagram)
venn.diagram(list("Day 10"=DEGs$GeneID[!is.na(DEGs$ASAG7_ASAG3.sig) & DEGs$ASAG7_ASAG3.sig== "yes"],
                  "Day 17 "=(DEGs.d17$GeneID[!is.na(DEGs.d17$ASAG8_ASAG4.sig) & DEGs.d17$ASAG8_ASAG4.sig == "yes"]), 
                  "Day 24"= (DEGs.d24$GeneID[!is.na(DEGs.d24$ASAG9_ASAG5.sig) & DEGs.d24$ASAG9_ASAG5.sig == "yes"]),
                  "Day 31" =(DEGs.d31$GeneID[!is.na(DEGs.d31$ASAG10_ASAG6.sig) & DEGs.d31$ASAG10_ASAG6.sig == "yes"]),
                  "Time"=res$GeneID[!is.na(res$DEG.Sig) & res$DEG.Sig=="yes"]),margin = 0.05,
             fill=c("red","blue","yellow","green","brown"),   height = 5000, width = 6000,cat.cex = 1.5, 
             alpha=c(0.5,0.5,0.5,0.5,0.5), cex=2, cat.fontface=4, fontfamily=3, filename="figures/DEGs.TimeSeries.and.FourTimePoint.VennDiagram.tiff")

venn.diagram(list("Day 10"=DEGs$GeneID[!is.na(DEGs$ASAG7_ASAG3.sig) & DEGs$ASAG7_ASAG3.sig== "yes"],
                  "Day 17 "=(DEGs.d17$GeneID[!is.na(DEGs.d17$ASAG8_ASAG4.sig) & DEGs.d17$ASAG8_ASAG4.sig == "yes"]), 
                  "Day 24"= (DEGs.d24$GeneID[!is.na(DEGs.d24$ASAG9_ASAG5.sig) & DEGs.d24$ASAG9_ASAG5.sig == "yes"]),
                  "Day 31" =(DEGs.d31$GeneID[!is.na(DEGs.d31$ASAG10_ASAG6.sig) & DEGs.d31$ASAG10_ASAG6.sig == "yes"]),
                  "Time"=res$GeneID[!is.na(res$DEG.Sig) & res$DEG.Sig=="yes"]), margin = 0.05,
             fill=c("red","blue","yellow","green","brown"),   height = 5000, width = 6000,cat.cex = 1.5, 
             alpha=c(0.5,0.5,0.5,0.5,0.5), cex=2, cat.fontface=4, fontfamily=3, filename="figures/DEGs.TimeSeries.and.FourTimePoints.VennDiagram.tiff")
