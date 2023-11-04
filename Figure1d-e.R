rm(list = ls())
gc()
setwd("/Users/cwn/workdir/project/IBD-px/result-all-CD")

siblingR <- read.table("2.1.1 wilcoxon_test_result-sib-R-genus.tsv", header=T, row.names = 1,sep="\t")
siblingR <- siblingR[which(siblingR$FDR<0.1),]
siblingR <- siblingR[which(siblingR$p.value<0.05),]

siblingA <- read.table("2.1.1 wilcoxon_test_result-sib-A-genus.tsv", header=T, row.names = 1,sep="\t")
siblingA <- siblingA[which(siblingA$FDR<0.1),]
siblingA <- siblingA[which(siblingA$p.value<0.05),]

nonR <- read.table("2.1.1 wilcoxon_test_result-non-R-genus.tsv", header=T, row.names = 1,sep="\t")
nonR <- nonR[which(nonR$FDR<0.1),]
nonR <- nonR[which(nonR$p.value<0.05),]

nonA <- read.table("2.1.1 wilcoxon_test_result-non-A-genus.tsv", header=T, row.names = 1,sep="\t")
nonA <- nonA[which(nonA$p.value<0.05),]
nonA <- nonA[which(nonA$p.value<0.05),]


intersect(siblingA$Taxonomy,siblingR$Taxonomy)

# library
library(VennDiagram)
library(ggplot2)
#Make the plot
venn.diagram(
  x = list(
    #Cousin_R=rownames(cousinR),Counsin_A=rownames(cousinA),
    nonRela_A=rownames(nonA),nonRela_R=rownames(nonR),
    Sibling_A=rownames(siblingA),Sibling_R=rownames(siblingR)
  ),
  #category.names = c("Booba (1995)" , "Nekfeu (663)" , "Brassens (471)"),
  filename = '2.1.1.2 venn-wilcox-genus-p-0.05.tiff',
  output = T ,
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff','#4dadd2'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3),alpha('#4dadd2',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135,-135),
  cat.dist = c(0.04, 0.04, 0.04,0.04),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff','#4dadd2'),
  #rotation = 1
)


df_inter <- get.venn.partitions(list(
  #Cousin_R=rownames(cousinR),Counsin_A=rownames(cousinA),
  nonRela_A=rownames(nonA),nonRela_R=rownames(nonR),
  Sibling_A=rownames(siblingA),Sibling_R=rownames(siblingR)))
for (i in 1:nrow(df_inter)) df_inter[i,'..values..'] <- paste(df_inter[[i,'..values..']], collapse = ',')
df_inter <- t(apply(df_inter,1,function(x){unlist(x)}))
write.table(df_inter, '2.1.1.2 venn-wilcox.tsv', row.names = FALSE, sep = '\t', quote = FALSE) 



