rm(list = ls())
gc()
setwd("/Users/cwn/workdir/project/IBD-px/result-all-CD")
library(ggthemes)
library(ggplot2)
library(patchwork)
library(dplyr)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
venn_result <- read.table("2.1.1.2 venn-wilcox-genus.tsv", header=T, row.names = NULL,sep="\t")
siblingR <- read.table("2.1.1 wilcoxon_test_result-Sib-R-genus.tsv", header=T, row.names = 1,sep="\t")
siblingA <- read.table("2.1.1 wilcoxon_test_result-Sib-A-genus.tsv", header=T, row.names = 1,sep="\t")

nonR <- read.table("2.1.1 wilcoxon_test_result-non-R-genus.tsv", header=T, row.names = 1,sep="\t")
nonA <- read.table("2.1.1 wilcoxon_test_result-non-A-genus.tsv", header=T, row.names = 1,sep="\t")

#
venn_result1 <- unlist(strsplit(venn_result[which(venn_result$Sibling_A==T & venn_result$Sibling_R==F &
                                                    venn_result$nonRela_A==T & venn_result$nonRela_R==F),]$..values..,
                                split = ','))
venn_result2 <- unlist(strsplit(venn_result[which(venn_result$Sibling_A==T & venn_result$Sibling_R==T &
                                                    venn_result$nonRela_A==T & venn_result$nonRela_R==T),]$..values..,
                                split = ','))
venn_result3 <- unlist(strsplit(venn_result[which(venn_result$Sibling_A==T & venn_result$Sibling_R==F &
                                                    venn_result$nonRela_A==T & venn_result$nonRela_R==T),]$..values..,
                                split = ','))
venn_result4 <- unlist(strsplit(venn_result[which(venn_result$Sibling_A==F & venn_result$Sibling_R==T &
                                                    venn_result$nonRela_A==T & venn_result$nonRela_R==T),]$..values..,
                                split = ','))
venn_result <- c(venn_result1,venn_result2,venn_result3,venn_result4)
tables_list <- list('CD-R vs Sibling'=siblingR,'CD-A vs Sibling'=siblingA,
                    'CD-A vs non-Relative'=nonA,'CD-R vs non-Relative'=nonR)

fc_final <- data.frame('Taxonomy'= rownames(nonA),row.names =rownames(nonA))
p_final <- data.frame('Taxonomy'= rownames(nonA),row.names =rownames(nonA))
for (i in tables_list) {
  table <- i[which(rownames(i) %in% venn_result),]
  fc<- data.frame(table[,3],row.names = rownames(table))
  p<- data.frame(table[,5],row.names = rownames(table))
  fc_final <- merge(fc_final,fc,by='row.names',all.y = T)
  fc_final <- tibble::column_to_rownames(fc_final,var='Row.names')
  p_final <- merge(p_final,p,by="row.names",all.y = T)
  p_final <- tibble::column_to_rownames(p_final,var='Row.names')
}
colnames(fc_final)[2:5] <- c('CD-R vs Sibling','CD-A vs Sibling','CD-A vs non-Relative','CD-R vs non-Relative') 
colnames(p_final)[2:5] <- c('CD-R vs Sibling','CD-A vs Sibling','CD-A vs non-Relative','CD-R vs non-Relative') 
p_final[is.na(p_final)] <- 1

#breaks
bk1 <- seq(-0.2,-0.02001,by=0.0001)
bk2 <- seq(-0.02,-0.0001,by=0.0001)
bk3 <- seq(0,0.01,by=0.00001)
bk4 <- seq(0.01001,0.2,by=0.001)

genus<-c()
for (i in fc_final$Taxonomy){
  j <- unlist(strsplit(i,split = ';'))[4]
  if(j=='__'){
    genus <- c(genus,i)
  }else{
    genus <- c(genus,unlist(strsplit(j,split = '__'))[2])
  }
}

annot <- data.frame('Order'=genus,row.names = rownames(fc_final))
table(annot$Order)
ann_colors = list(
  Order = c(Bacteroidales="#00818A", 
             Enterobacterales="#BEAED4",
             Lachnospirales="#737b81",
             Monoglobales="#FFFF99",
             Oscillospirales="#386CB0",
             `Peptostreptococcales-Tissierellales`="#ca3e47",
             `Veillonellales-Selenomonadales`="#F7B32D"))

p <- pheatmap(as.matrix((fc_final[,2:5])),scale = "none",cluster_row = T, cluster_col = T, border=NA,
              display_numbers = t(pmt),fontsize_number = 15, number_color = "black",
              cellwidth = 60, cellheight =20,legend=TRUE,annotation_names_col = TRUE,
              color = c(colorRampPalette(colors = c("#0072B2","#339cd6"))(length(bk1)),
                        colorRampPalette(colors = c("#339cd6","white"))(length(bk2)),
                        colorRampPalette(colors = c("white","#f2b88a"))(length(bk3)),
                        colorRampPalette(colors = c("#f2b88a","#D55E00"))(length(bk4))),
              breaks=c(bk1,bk2,bk3,bk4),
              annotation_row = annot,
              show_colnames=T,
              annotation_colors = ann_colors)

pdf("2.2.2.1 wilcox-abundance-heatmap-genus.pdf",width = 20,height = 15)
p
dev.off()
    

