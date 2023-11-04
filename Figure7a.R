rm(list = ls())
setwd("/Users/cwn/workdir/project/IBD-px/validation/HMP/new_result_hbi")
library(pheatmap)
#shortname
short_name <- function(taxa){
  final<-c()
  for (i in taxa) {
    x = unlist(strsplit(i,"g__"))[2]
    if(is.na(x)){
      x = unlist(strsplit(i,"f__"))[2]
      if(is.na(x)){
        x = unlist(strsplit(i,"o__"))[2]
        if(is.na(x)){
          x = unlist(strsplit(i,"c__"))[2]
          if(is.na(x)){
            x = unlist(strsplit(i,"p__"))[2]
            if(is.na(x)){
              x = unlist(strsplit(i,"d__"))[2]
              if(is.na(x)){
                final <- c(final,"Unassign")
              }else{
                final <- c(final,x)
              }
            }else{
              final <- c(final,x)
            }
          }else{
            final <- c(final,x)
          }
        }else{
          final <- c(final,x)
        }
      }else{
        final <- c(final,x)
      }
    }else{
      final <- c(final,x)
    }
  }
  final <- make.unique(final, sep = "_")#去重
  # print(final)
}
#overlap
A_overlap <- c("Blautia","Lachnospiraceae;__", "Romboutsia", "Anaerostipes","Butyricicoccus","Alistipes", "Lachnospira","Odoribacter","Oscillospirales;__;__",
               "Ruminococcaceae;__","Agathobacter","Dorea","[Eubacterium]_coprostanoligenes_group",
               "Lachnospiraceae_NK4A136_group","Fusicatenibacter","Roseburia",
               "Monoglobus","NK4A214_group","Prevotellaceae;__","UCG-002","Prevotella","Subdoligranulum","Enterobacteriaceae;__","Escherichia","Veillonella","Clostridioides","Proteus")
R_overlap <- c("Faecalibacterium",
               "Dorea","Fusicatenibacter")
both_overlap <- intersect(A_overlap,R_overlap)

i="Remission"
j="genus"
table_path=paste0('5.1 masslin2-',i,'-',j,'-select2.tsv')
maaslin2_all_results = read.table(table_path, header = T,comment.char = "$", row.names= 1, sep="\t")
maaslin2_results_R<-maaslin2_all_results %>% filter(metadata %in% c('Group')) 
maaslin2_results_R$qval<-p.adjust(maaslin2_results_R$pval, method = 'BH')
maaslin2_results_R <- maaslin2_results_R[order(abs(maaslin2_results_R$coef),decreasing=T),]
maaslin2_results_R$order <- seq(nrow(maaslin2_results_R))
#选择|coef|前25，且p<0.05的属
final_R<- maaslin2_results_R[which(maaslin2_results_R$pval<0.05 | maaslin2_results_R$order<=25),]


i="Active"
table_path=paste0('5.1 masslin2-',i,'-',j,'-select2.tsv')
maaslin2_all_results = read.table(table_path, header = T,comment.char = "$", row.names= 1, sep="\t")
maaslin2_results_A<-maaslin2_all_results %>% filter(metadata %in% c('Group')) 
maaslin2_results_A$qval<-p.adjust(maaslin2_results_A$pval, method = 'BH')

maaslin2_results_A <- maaslin2_results_A[order(abs(maaslin2_results_A$coef),decreasing=T),]
maaslin2_results_A$order <- seq(nrow(maaslin2_results_A))
#选择|coef|前25，且p<0.05的属
final_A<- maaslin2_results_A[which(maaslin2_results_A$pval<0.05 | maaslin2_results_A$order<=25),]


#合并
final_join <- dplyr::left_join(maaslin2_results_R,maaslin2_results_A,by="feature")
final_join <- final_join[which(final_join$feature %in% union(final_A$feature,final_R$feature)),]
final_join$name <- short_name(final_join$feature)
final_join <- final_join[,c(22,4,6,11,14,16,21)]
final_join$coef.x <- -final_join$coef.x
final_join$coef.y <- -final_join$coef.y
final_join$overlap <- NA
final_join[which(final_join$name %in% R_overlap),8]<-"CD-R"
final_join[which(final_join$name %in% A_overlap),8]<-"CD-A"
final_join[which(final_join$name %in% both_overlap),8]<-"both"
colnames(final_join)[1:7] <- c("Genus","coef_r","p_r","order_r","coef_a","p_a","order_a")
# final_join$order_r <- 26-final_join$order_r
# final_join$order_r[final_join$order_r<0]<-NA
# final_join$order_a <- 26-final_join$order_a
# final_join$order_a[final_join$order_a<0]<-NA

#heatmap
pmt <- final_join[c(3,6)]
rownames(pmt) <- final_join$Genus
pmt <- round(pmt, digits = 3)
if (!is.null(pmt)){
  sssmt<- pmt< 0.001
  ssmt <- pmt> 0.001& pmt <0.01
    smt <- pmt >0.01& pmt <0.05
    pmt[sssmt] <-'***'
    pmt[ssmt] <-'**'
    pmt[smt] <- '*'
    pmt[!ssmt&!smt&!sssmt]<- ''
} else {
  pmt <- F
}
min(final_join[,c(2,5)],na.rm = T)
max(final_join[,c(2,5)],na.rm = T)

bk1 <- seq(-4.35,-0.0001,by=0.1) #FC
bk2 <- seq(0,2.72,by=0.1)  #FC

annot <- data.frame('overlap'=final_join$overlap,
                    # 'overlap_a'=final_join$overlap_a,
                    # 'order_a'=final_join$order_a,
                    # 'order_r'=final_join$order_r,
                    row.names = rownames(final_join))
ann_colors = list(
  overlap = c(`CD-R`= "#8294C4",
              `CD-A`= "#E49393",
              `both` = "#38ada9")
  # order_a = c("white","#f19066"),
  # order_r = c("white","#f6b93b")
  )
rownames(final_join) <- final_join$Genus
p <- pheatmap(as.matrix(final_join[,c(2,5)]),scale = "none",cluster_row = T, cluster_col = F, 
              treeheight_row = 8, border=NA,
              display_numbers = pmt,
              fontsize_number = 10,
              # number_format = "%.3f",
              number_color = "black",
              cellwidth = 40 , cellheight =10,legend=TRUE,
              annotation_names_col = TRUE,
              color = c(colorRampPalette(colors = c("#000a7a","white"))(length(bk1)),
                        colorRampPalette(colors = c("white","#be3830"))(length(bk2))),
              breaks=c(bk1,bk2),
              annotation_colors = ann_colors,
              annotation_row = annot,
              annotation_row_width=20,
              show_colnames=T,
              angle_col=0,
              fontsize_row=10,
              fontsize_col = 10)

pdf("5.2 heatmap-genus.pdf",width = 15,height = 15)
p
dev.off()

# write.table(final_join,file = '5.2 heatmap-masslin2-genus.tsv',row.names = T,col.names = TRUE, sep = '\t', quote = FALSE, na = '')


