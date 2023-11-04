rm(list = ls())
setwd("/Users/cwn/workdir/project/IBD-px/validation/HMP/new_result_hbi")
library(pheatmap)
#shortname
short_name <- function(taxa){
  final<-c()
  for (i in taxa) {
    x = unlist(strsplit(i,"..",fixed = T))[1]
    final <- c(final,x)
  }
  final <- make.unique(final, sep = "_")#去重
  # print(final)
}
short_name1 <- function(taxa){
  final<-c()
  for (i in taxa) {
    x = unlist(strsplit(i,"..",fixed = T))[2]
    final <- c(final,x)
  }
  final <- make.unique(final, sep = "_")#去重
  # print(final)
}
#读入discovery cohort表格-AR
DC_A <- read.table('../../../result-all-CD/3.2 pathway-wilcoxon_test_result-Sib-A.tsv', header = T,comment.char = "$", row.names= 1, sep="\t",check.names = F)
DC_A <- DC_A[which(DC_A$FDR<0.05),]
rownames(DC_A) <- sapply(rownames(DC_A),function(x){gsub('-','.',x)})


i="Active"
j="humann3"

table_path=paste0('5.1 masslin2-',i,'-',j,'.tsv')
maaslin2_all_results = read.table(table_path, header = T,comment.char = "$", row.names= 1, sep="\t")
maaslin2_results_A<-maaslin2_all_results %>% filter(metadata %in% c('Group')) 
maaslin2_results_A$qval<-p.adjust(maaslin2_results_A$pval, method = 'BH')
maaslin2_results_A <- maaslin2_results_A[order(abs(maaslin2_results_A$coef),decreasing=T),]
maaslin2_results_A$order <- seq(nrow(maaslin2_results_A))
maaslin2_results_A$PWY<-short_name(maaslin2_results_A$feature)
maaslin2_results_A$description<-short_name1(maaslin2_results_A$feature)

final_join <- maaslin2_results_A
final_join$coef <- -final_join$coef
final_join$overlap <- NA
final_join[which(final_join$PWY %in% rownames(DC_A)),14]<-"CD"
final_join$sig <- NA
final_join[which(final_join$pval <0.05),15]<-"sig"
final_join$label <- paste(final_join$sig,final_join$overlap,sep = '-')
table(final_join$label)

cut_off_pvalue=0.05
cut_off_FC=0
final <- final_join
#注释
annot = read.table("../../pwy_annotation.csv", header = T,comment.char = "$", row.names= 1, sep=",")
annot$feature <- sapply(rownames(annot),function(x){gsub('[^A-Za-z0-9]','.',x)})
final <- dplyr::left_join(final,annot,by='feature')
a<-final[which(final$pval<0.05),]

# final[which(final$label %in% c('NA-CD','NA-NA')),18]<-'Others'
# Class = c(`Amide, Amidine, Amine, and Polyamine Degradation`="#4d1344",
#             `Amino Acid Biosynthesis/Degradation`="#bd4a4f",
#             `Aromatic Compound Degradation/Biosynthesis`="#d77150",
#             `Carbohydrates Biosynthesis/Degradation`="#eaaa6d",
#             `Cell Structure Biosynthesis`="#ecf594",
#             `Cofactor, Prosthetic Group, Electron Carrier, and Vitamin Biosynthesis`="#aed214",
#             `Fatty Acid and Lipid Biosynthesis/Degradation`="#1abc9c",
#             `Fermentation to Short-Chain Fatty Acids`="#4ea06b",
#             `Nucleoside and Nucleotide Biosynthesis/Degradation`="#2980b9",
#             `Secondary Metabolite Biosynthesis/Degradation`="#8e44ad",
#             `Energy Metabolism`="#564191",
#             `Others`="#9c9a9b")

p <- ggplot(final, aes(x = coef, y = -log10(pval), color = label)) +
  geom_point(aes(size = 1.5),alpha = 1) +
  scale_colour_manual(values  = c('#535c68',"#535c68", '#e74c3c','#2980b9'), limits = c('NA-CD','NA-NA', 'sig-CD', 'sig-NA')) +
  # scale_colour_manual(values  = Class, limits = names(Class))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  # 辅助线
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
  # geom_vline(xintercept = c(0), color = 'gray', size = 0.3) +
  # geom_hline(yintercept = -log(0.01, 10), color = 'gray', size = 0.3) +
  # xlim(-5,5) +
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  theme(axis.text = element_text(size = 16,colour = 'black',family = 'Arial'), 
        # axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18,colour = 'black',family = 'Arial'),
        axis.title.x  = element_text(size = 18,colour = 'black',family = 'Arial'),
        panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title=element_blank(),
        legend.position="none",
        panel.border = element_rect(fill=NA,color="black", linewidth=2, linetype="solid"))+
  labs(x = '\nCoef', y = '-Log10(p)', color = '')
p
table(final$Surperclass)
cd <- subset(final, label == 'sig-CD')
cd <- cd[c(1,2,3,7,10),]
# scfa <- subset(final, Surperclass %in% c('Fermentation to Short-Chain Fatty Acids',
#                                          'Carbohydrates Biosynthesis/Degradation',
#                                          'Energy Metabolism'))
#一种样式，借助 ggrepel 包中的函数
options(ggrepel.max.overlaps = Inf)
p1 <- p + theme(legend.position = 'right') +
  geom_text_repel(data = cd, aes(x = coef, y = -log10(pval), label = cd$Description),
                  size = 4,family = 'Arial',box.padding = unit(1, 'lines'), segment.color = 'black', show.legend = FALSE)
p1

hmp <- subset(final, label == 'sig-NA')
hmp <- hmp[c(1,5,11,12,13),]
#一种样式，借助 ggrepel 包中的函数
options(ggrepel.max.overlaps = Inf)
p2 <- p1 + theme(legend.position = 'right') +
  geom_text_repel(data = hmp, aes(x = coef, y = -log10(pval), label = hmp$Description),
                  size = 4,family = 'Arial',box.padding = unit(0.7, 'lines'), segment.color = 'black', show.legend = FALSE)
p2

ggsave('5.2 volcano-masslin2-A.pdf', p2, width = 9, height = 6,device = cairo_pdf)

final <- final[,c(12,13,4,6,11)]
final <- final[which(final$pval<0.05),]
write.table(final,file = '5.2 volcano-masslin2-humann3-A.tsv',row.names = T,col.names = TRUE, sep = '\t', quote = FALSE, na = '')
 