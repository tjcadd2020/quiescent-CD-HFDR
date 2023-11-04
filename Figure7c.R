rm(list = ls())
setwd("/Users/cwn/workdir/project/IBD-px/validation/HMP/new_result_hbi")
library(ggplot2)
library(ggrepel)
library(extrafont)
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
DC_R <- read.table('../../../result-all-CD/3.2 pathway-wilcoxon_test_result-Sib-R.tsv', header = T,comment.char = "$", row.names= 1, sep="\t",check.names = T)
DC_R <- DC_R[which(DC_R$FDR<0.05),]
rownames(DC_R) <- sapply(rownames(DC_R),function(x){gsub('-','.',x)})


i="Remission"
j="humann3"
table_path=paste0('5.1 masslin2-',i,'-',j,'.tsv')
maaslin2_all_results = read.table(table_path, header = T,comment.char = "$", row.names= 1, sep="\t")
maaslin2_results_R<-maaslin2_all_results %>% filter(metadata %in% c('Group')) 
maaslin2_results_R$qval<-p.adjust(maaslin2_results_R$pval, method = 'BH')
maaslin2_results_R <- maaslin2_results_R[order(abs(maaslin2_results_R$coef),decreasing=T),]
maaslin2_results_R$order <- seq(nrow(maaslin2_results_R))
maaslin2_results_R$PWY<-short_name(maaslin2_results_R$feature)
maaslin2_results_R$description<-short_name1(maaslin2_results_R$feature)


final_join <- maaslin2_results_R
final_join$coef <- -final_join$coef
final_join$overlap <- NA
final_join[which(final_join$PWY %in% rownames(DC_R)),14]<-"CD"
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

p <- ggplot(final, aes(x = coef, y = -log10(pval), color = label)) +
  geom_point(aes(size = 1.5),alpha = 1) +
  scale_colour_manual(values  = c('#535c68',"#535c68", '#e74c3c','#2980b9'), limits = c('NA-CD','NA-NA', 'sig-CD', 'sig-NA')) +
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

cd <- subset(final, label == 'sig-CD')
#一种样式，借助 ggrepel 包中的函数
options(ggrepel.max.overlaps = Inf)
p1 <- p + theme(legend.position = 'right') +
  geom_text_repel(data = cd, aes(x = coef, y = -log10(pval), label = cd$Description),
                  size = 4,family = 'Arial',box.padding = unit(1, 'lines'), segment.color = 'black', show.legend = FALSE)
p1

hmp <- subset(final, label == 'sig-NA')
hmp <- hmp[c(2,4,6,7,8,9),]
#一种样式，借助 ggrepel 包中的函数
options(ggrepel.max.overlaps = Inf)
p2 <- p1 + theme(legend.position = 'right') +
  geom_text_repel(data = hmp, aes(x = coef, y = -log10(pval), label = hmp$Description),
                  size = 4,family = 'Arial',box.padding = unit(0.7, 'lines'), segment.color = 'black', show.legend = FALSE)
p2

ggsave('5.2 volcano-masslin2-R.pdf', p2, width = 9, height = 6,device = cairo_pdf)

final <- final[,c(12,13,4,6,11)]
final <- final[which(final$pval<0.05),]
write.table(final,file = '5.2 volcano-masslin2-humann3-R.tsv',row.names = T,col.names = TRUE, sep = '\t', quote = FALSE, na = '')

