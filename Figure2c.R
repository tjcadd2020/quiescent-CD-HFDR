rm(list = ls())
gc()
setwd("/Users/cwn/workdir/project/IBD-px/result-all-CD")
library(ggplot2)
library(dplyr)
library(ggpubr)
venn_result <- read.table("2.1.1.2 venn-wilcox-genus.tsv", header=T, row.names = NULL,sep="\t")
siblingR <- read.table("2.1.1 wilcoxon_test_result-Sib-R-genus.tsv", header=T, row.names = 1,sep="\t")
siblingA <- read.table("2.1.1 wilcoxon_test_result-Sib-A-genus.tsv", header=T, row.names = 1,sep="\t")

table_path='genus-rela.tsv'
table0 = read.table(table_path, header = T,comment.char = "$", row.names= 1, sep="\t")

metadata <- read.table("metadata-sib-nonsib-cou.tsv", header=T, row.names = 1,sep="\t")


venn_result1 <- unlist(strsplit(venn_result[which(venn_result$Sibling_A==TRUE & venn_result$Sibling_R==TRUE &
                                                    venn_result$nonRela_A==TRUE & venn_result$nonRela_R==TRUE),]$..values..,
                                split = ','))
venn_result1 <- c("d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Anaerostipes","d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Alistipes")
table <- table0[which(rownames(table0) %in% venn_result1),]
rownames(table) <- c('Anaerostipes','Alistipes')
# rownames(table) <- c('Fusicatenibacter','Dorea')
group <- read.table("metadata-nonRelative-R.csv", header=T, row.names = 1,sep=",")

df <- merge(group,t(table),by='row.names')
#plot
df$Disease <- as.factor(df$Disease)
p <- ggpaired(df, x="Disease", y="Anaerostipes", 
         point.size=2,color="Disease",width=0.5,
         add="jitter",
         line.color = "gray", line.size = 0.5,
         id = 'Row.names',
         palette=c('#C73B78', '#93B4F1'),
         xlab=" ", 
         ylab="Anaerostipes", title = "CD-R vs non-Relative",
         legend.title=" ",show.legend = F) + 
  stat_compare_means(method = "wilcox.test",paired = FALSE, comparisons=list(c("CD", "non-Relative"))) +
  theme(legend.position = 'none')

ggsave(p, filename = '2.2.2.2 Anaerostipes.pdf',
       width = 3, height = 4, useDingbats=FALSE)
 
