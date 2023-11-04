rm(list = ls())
gc()
setwd("/Users/cwn/workdir/project/IBD-px/result-all-CD")
library(philentropy)
library(vegan)
library(ggplot2)
library('ade4')
library(extrafont)
library(EcolUtils)

metadata <- read.delim('metadata-sib-nonsib.tsv',sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd <- read.delim('ASV.tsv',sep = '\t',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd <- feat_abd[,which(colnames(feat_abd) %in% rownames(metadata))]


##BC
dist = vegdist(t(feat_abd), method = 'bray', correction = 'lingoes', na.rm = TRUE)

##########
##anova
#pairs distance
family_sib <- levels(as.factor(metadata[which(metadata$Disease=='Sibling'),]$familyID))
metadata1 <- metadata[which(metadata$Disease!='non-Relative'),]
siblingR <- metadata1[which(metadata1$familyID %in% metadata1[which(metadata1$State=='Remisson'),]$familyID),]
siblingA <- metadata1[which(metadata1$familyID %in% metadata1[which(metadata1$State=='Activate'),]$familyID),]
siblingR[which(is.na(siblingR$State)),'State']<-'Sibling-R'
siblingR$State1<-'Remission'
siblingA[which(is.na(siblingA$State)),'State']<-'Sibling-A'
siblingA$State1<-'Active'

metadata1<-rbind(siblingA,siblingR)
non_relative_meta <- metadata[which(metadata$Disease=='non-Relative'),]
non_relative_meta$State <- 'non-Relative'
non_relative_meta$State1 <- 'non-Relative'
metadata1<-rbind(metadata1,non_relative_meta)
metadata1[which(metadata1$State=='Activate'),39] <- 'CD-A'
metadata1[which(metadata1$State=='Remisson'),39] <- 'CD-R'
metadata1[which(metadata1$State=='CD-A'),40] <- 'Active'
metadata1[which(metadata1$State=='CD-R'),40] <- 'Remission'
metadata1[which(metadata1$State %in% c('Sibling-A','Sibling-R','non-Relative')),40] <- 'Healthy'
write.table(metadata1,file = 'metadata-sib-non-0525.tsv',row.names = TRUE,col.names = TRUE, sep = '\t', quote = FALSE, na = '')

table1 <- feat_abd[,which(colnames(feat_abd) %in% rownames(metadata1))]
adonis_result <- adonis(t(table1)~State1, metadata1, method  = 'bray', permutations = 999)
adonis_result$aov.tab$R2[1]
adonis_result$aov.tab$`Pr(>F)`[1]


adonis.pair(dist.mat = vegdist(t(table1),method='bray'),Factor = as.factor(metadata1$State),nper = 1000,corr.method = 'fdr')
adonis.pair(dist.mat = dist,Factor = as.factor(metadata1$State),nper = 1000,corr.method = 'fdr')


##pcoa
df.plot1 <- pco.results$li
df.plot1 <- df.plot1[which(rownames(df.plot1) %in% rownames(metadata1)),]
df.plot1 <- merge(df.plot1,metadata1,by='row.names')

pcoa_plot<-ggplot(df.plot1,aes(x=A1,y=A2,colour=State,shape=State1))+
  scale_colour_manual(values=c('CD-A' = '#D44F4D', 'CD-R' = '#279AE7', 'Sibling-A' = '#C1D054','Sibling-R' = '#4bc2c5','non-Relative' = '#999E9D'))+
  #scale_colour_manual(values=c('Remisson' = '#D44F4D', 'Activate' = '#279AE7'))+
  scale_size_area()+
  xlab(axis.1.title)+
  ylab(axis.2.title)+
  stat_ellipse(geom = "polygon",level=0.95,linetype=2,type="norm",alpha = 0.1, aes(fill = State),size=1)+
  scale_fill_manual(values=c('CD-A' = '#D44F4D', 'CD-R' = '#279AE7', 'Sibling-A' = '#C1D054','Sibling-R' = '#4bc2c5','non-Relative' = '#999E9D'))+
  geom_point(alpha =.7, size=2)+
  theme_classic()+
  # ylim(-0.5,0.5)+
  geom_vline(xintercept = 0, color = 'black', size = 0.6, linetype = 4)+
  geom_hline(yintercept = 0, color = 'black', size = 0.6, linetype = 4)+
  theme(axis.text = element_text(size = 12,color = 'black',family = 'Arial',face = 'bold'), 
        axis.title.x = element_text(size = 14,color = 'black',family = 'Arial',face = 'bold'),
        axis.title.y  = element_text(size = 14,color = 'black',family = 'Arial',face = 'bold'),
        panel.grid = element_line(color = 'gray', linetype = 1, size = .6),
        panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title=element_blank(),aspect.ratio=1)

ggsave("1.1 pcoa_BC-AR1.pdf",pcoa_plot,width=6 ,height=7,device = cairo_pdf)
