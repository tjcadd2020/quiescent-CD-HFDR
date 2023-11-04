library(tidyverse)
library(RColorBrewer)
# install.packages("aplot")
# BiocManager::install("ggtree")
# library(ggtree)
# library(aplot)
library(grid)
# library(ggplotify)
# library(cowplot)
# devtools::install_github("junjunlab/jjAnno")
# library(jjAnno)
library(ggprism)
library(vegan)

rm(list = ls())
gc()
setwd("/Users/cwn/workdir/project/IBD-px/result-all-CD")
metadata <- read.delim('metadata-sib-non-0525.tsv',sep = '\t', row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
otu <- read.delim(paste0('genus.tsv'),sep = '\t',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

otu <-decostand(otu,'total',2)#按列标准化otu

#colSums(otu)#若想后面做成相对丰度的差异比较，可把这两行释放出来即可
dat =otu[,which(colnames(otu) %in% rownames(metadata))]
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
              final <- c(final,x)
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
}
name_all<- short_name(rownames(dat))
name_all[is.na(name_all)] <- "Others1"
# name_all[name_all=="Bacteria;__;__;__;__;__"] <- "Others"
rownames(dat) <- name_all

#根据行求和结果对数据排序
order<-sort(rowSums(dat[,1:ncol(dat)]),index.return=TRUE,decreasing=T)   
#根据列求和结果对表格排序
cc<-dat[order$ix,]
head(cc, n = 3)
##只展示排名前10的物种，之后的算作Others(根据需求改数字)
dd<-rbind(colSums(cc[21:as.numeric(length(rownames(cc))),]),cc[20:1,])
dd[1,] <- dd[1,]+dd[7,]
dd <- dd[c(1:6,8:21),]
head(dd, n = 3)
rownames(dd)[1]<-"Others"
# head(dd, n = 3)
##再与metadata合并

bb<-merge(t(dd),metadata[,c('Disease','State')],
          by = "row.names")
##宽数据变长数据
kk<-tidyr::gather(bb,Phylum,Abundance,-c(Disease,Row.names,State))
##根据Group,Phylum分组运算
hh <- kk %>%
  group_by(State,Phylum) %>%
  dplyr :: summarise(Abundance=sum(Abundance))
head(hh, n = 3)
##更改因子向量的levels
hh$Phylum = factor(hh$Phylum,order = T,levels = row.names(dd))
yanse <-c("#999999","#F781BF","#A65628","#FFFF33","#FF7F00","#984EA3",
          "#4DAF4A","#377EB8","#74D944","#E41A1C","#DA5724","#CE50CA",
          "#D3D93E","#C0717C","#CBD588","#D7C1B1","#5F7FC7","#673770",
          "#3F4921","#CD9BCD","#38333E","#689030","#AD6F3B")#要确保颜色够用，否则会报错

##按物种丰度排序好的堆积柱形图
p1 <- ggplot(hh,aes(x = State,y = Abundance,fill = Phylum)) + 
  geom_bar(position="fill",stat = "identity",width = 0.6) +
  scale_fill_manual(values =c("#999999","#FF7F00","#CBD588","#689030","#4DAF4A","#377EB8","#5F7FC7","#984EA3","#673770",brewer.pal(11, "Spectral"))) +
  labs(x='',y='Abundance')+
  scale_x_discrete(limits = c("CD-A","CD-R","Sibling-A","Sibling-R","non-Relative"))+
  guides(fill=guide_legend(reverse = TRUE))+ 
  ggprism::theme_prism()+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 45))

ggsave('1.5 taxbar-genus.pdf',p1,width = 7.2,
       height = 7.5,dpi = 300)
