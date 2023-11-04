rm(list = ls())
gc()
setwd("/Users/cwn/workdir/project/IBD-px/result-all-CD")
library(vegan)
library(ggplot2)

metadata <- read.delim('metadata-sib-non-0525.tsv',sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
metadata$'sample-id' <- rownames(metadata)
feat_abd <- read.delim('ASV.tsv',sep = '\t',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
feat_abd <- t(feat_abd) 

##alpha diversity
##shannon
shannon <- diversity(feat_abd,index = 'shannon')
##richness
richness <- apply(feat_abd>0,1,sum)
##simpson
simp <- diversity(feat_abd, "simpson")
##invsimpson
invsimp <- diversity(feat_abd, "inv")
##Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(feat_abd, 2) - 1
## Fisher alpha
fisher <- fisher.alpha(feat_abd)
## evenness
S<-specnumber(feat_abd)#统计每个样品内的物种个数
evenness<-shannon/log(S)
Simpson_evenness<-invsimp/S
##faith pd
faith_pd <- read.delim('result-dada2-core-metrics-results/faith_pd_vector/alpha-diversity.tsv',sep = '\t',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
faith_pd$'sample-id' <- rownames(faith_pd)
faith_pd <- left_join(metadata,faith_pd)
###
sib <- faith_pd[which(faith_pd$Disease!='non-Relative'),]
non_rela <- faith_pd[which(faith_pd$Disease=='non-Relative'),]
D <- sib[which(sib$Disease=='CD'),42]
names(D) <- sib[which(sib$Disease=='CD'),1]
C <- sib[which(sib$Disease=='Sibling'),42]
names(C) <- sib[which(sib$Disease=='Sibling'),1]
u <- cbind(C,D)
colnames(u)<-c("Control","CD")
write.csv(u, 'result-dada2-core-metrics-results/alpha_diversity/faith.csv', quote = F)


naCD <- sib[which(sib$State=='CD-R'),]
nau <- u[which(rownames(u) %in% naCD$familyID),]
colnames(nau)<-c("Control","CD-R")
write.csv(nau, 'result-dada2-core-metrics-results/alpha_diversity/faith-inactive.csv', quote = F)

aCD <- sib[which(sib$State=='CD-A'),]
au <- u[which(rownames(u) %in% aCD$familyID),]
colnames(au)<-c("Control","CD-A")
write.csv(au,'result-dada2-core-metrics-results/alpha_diversity/faith-active.csv', quote = F)

#non-relative
non <- faith_pd[which(faith_pd$State=="non-Relative"),]
write.csv(non,'result-dada2-core-metrics-results/alpha_diversity/faith-non.csv', quote = F)

## Plot all
pairs(cbind(shannon, simp, invsimp, unbias.simp, fisher), pch="+", col="blue")

# ## beta diversity defined as gamma/alpha - 1:
# data(dune)
# data(dune.env)
# alpha <- with(metadata, tapply(specnumber(feat_abd), Disease, mean))
# gamma <- with(metadata, specnumber(feat_abd, Disease))
# gamma/alpha - 1

#output
report=cbind(shannon, richness, simp, invsimp, unbias.simp, fisher, evenness, Simpson_evenness,faith_pd)
write.table(report,file="alpha.txt",sep="\t",col.names=NA)


library(dplyr)
report <- tibble::rownames_to_column(data.frame(report))
colnames(report)[1]="sample-id"
new_table<-left_join(report,metadata,by="sample-id")

for(i in 2:10){
  D <- new_table[which(new_table$Disease=='CD'),i]
  names(D) <- new_table[which(new_table$Disease=='CD'),11]
  C <- new_table[which(new_table$Disease=='control'),i]
  names(C) <- new_table[which(new_table$Disease=='control'),11]
  u <- cbind(C,D)
  colnames(u)<-c("CD","Control")
  write.csv(u, paste0('alpha_diversity/',colnames(new_table)[i],'.csv'), quote = F)
  
  
  naCD <- metadata[which(metadata[[9]]==1),]
  nau <- u[which(rownames(u) %in% naCD$familyID),]
  colnames(nau)<-c("CD-R","Control")
  write.csv(nau, paste0('alpha_diversity/',colnames(new_table)[i],'-inactive.csv'), quote = F)
  
  aCD <- metadata[which(metadata[[9]]>1),]
  au <- u[which(rownames(u) %in% aCD$familyID),]
  colnames(au)<-c("CD-A","Control")
  write.csv(au, paste0('alpha_diversity/',colnames(new_table)[i],'-active.csv'), quote = F)
}

