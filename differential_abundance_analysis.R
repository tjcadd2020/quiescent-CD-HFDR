package_list = c("ggplot2","ggpubr","dplyr")
###CHECK PACKAGES
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

rm(list = ls())
gc()
setwd("/Users/cwn/workdir/project/IBD-px/result-all-CD")
metadata <- read.table("metadata-sib-nonsib-cou.tsv", header=T, row.names = 1,sep="\t")
#pairs distance
family_sib <- levels(as.factor(metadata[which(metadata$Disease=='Sibling'),]$familyID))
sibling <- metadata[which(metadata$kinship=='sibship' & metadata$Disease != 'non-Relative'),]
family_cou <- levels(as.factor(metadata[which(metadata$Disease=='Cousin'),]$familyID))
cousin <- metadata[which(metadata$kinship=='consin' & metadata$Disease != 'non-Relative'),]
nonR <- metadata[which(metadata$Disease %in% c('CD','non-Relative')),]

table_path='genus-rela.tsv'
table0 = read.table(table_path, header = T,comment.char = "$", row.names= 1, sep="\t")

#sibling-R
siblingR <- sibling[which(sibling$familyID %in% sibling[which(sibling$State=='Remisson'),]$familyID),]
inter <- intersect(rownames(siblingR),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(siblingR,table1)

#sibling-A
siblingA <- sibling[which(sibling$familyID %in% sibling[which(sibling$State=='Activate'),]$familyID),]
inter <- intersect(rownames(siblingA),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(siblingA,table1)

#cousin-R
cousinR <- cousin[which(cousin$familyID %in% cousin[which(cousin$State=='Remisson'),]$familyID),]
inter <- intersect(rownames(cousinR),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(cousinR,table1)

#cousin-A
cousinA <- cousin[which(cousin$familyID %in% cousin[which(cousin$State=='Activate'),]$familyID),]
inter <- intersect(rownames(cousinA),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(cousinA,table1)

#nonReal-R
nonRR <- nonR[which(nonR$State=='Remisson' | is.na(nonR$State)), ]
inter <- intersect(rownames(nonRR),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(nonRR,table1)

#nonReal-A
nonRA <- nonR[which(nonR$State=='Activate' | is.na(nonR$State)), ]
inter <- intersect(rownames(nonRA),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(nonRA,table1)

#CD-R & CD-A
RA <- metadata[complete.cases(metadata$State),]
inter <- intersect(rownames(RA),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(RA,table1)

#Sibling-A & Sibling-R
siblingR <- sibling[which(sibling$familyID %in% sibling[which(sibling$State=='Remisson'),]$familyID),]
siblingA <- sibling[which(sibling$familyID %in% sibling[which(sibling$State=='Activate'),]$familyID),]
RR <- siblingR[-which(complete.cases(siblingR$State)),]
RA <- siblingA[-which(complete.cases(siblingA$State)),]
RR$State <- 'Remisson'
RA$State <- 'Activate'
finalsibling <- rbind(RR,RA)
inter <- intersect(rownames(finalsibling),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(finalsibling,table1)

#Sibling & non-Relative
sib <- sibling[-which(complete.cases(sibling$State)),]
non <- nonR[which(is.na(nonR$State)), ]
sibnon <- rbind(sib,non)
inter <- intersect(rownames(sibnon),colnames(table0))
table1 <- table0[,which(colnames(table0) %in% inter)]
table1 <- t(table1)
dat <- cbind(sibnon,table1)

#head(dat)
##非参数检验,paired
final=data.frame("ASV_ID"=as.character())

for (i in colnames(dat)[(ncol(metadata)+1):length(dat)]){
  sub_dat <- dat[,c(3,which(colnames(dat)==i))]
  sub_dat[is.na(sub_dat)]<-0
  case <- sub_dat[which(sub_dat$Disease == "CD" ),2]
  control <- sub_dat[which(sub_dat$Disease == "Sibling" ),2]
  wilcox_test <- wilcox.test(case, control,exact=FALSE, paired = TRUE, alternative = 'two.sided')
  control_mean <- mean(control)
  case_mean <- mean(case)
  FC_mean <- mean(case-control)
  table=data.frame("ASV_ID"=i,"case_mean"=case_mean,"control_mean"=control_mean,"FC"=FC_mean,"p.value"=wilcox_test$p.value)
  final<-rbind(final,table)
}

##非参数检验，unpaired
final=data.frame("ASV_ID"=as.character())
for (i in colnames(dat)[(ncol(metadata)+1):ncol(dat)]){
    sub_dat <- dat[,c(3,which(colnames(dat)==i))]
    sub_dat[is.na(sub_dat)]<-0
    wilcox_test <- wilcox.test(sub_dat[[i]]~Disease, data = sub_dat,exact=FALSE, paired = FALSE, alternative = 'two.sided')
    case_mean <- mean(sub_dat[which(sub_dat$Disease=='CD'),2])
    control_mean <- mean(sub_dat[which(sub_dat$Disease=='non-Relative'),2])
    table=data.frame("ASV_ID"=i,"case_mean"=case_mean,"control_mean"=control_mean,"FC"=case_mean-control_mean,"p.value"=wilcox_test$p.value)
    final<-rbind(final,table)
}

##非参数检验，unpaired,Active&Remission
final=data.frame("ASV_ID"=as.character())
for (i in colnames(dat)[(ncol(metadata)+1):ncol(dat)]){
  sub_dat <- dat[,c(39,which(colnames(dat)==i))]
  sub_dat[is.na(sub_dat)]<-0
  wilcox_test <- wilcox.test(sub_dat[[i]]~State, data = sub_dat,exact=FALSE, paired = FALSE, alternative = 'two.sided')
  case_mean <- mean(sub_dat[which(sub_dat$State=='Activate'),2])
  control_mean <- mean(sub_dat[which(sub_dat$State=='Remisson'),2])
  table=data.frame("ASV_ID"=i,"active_mean"=case_mean,"remission_mean"=control_mean,"FC(division)"=case_mean/control_mean,"FC(minus)"=case_mean-control_mean,"p.value"=wilcox_test$p.value)
  final<-rbind(final,table)
}

##非参数检验，unpaired,sibActive&sibRemission
final=data.frame("ASV_ID"=as.character())
for (i in colnames(dat)[(ncol(metadata)+1):ncol(dat)]){
  sub_dat <- dat[,c(39,which(colnames(dat)==i))]
  sub_dat[is.na(sub_dat)]<-0
  wilcox_test <- wilcox.test(sub_dat[[i]]~State, data = sub_dat,exact=FALSE, paired = FALSE, alternative = 'two.sided')
  case_mean <- mean(sub_dat[which(sub_dat$State=='Activate'),2])
  control_mean <- mean(sub_dat[which(sub_dat$State=='Remisson'),2])
  table=data.frame("ASV_ID"=i,"active_mean"=case_mean,"remission_mean"=control_mean,"FC"=case_mean-control_mean,"p.value"=wilcox_test$p.value)
  final<-rbind(final,table)
}

##非参数检验，unpaired,sib&non
final=data.frame("ASV_ID"=as.character())
for (i in colnames(dat)[(ncol(metadata)+1):ncol(dat)]){
  sub_dat <- dat[,c(3,which(colnames(dat)==i))]
  sub_dat[is.na(sub_dat)]<-0
  wilcox_test <- wilcox.test(sub_dat[[i]]~Disease, data = sub_dat,exact=FALSE, paired = FALSE, alternative = 'two.sided')
  case_mean <- mean(sub_dat[which(sub_dat$Disease=='Sibling'),2])
  control_mean <- mean(sub_dat[which(sub_dat$Disease=='non-Relative'),2])
  table=data.frame("ASV_ID"=i,"sibling_mean"=case_mean,"non-relative_mean"=control_mean,"FC"=case_mean-control_mean,"p.value"=wilcox_test$p.value)
  final<-rbind(final,table)
}

##整理结果，计算FDR
# final<-final[complete.cases(final),]
# final <- final[which(final$p.value!=1),]
final <- final[complete.cases(final$p.value),]
final$FDR <- p.adjust(final$p.value, method = 'fdr')
final_sign <- final[which(final$FDR< 0.1),]
final_sign1 <- final[which(final$p.value< 0.05),]
write.table(final,file = '2.1.1 wilcoxon_test_result-R&A-genus.tsv',row.names = FALSE,col.names = TRUE, sep = '\t', quote = FALSE, na = '')



