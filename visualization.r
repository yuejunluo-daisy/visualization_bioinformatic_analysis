##########
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  library(e1071)
  library(parallel)
  library(preprocessCore)

  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- read.table(mixture_file, header=T, sep="\t", row.names=1,check.names=F)

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  if(max(Y) < 50) {Y <- 2^Y}

  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  X <- (X - mean(X)) / sd(as.vector(X))

  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")

##########
#####Network#####
library(openxlsx)
library(reshape2)
library(corrplot)
library(stringr)

setwd("") 

typename<-gsub(".xlsx","",group_file) 
name<-gsub(".xlsx","",data_file)

  outTab=data.frame()
  
  rt=read.xlsx(paste(name,".xlsx",sep=""),sheet = 1, startRow = 1, colNames = TRUE, rowNames = T,na.strings = "NA")
  head(rt)
  head(rt[1:5,1:5])
  str(rt)
  outTab=data.frame()
  
  #rt1=log2(rt[,3:ncol(rt)]+1)
  #rt=cbind(rt[,1:2],rt1)
  
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    outTab=rbind(outTab,cbind(ID=i,
                              HR=coxSummary$coefficients[,"exp(coe)"],
                              CI_up = coxSummary$conf.int[,"upper .95"],
                              CI_low = coxSummary$conf.int[,"lower .95"],
                              log_rank_p = coxSummary$coefficients[,"Pr(>z)"],
                              fraction = coxSummary$coefficients[,"z"]))
  }
  head(outTab)
  write.table(outTab,file=paste0(name,".csv",sep=""),sep=",",row.names=F,quote=F)

##########
#####Riskscore#####
library(openxlsx)
library(survival)
library("survminer")
library(dplyr)

rtee=read.xlsx(paste(name,".xlsx",sep=""),sheet = 1, startRow = 1, colNames = TRUE, rowNames = T,na.strings = "NA")
rtee<-as.data.frame(rtee)

head(rtee[1:5,1:ncol(rtee)])

rtff<-select(rtee,futime,fustat,riskscore)
head(rtff[1:3,1:3])
for(i in 3:ncol(rtff)){
  rt111<-cbind(rtff[,0],rtff[,1], rtff[,2], rtff[,i])
  rt111<-data.frame(rt111)
  names(rt111) <- c('time','status','gene')
  
  cut2 <- coxph(Surv(time, status) ~ gene, data= rt111)
  cut1 <- cutp(cut2)$gene
  cut1
  
  data.table::setorder(cut1, "p")
  cut1[]
  datac1=cut1[]
  
  ##print(cutp(coxph(Surv(time, status) ~ gene, data=rt111)$gene))
  
  
  cutoffpoint_tcga<-datac1$gene[1]
  cutoffpoint_tcga
  
  rt111$gene <- factor((rt111$gene < cutoffpoint_tcga)*1) 
  Nlow<-sum(rt111$gene==1)
  Nlow
  Nhig<-sum(rt111$gene==0)
  Nhig
  
  sdf <- survdiff(Surv(time,status)~gene, data=rt111)
  
  pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  pvalue
  pvalue = round(pvalue,6)
  pvalue
  
  fit <- survfit(Surv(time, status) ~ gene, rt111)
  summary(fit)   
  
  genename<-colnames(rtff)[i]  
  genename
  
##########
#####survivecurve#####
library(survival)
library("survminer")
library('survMisc')

paste(name,".xlsx",sep="")

rt=read.xlsx(paste(name,".xlsx",sep=""),sheet = 1, startRow = 1, colNames = TRUE, rowNames = T,na.strings = "NA")
for(i in 3:ncol(rt)){
  rt111<-cbind(rt[,0],rt[,1], rt[,2], rt[,i])
  rt111<-data.frame(rt111)
  names(rt111) <- c('time','status','gene')
  cut2 <- coxph(Surv(time, status) ~ gene, data= rt111)
  cut1 <- cutp(cut2)$gene
  datac1=cut1[]
  cutoffpoint<-datac1$gene[1]
  rt111$gene <- factor((rt111$gene < cutoffpoint)*1) 
  Nlow<-sum(rt111$gene==1)
  Nhig<-sum(rt111$gene==0)
  sdf <- survdiff(Surv(time,status)~gene, data=rt111)
  pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  pvalue = round(pvalue,6)
  fit <- survfit(Surv(time, status) ~ gene, rt111)
  summary(fit)     
  genename<-colnames(rt)[i]  
  genename
  
##########
#####GSVA#####
library(GSVA)
library(limma)

rt=read.table(inputFile,sep="t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

geneset<-read.table("metagene.txt",header=T, sep="t",check.names=F)

##########
#####heatmap#####
library(pheatmap)
library(RColorBrewer)

data=read.table(input_File, header=T, sep="\t", check.names=F, row.names=1)
clinc_group=read.table(clinc_group, header=T, sep="\t", check.names=F, row.names=1)
gene_group=read.table(gene_group, header=T, sep="\t", check.names=F, row.names=1)

exp_data=as.data.frame(t(data))
sameSample=intersect(row.names(exp_data), row.names(clinc_group))
expData=exp_data[sameSample,]
cli=clinc_group[sameSample,]

data=cbind(expData, cli)
data=data[order(data$riskscore),]
Type=data[,(ncol(data)+1-ncol(cli)):ncol(data)]
geneType=as.data.frame(gene_group)
data=t(data[,1:ncol(expData)])

pdf(paste(genename," heatmap in dataset",".pdf"),height=8, width=8)
col<- colorRampPalette(c("#313695","#4575B4")
pheatmap(data,
         annotation=Type,
         annotation_row=geneType)
dev.off()