rm(list = ls())   
options(stringsAsFactors = F)
library("impute")
require(GEOquery)
require(Biobase)
eset <- getGEO("GSE88929",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
beta.m2 <- beta.m[,39:108]
pD2 <- pD.all[39:108,]
beta.mm <- beta.m2[,c(-11,-62)]
pDmm <- pD2[c(-11,-62),]
beta3 <- impute.knn(beta.mm)
betadata3 <- beta3$data
betadata3=betadata3+0.00001
pd3 <- pDmm[,c(2,10)]
colnames(pd3) <- c("Sample","Group")
identical(colnames(betadata3),rownames(pd3))
library(ChAMP)
myLoad3=champ.filter(beta = betadata3,pd=pd3)
#myLoad
save(myLoad3,file = "stepfinal.Rdata")
if(F){
  myNorm <- champ.norm(beta=myLoad3$beta,arraytype="450K",cores=3)
  
  #champ.SVD(beta=myNorm,pd=myLoad1$pd)
  dim(myNorm) 
  pD=myLoad3$pd
  dim(pD)
  save(myNorm,pD,file = 'step2_norm_final.Rdata')
}
load(file = 'step2_norm_final.Rdata')
beta.m=myNorm
library(stringr)
group_list=ifelse(str_detect(pD$Group,"sample group: IGDM"),"GDM","control")
group_list=factor(group_list,levels = c("control","GDM"))
dim(beta.m) 
grset=makeGenomicRatioSetFromMatrix(beta.m,what="Beta")
M = getM(grset)
dmp <- dmpFinder(M, pheno=group_list, type="categorical")
dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]
dim(dmpDiff)
myDMP <- champ.DMP(beta = beta.m,pheno=group_list)
head(myDMP[[1]])
shen <- myDMP[[1]]
table(group_list)
shen <- shen[shen$gene !="",]
deltabeta_t <- 0.01
shen$change <- ifelse(abs(shen$deltaBeta)>deltabeta_t,ifelse(shen$deltaBeta>deltabeta_t,"up","down"),"not")
down_gene <- shen$gene[shen$change=="down"]
down_gene <- as.character(down_gene)
up_gene <- shen$gene[shen$change=="up"]
up_gene <- as.character(up_gene)
save(up_gene,down_gene,file = "cgi.Rdata")