rm(list = ls())
options(stringsAsFactors = F)

library(GEOquery)
library(AnnoProbe)
eset <- getGEO("GSE112168",destdir = ".",AnnotGPL = F,getGPL = F)
str(eset)
names(eset)

pd <- pData(eset[[1]])

count_files = dir("GSE112168_RAW/",pattern = "*.csv.gz$",recursive = T)
ex = function(x){
  result <- read.csv(file.path("GSE112168_RAW/",x),row.names = 1)
  return(result)
}
exp = lapply(count_files,ex)
exp <- do.call(cbind,exp)

group_list=c(rep("control",6),rep("GDM",6))
group_list=factor(group_list,levels = c("control","GDM"))
#deseq2----
cancer_type="GDM_112168"
library(DESeq2)

colData <- data.frame(row.names =colnames(exp), 
                      condition=group_list)
if(!file.exists(paste0(cancer_type,"dd.Rdata"))){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = paste0(cancer_type,"dd.Rdata"))
}
load(paste0(cancer_type,"dd.Rdata"))


res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$pvalue),] 
DEG <- as.data.frame(resOrdered)
head(DEG)
DEG=DEG[1:1739,]
table(is.na(DEG))

logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
logFC_cutoff <- 1
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)

DEG$change=="UP"
up_mirna=rownames(DEG)[DEG$change=="UP"]
down_mirna=rownames(DEG)[DEG$change=="DOWN"]


