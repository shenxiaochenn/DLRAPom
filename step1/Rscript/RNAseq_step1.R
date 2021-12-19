rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(tinyarray)
library(AnnoProbe)
eset <- getGEO("GSE154377",destdir = ".",AnnotGPL = F,getGPL = F)
str(eset)
names(eset)
b <- eset[[1]]
phe=pData(b)# clinical information

count_files = dir("GSE154377_RAW/",pattern = "*.txt.gz$",recursive = T)
ex = function(x){
  result <- read.table(file.path("GSE154377_RAW/",x),row.names = 1,sep = "\t")
  return(result)
}
exp = lapply(count_files,ex)
exp <- do.call(cbind,exp)
dim(exp)# matrix
dim(phe)
shen <- rownames(phe)
colnames(exp) <- shen
clinical_data <- phe[1:77,]
exp_data <- exp[,1:77]
library(stringr)
group_list=ifelse(str_detect(clinical_data$title,"Normal"),"control","GDM")
group_list = factor(group_list,
                    levels = c("control","GDM"))
save(exp_data,clinical_data,group_list,file = "154377_step1.Rdata")
