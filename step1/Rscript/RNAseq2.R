rm(list = ls())
options(stringsAsFactors = F)
exp <- read.table("GSE150621_read_counts.txt.gz",sep = "\t",comment.char = "!",header = T)
library(GEOquery)

library(AnnoProbe)
eset <- getGEO("GSE150621",destdir = ".",AnnotGPL = F,getGPL = F)
str(eset)
names(eset)
b <- eset[[1]]
phe=pData(b)# clinical infoemation
exp_newdata <- exp[!duplicated(exp$refGene),]
rownames(exp_newdata) <- exp_newdata$refGene
exp_newdata <- exp_newdata[,-1]
exp = exp_newdata[apply(exp_newdata, 1, function(x) sum(x > 1) > 7), ]
dim(exp)
exp[1:4,1:4]
group_list=as.factor(c(rep("control",8),rep("GDM",6)))
library(edgeR)
express_cpm <- log2(cpm(exp)+1)
express_cpm[1:6,1:6]
dat <- express_cpm
dat1 <- express_cpm[,-c(1,3,5,7,8,13)]
tmp <- sort(apply(dat1,1, mad),decreasing = T)[1:500]
exprSet <-dat1[names(tmp),]
library(corrplot)
dim(exprSet)


M <- cor(exprSet)
g <- corrplot(M,order = "AOE",addCoef.col = "white")

corrplot(M,order = "AOE",type="upper",tl.pos = "d")
corrplot(M,add=TRUE, type="lower", method="number",order="AOE",diag=FALSE,tl.pos="n", cl.pos="n")


group_list1 <- as.factor(c(rep("control",3),rep("GDM",5)))
anno <- data.frame(sampleType=group_list1)
rownames(anno) <- colnames(exprSet)
anno
p <- pheatmap::pheatmap(M,display_numbers = T,annotation_col = anno,fontsize = 11,cellheight = 28,cellwidth = 28,annotation_row = anno,fondsize=15)
p
ggsave("map150621.pdf",p,dpi=600)

#deseq2----
exp=exp[,-c(1,3,5,7,8,13)]
group_list=group_list1
cancer_type="GSE150621"
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


#logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
logFC_cutoff <- 1.5
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
deg <- na.omit(DEG)
table(deg$change)
length(deg$change)
cg = rownames(deg)[deg$change !="NOT"]
library(tinyarray)
h1 = draw_heatmap(dat1[cg,],group_list1,scale_before = T,annotation_legend = T)
h1
ggsave("heatmap150621.pdf",h1,dpi=600)
library(tinyarray)
v=draw_volcano(deg,logFC_cutoff = 1.5)+theme(text = element_text(size=15,face = 'bold'))
ggsave("volcano150621.pdf",v,dpi=600)
library(patchwork)
h1/v
p2 <- cowplot::plot_grid(p,h1, v, nrow = 1, labels = LETTERS[1:3])
p2
p
gg150621 <- deg[cg,]
save(gg150621,file = "gg150621.Rdata")
