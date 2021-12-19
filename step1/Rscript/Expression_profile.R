rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(AnnoProbe)
eset <- getGEO("GSE87295",destdir = ".",AnnotGPL = F,getGPL = F)
exp <- exprs(eset[[1]])
exp[1:4,1:4]
exp = log2(exp+1)
pd <- pData(eset[[1]])
p = identical(rownames(pd),colnames(exp));p
gpl <- eset[[1]]@annotation
library(stringr)
group_list=c(rep("GDM",times=5),rep("control",times=5))
group_list
group_list = factor(group_list,
                    levels = c("control","GDM"))
a <- read.csv(file = "GPL10558-50081.txt",comment.char = "#",header = T,check.names = FALSE,sep = "\t") 
b <- a[,c("ID","Symbol")]
colnames(b) <- c("probe_id","symbol")
c <- b[!b$symbol=="permuted_negative",]
ids <- c
dat1 <- exp[,-c(6,7,8)]
tmp <- sort(apply(dat1,1, mad),decreasing = T)[1:500]
exprSet <-dat1[names(tmp),]

library(corrplot)
dim(exprSet)

M <- cor(exprSet)
g <- corrplot(M,order = "AOE",addCoef.col = "white")

group_list=c(rep("GDM",times=5),rep("control",times=2))
group_list
group_list = factor(group_list,
                    levels = c("control","GDM"))
anno <- data.frame(sampleType=group_list)
rownames(anno) <- colnames(exprSet)
anno
p <- pheatmap::pheatmap(M,display_numbers = T,annotation_col = anno,fontsize = 11,cellheight = 28,cellwidth = 28,annotation_row = anno)
p
ggsave("map87295.pdf",p,dpi=600)

library(limma)
design=model.matrix(~group_list)
fit=lmFit(dat1,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)

deg <- inner_join(deg,ids,by="probe_id")

deg <- deg[!duplicated(deg$symbol),]

logFC_t=1
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)
table(deg$change) 
cg = deg$probe_id[deg$change !="stable"]
library(tinyarray)
h1 = draw_heatmap(dat1[cg,],group_list,scale_before = T,annotation_legend = T)
h1
ggsave("heatmap87295.pdf",h1,dpi=600)
colnames(deg)[1]='log2FoldChange'
colnames(deg)[4]='pvalue'
deg$AveExpr=deg$log2FoldChange

v=draw_volcano(deg,logFC_cutoff = 1)+theme(text = element_text(size=15,face = 'bold'))
v
ggsave("volcano87295.pdf",v,dpi=600)


dd <- get_deg(exp = dat1,group_list = group_list,ids = ids)
table(dd$change)

cg = rownames(dd)[dd$change !="stab"]
h1 = draw_heatmap(express_cpm[cg,],group_list,scale_before = T,annotation_legend = T) 

ddd <- get_deg_all(exp = dat1,group_list = group_list,ids = ids)
aa <- ddd$cgs
ddd$plots
save(aa,file = "87295.Rdata")