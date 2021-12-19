rm(list = ls())
options(stringsAsFactors = F)
load("154377_step1.Rdata")

dim(exp_data)
exp_data <- exp_data[-(1:4), ]
exp=exp_data[apply(exp_data,1,function(x) sum(x>9)>30), ]
aa=c(16,26,27,37,2,35,22,17,12,18,23,24,39,29,19,15,20,45:71,73:77)
aa=sort(aa)
exp1=exp[,aa]
clinical_data1 <- clinical_data[aa,]
library(edgeR)
express_cpm <- log2(cpm(exp1)+1)
express_cpm[1:6,1:6]
group_list=c(rep("control",17),rep("GDM",32))
group_list=factor(group_list,levels = c("control","GDM"))
dat=as.data.frame(t(express_cpm))
library(Rtsne)
dat=unique(dat)
tsne_out <- Rtsne(dat,dims=2,
                  perplexity=10,theta=0.5) # Run TSNE
summary(tsne_out)
a <- tsne_out$Y
library(ggplot2)
tsne_res <- as.data.frame(a)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
library(ggrepel)
group=group_list
p1 <- ggplot(tsne_res,aes(tSNE1,tSNE2,color=group)) + 
  geom_point() + theme_bw() + 
  geom_hline(yintercept = 0,lty=2,col="red") + 
  geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  stat_ellipse(aes(fill=group),type = "norm", geom = "polygon",alpha= 0,color=NA,linetype=2)+
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "tSNE plot",color="group")+
  stat_ellipse(aes(fill=group_list),
               type = "norm", geom = "polygon",alpha= 0.2,color=NA)+
  theme(text = element_text(size=15,face = 'bold'))


p1
ggsave("tsne.pdf",p1,dpi=600)
#deseq2----
cancer_type="GDM_154377"
library(DESeq2)

colData <- data.frame(row.names =colnames(exp1), 
                      condition=group_list)
if(!file.exists(paste0(cancer_type,"dd.Rdata"))){
  dds <- DESeqDataSetFromMatrix(
    countData = exp1,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = paste0(cancer_type,"dd.Rdata"))
}
load(paste0(cancer_type,"dd.Rdata"))

res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$pvalue),] # sort by P value
DEG <- as.data.frame(resOrdered)
head(DEG)
logFC_cutoff <- 1
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
library(tinyarray)
cg = rownames(DEG)[DEG$change !="NOT"]
h1 = draw_heatmap(express_cpm[cg,],group_list,scale_before = T,annotation_legend = T) 

v=draw_volcano(DEG,logFC_cutoff = 1);v
v1 <- v+theme(text = element_text(size=15,face = 'bold'))
ggsave("volcano154377.pdf",v1,dpi=600)
library(patchwork)
h1/v
h2 <- h1+theme(text = element_text(size=15,face = 'bold'))
h2
ggsave("heatmap154377.pdf",h2,dpi=600)
gg_154377 <- DEG[cg,]
save(gg_154377,file = "154377.Rdata")
