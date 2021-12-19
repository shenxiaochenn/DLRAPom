rm(list = ls())
options(stringsAsFactors = F)

load("cgi.Rdata")
load("154377.Rdata")
load("gg150621.Rdata")
load("87295.Rdata")
load("WGCNA_sym.Rdata")

up154377 <- rownames(gg_154377)[gg_154377$change=="UP"]
down154377 <- rownames(gg_154377)[gg_154377$change=="DOWN"]
library(clusterProfiler)
library(org.Hs.eg.db)
up154377 <- bitr(up154377,fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
down154377 <- bitr(down154377,fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
up150621 <- rownames(gg150621)[gg150621$change=="UP"]
down150621 <- rownames(gg150621)[gg150621$change=="DOWN"]
up87295 <-aa$deg$up$upgenes
down87925 <- aa$deg$down$downgenes

up_cg_gene <- up_gene
down_cg_gene <- down_gene

library(tinyarray)
library(ggVennDiagram)


x1_up <- list(WGCNA=gene_symbol$SYMBOL,
              GSE87295=up87295,
              GSE154377=up154377$SYMBOL,
              GSE150621=up150621,
              GSE88929=down_cg_gene)
x1_down <- list(WGCNA=gene_symbol$SYMBOL,
                GSE87295=down87925,
                GSE154377=down154377$SYMBOL,
                GSE150621=down150621,
                GSE88929=up_cg_gene)
library(ggplot2)
p <- ggVennDiagram(x1_up,edge_size = 1,label = 'count',edge_lty = "dashed",set_size = 10,label_size = 8) + scale_fill_distiller(palette = "RdBu")+scale_x_continuous(expand = expansion(mult = .1))
p2 <- p+theme(legend.title = element_text(size = 20),legend.text =element_text(size = 20) )                                                                                                  
p2
ggsave("veen1.pdf",p2,dpi=600,units = 'in',width = 15,height = 15)
p <- ggVennDiagram(x1_down,edge_size = 1,label = 'count',edge_lty = "dashed",set_size = 6,label_size = 5) + scale_fill_distiller(palette = "RdBu")+scale_x_continuous(expand = expansion(mult = .1))
p1 <- p+theme(legend.title = element_text(size = 15),legend.text =element_text(size = 15) )  
p1
ggsave("veen2.pdf",p1,dpi=600)
