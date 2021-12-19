rm(list = ls())
options(stringsAsFactors = F)
load("forWGCNA.Rdata")
library(edgeR)
exp <- log2(cpm(exp1)+1)
exp <- t(exp)
library(WGCNA)
gsg = goodSamplesGenes(exp, verbose = 3)
gsg$allOK
if(T){
  sampleTree = hclust(dist(exp), method = "average")
  pdf(file = "sampleClustering.pdf", width = 12, height = 9)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  
  abline(h = 200, col = "red") 
  dev.off()
}


clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust) 
keepSamples = (clust==1)
datExpr = exp[keepSamples, ]
datExpr=as.data.frame(datExpr)
allTraits=as.data.frame(clinical_data1[,35])
rownames(allTraits) <- rownames(clinical_data1)
colnames(allTraits) <- "group"
allTraits$name=rownames(allTraits)
traitRows = match(rownames(datExpr), rownames(allTraits))
allTraits=allTraits[traitRows,]

allTraits <- allTraits[,1]
allTraits <- as.data.frame(allTraits)
row.names(allTraits) <- rownames(datExpr)         
colnames(allTraits) <- "group"
allTraits$group <- c(rep(0,16),rep(1,30))
identical(rownames(allTraits),rownames(datExpr))

datTraits=allTraits
sampleTree2 = hclust(dist(datExpr), method = "average")

traitColors = numbers2colors(datTraits, signed = FALSE)

if(T){
  pdf(file = "Sample dendrogram and trait heatmap.pdf", width = 18, height = 10)
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
}

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

if(T){
  pdf(file = "Soft threshold.pdf", width = 18, height = 10)
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  
  abline(h=0.92,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}

net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 60,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "GDM_TOM",
                       verbose = 3)

table(net$colors)

mergedColors = labels2colors(net$colors)

if(T){
  pdf(file = "DendroAndColors.pdf", width = 18, height = 10)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]


nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p") 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


if(T){
  pdf(file = "labeledHeatmap_new.pdf", width = 18, height = 10)
  
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 15, 3, 3))
  
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = "GDM",
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 1.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"),
                 cex.lab = 2,
                 #cex.legendLabel = 10,
                 cex.main=3)
  dev.off()
}


modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


weight = as.data.frame(datTraits$group)
names(weight) = "GDM"

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"))

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

if(T){
  module = "greenyellow"
  pdf(file = paste0(module,"-MM-GS-scatterplotnew.pdf"), width = 16, height = 10)
  column = match(module, modNames) 
  moduleGenes = moduleColors==module 
  par(mfrow = c(1.6
                ,1))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for GDM",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 2, cex.lab = 1.5, cex.axis = 2, col = module)
  dev.off()
}

if(T){
  nSelect = 500
  set.seed(10)
  
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6)
  select = sample(nGenes, size = nSelect)
  selectTOM = dissTOM[select, select]
  
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select]
  
  pdf(file = paste0("Sub500-netheatmap1.pdf"), width = 5, height = 5)
  plotDiss = selectTOM^7
  diag(plotDiss) = NA 
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),cex.main=20)
  dev.off()
}


MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

weight = as.data.frame(datTraits$group)
names(weight) = "GDM"

MET = orderMEs(cbind(MEs, weight))


if(T){
  pdf(file = "Eigengene-dengro-heatmap1.pdf", width = 10, height = 15)
  par(cex = 2)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(4,4,1,2), cex.lab = 1, xLabelsAngle
                        = 90)
  dev.off()
}


TOM = TOMsimilarityFromExpr(datExpr, power = 6)

module = "greenyellow"

datExpr=as.data.frame(datExpr)
probes = names(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule]

modProbes_greenyellow=modProbes


modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

library(clusterProfiler)
library(org.Hs.eg.db)
gene_symbol <- bitr(modProbes_greenyellow,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
save(gene_symbol,file = "WGCNA_sym.Rdata")