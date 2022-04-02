---
title: WGCNA
date: 2022-03-30 10:46:54
index_img: /img/download1.jpg 
categories:
- bioinformation
tags:
---
参考
https://www.zhouxiaozhao.cn/2020/08/03/WGCNA1/
https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
https://cloud.tencent.com/developer/article/1516749
### 1.定义
加权基因共表达网络分析 (WGCNA, Weighted correlation network analysis)是用来描述不同样品之间基因关联模式的系统生物学方法，可以用来鉴定高度协同变化的基因集, 并根据基因集的内连性和基因集与表型之间的关联鉴定候补生物标记基因或治疗靶点。相比于只关注差异表达的基因，WGCNA利用数千或近万个变化最大的基因或全部基因的信息识别感兴趣的基因集，并与表型进行显著性关联分析。一是充分利用了信息，二是把数千个基因与表型的关联转换为数个基因集与表型的关联，免去了多重假设检验校正的问题。
共表达网络：
点代表基因，边代表基因表达相关性。加权是指对相关性值进行冥次运算(冥次的值也就是软阈值 (power, pickSoftThreshold这个函数所做的就是确定合适的power))。无向网络的边属性计算方式为abs(cor(genex, geney)) ^ power；有向网络的边属性计算方式为(1+cor(genex, geney)/2) ^ power; sign hybrid的边属性计算方式为cor(genex, geney)^power if cor>0 else 0。这种处理方式强化了强相关，弱化了弱相关或负相关，使得相关性数值更符合无标度网络特征，更具有生物意义。如果没有合适的power，一般是由于部分样品与其它样品因为某种原因差别太大导致的，可根据具体问题移除部分样品或查看后面的经验值。

Module(模块）
高度內连的基因集。在无向网络中，模块内是高度相关的基因。在有向网络中，模块内是高度正相关的基因。把基因聚类成模块后，可以对每个模块进行三个层次的分析：1. 功能富集分析查看其功能特征是否与研究目的相符；2. 模块与性状进行关联分析，找出与关注性状相关度最高的模块；3. 模块与样本进行关联分析，找到样品特异高表达的模块。

Hub gene: 关键基因 (连接度最多或连接多个模块的基因)。
Adjacency matrix (邻接矩阵)：
基因和基因之间的加权相关性值构成的矩阵。
TOM (Topological overlap matrix)：
把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得的新距离矩阵，这个信息可拿来构建网络或绘制TOM图。
模块内部基因连接度分析
● Intramodular connectivity KIM-模块内部连接度IC:某个模块中的基因与该模块中其他基因的关联程度（共表达程度）。可用来衡量模块身份（module membership,MM).
● Module Membership MM,or Epigengene-based connectivity KME:模块身份，用一个基因在所有样本中的表达语与某个模块特征值的表达谱的相关性，来衡量这个基因在这个模块中的身份。
● KME值接近0,说明这个基因不是该模块的成员：KME接近1或者－1,说明这个基因与该模块密切相关（正相关或者负相关）。 可以对所有基因计算相对某个模块的KME值，并不一定要是该模块的成员。 KME与KIM高度相关。某个模块中KIM值高的hub基因一定与该模块的KME也很高。 KME与KIM的区别：IC衡量基因在特定模块中的身份，MM衡量基因在全局网络中的位置。
● 筛选关键基因： TOM值（模块调控系表中的weight值）大于阈值（默认是0.15)的两个基因才认为是相关的，然后计算每个基因的连接度。即先筛选有足够强度的关系，然后计算连接度。
● 模块内部高连接度的基因，模块内排名前30或者10%(KME或KIM). 筛选关键基因：将该基因模块身份MM相对于基因显著性GS做散点图，选择右上角MM和GS均高的基因进一步分析。 基因显著性值（Gene significance,GS)因变量水平的相关系数。衡量基因与表型性状的关联程度，GS越高，说明与表型越相关，越具有生物学意义。GS可以为正值或负值（正相关或负相关） Cytoscape中一般用weight值（TOM值）来绘制网络图。
过程
1. 构建基因共表达网络：使用加权的表达相关性。
2. 识别基因集：基于加权相关性，进行层级聚类分析，并根据设定标准切分聚类结果，获得不同的基因模块，用聚类树的分枝和不同颜色表示。
3. 如果有表型信息，计算基因模块与表型的相关性，鉴定性状相关的模块。
4. 研究模型之间的关系，从系统层面查看不同模型的互作网络。
5. 从关键模型中选择感兴趣的驱动基因，或根据模型中已知基因的功能推测未知基因的功能。
6. 导出TOM矩阵，绘制相关性图。
### 2.代码
~~~
#BiocManager::install("WGCNA",force = TRUE)
# BiocManager::install("GO.db")
library(WGCNA)
library(reshape2)
library(stringr)  
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()

data_type <- "Gene Expression Quantification"
data_catetory <- "Transcriptome Profiling"
expr <- GDCquery(project = "TCGA-OV",
                 data.category = data_catetory,
                 data.type = data_type,
                 workflow.type = "HTSeq - FPKM")
GDCdownload(expr,method = "api",files.per.chunk = 100)
EXPR2 <- GDCprepare(query = expr)
exprmat <- assay(EXPR2)
exprlog <- log2(exprmat+1)

# WGCNA_matrix = t(exprmat[order(apply(exprmat,1,mad), decreasing = T)[1:5000],])
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(exprlog,1,mad)
dataExprVar <- exprlog[which(m.mad > 
                               max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

## 转换为样品在行，基因在列的矩阵
datExpr <- as.data.frame(t(dataExprVar))
# datExpr0 <- WGCNA_matrix  ## top 5000 mad genes
# datExpr <- datExpr0 

gsg = goodSamplesGenes(datExpr, verbose = 3)


if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

## 查看是否有离群样品
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="sample",label=F)

abline(h = 170, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 170, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr0 = datExpr[keepSamples, ]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# 还差临床信息 ------------------------------------------------------------------
sampleName <- rownames(datExpr0)
query <-GDCquery(project = "TCGA-OV", 
                 data.category = "Clinical",
                 data.type = "Clinical Supplement", 
                 data.format = "BCR Biotab")

GDCdownload(query)
clinicalov <-GDCprepare(query)
clinical <- clinicalov$clinical_patient_ov
clinical <- 
  clinical[match(substr(sampleName,1,12),clinical$bcr_patient_barcode),c(2,6,33) ] %>% 
  as.data.frame()

rownames(clinical) <- sampleName

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
# 筛选标准。R-square=0.85
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
power = sft$powerEstimate
power
nGenes = ncol(datExpr0)

nSamples = nrow(datExpr0)



##21换成了0.8power改为5
project <- readRDS("~/hrdect/HRDzong/WGCNA/project.rds")
for( i in 1:nrow(project)){
  pdf(file=paste("/home/lhm/hrdect/HRDzong/WGCNA/sftpicture/",project[i,],".pdf",sep = ""),
      width = 9,height = 5)
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
  sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab=paste("Soft Threshold (power)","\n",project[i,]),ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
  abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab=paste("Soft Threshold (power)","\n",project[i,]),ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}



power = sft$powerEstimate
power
nGenes = ncol(datExpr0)

nSamples = nrow(datExpr0)

# This step is the bedrock of all network analyses using the WGCNA methodology. We present three different ways of
# constructing a network and identifying modules:
#   a. Using a convenient 1-step network construction and module detection function, suitable for users wishing to arrive
# at the result with minimum effort;
# b. Step-by-step network construction and module detection for users who would like to experiment with customized/alternate methods;
# c. An automatic block-wise network construction and module detection method for users who wish to analyze data
# sets too large to be analyzed all in one.


project <- readRDS("~/hrdect/HRDzong/WGCNA/project.rds")
for (i in 1:nrow(project)){
  i=12
  load(paste("~/hrdect/HRDzong/WGCNA/first_step/",project[i,],".Rdata",sep = ""))
  load(paste("~/hrdect/HRDzong/WGCNA/sft/",project[i,],"sft.Rdata",sep = ""))
  nGenes = ncol(datExpr0)
  cor <- WGCNA::cor
 if(i==21){
   print(T)
   net = blockwiseModules(datExpr0, power =5 ,maxBlockSize = nGenes,
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,loadTOM =TRUE,
                          saveTOMFileBase = "TOM",
                          verbose = 3)}
 else{
   print("power")
   net = blockwiseModules(datExpr0, power = sft$powerEstimate,maxBlockSize = nGenes,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,loadTOM =TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
 }
  save(net,file =paste("~/hrdect/HRDzong/WGCNA/netdata2/",project[i,],".Rdata",sep = ""))

}###修改了maxBlockSize

cor <- WGCNA::cor
net = blockwiseModules(datExpr0, power = 5,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,loadTOM =TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)



table(net$colors)

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#绘制模块之间相关性图
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
##绘制模块之间的相关性
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
# MEs = net$MEs;
geneTree = net$dendrograms[[1]];

nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)


datTraits <- clinical
datTraits$HH <- rep(c("1","0"),each=189)
datTraits$HY <- rep(c("0","1"),each=189)
datTraits <- datTraits[,-c(1,2)]



moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# Since we have a moderately large number of modules and traits, a suitable graphical representation will help in
# reading the table. We color code each association by the correlation value:
#   sizeGrWindow(10,6)
# Will display correlations and their p-values\



textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，
# 所有的连续型性状也可以跟基因的表达值算出相关系数。
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要
# # Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$HH);
names(weight) = "HH"
# names (colors) of the modules
# 首先计算模块与基因的相关性矩阵 
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
# 再计算性状与基因的相关性矩阵 
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")
# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析 

module = "blueviolet"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Luminal",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#网络的可视化


# 主要是可视化 TOM矩阵，WGCNA的标准配图
# 然后可视化不同 模块 的相关性 热图
# 不同模块的层次聚类图
# 还有模块诊断，主要是 intramodular connectivity
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOMsimilarityFromExpr(datExpr0, power = 5); 
plotTOM = dissTOM^7; 
diag(plotTOM) = NA; 
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
png("step7-Network-heatmap.png",width = 800,height = 600)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

# ###最后画模块和性状的关系  ---------------------------------------------------------
##性状
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr0, moduleColors)$eigengenes
## 只有连续型性状才能只有计算
## 这里把是否属 Luminal 表型这个变量0,1进行数值化
# design=model.matrix(~0+ datTraits$HH)
# colnames(design)=levels(datTraits$HH)
Luminal = as.data.frame(design[,2]);
names(Luminal) = "Luminal"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, Luminal))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);

par(cex = 0.9)
png("step7-Eigengene-dendrogram.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
## 模块的进化树
png("step7-Eigengene-dendrogram-hclust.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
## 性状与模块热

png("step7-Eigengene-adjacency-heatmap.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


# ##### -------------------------------------------------------------------
##step8:提取指定模块的基因名
## step 8 
# 主要是关心具体某个模块内部的基因
# Select module
module = "blueviolet";
# Select module probes
probes = colnames(datExpr0) 
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
head(modProbes)

# 如果使用WGCNA包自带的热图就很丑。
which.module="blueviolet";
dat=datExpr0[,moduleColors==which.module ] 
plotMat(t(scale(dat)),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
datExpr0[1:4,1:4]
dat=t(datExpr0[,moduleColors==which.module ] )
library(pheatmap)
pheatmap(dat ,show_colnames =F,show_rownames = F) 
n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
group_list=datTraits$HH
ac=data.frame(g=group_list)
rownames(ac)=colnames(n) 
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac )
# 可以很清晰的看到，所有的形状相关的模块基因
# 其实未必就不是差异表达基因。

# -- ----------------------------------------------------------------------

# Step9: 模块的导出
# 主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape,VisANT等网络可视化软件。

# -- ----------------------------------------------------------------------

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr0, power = 5); 
# Select module
module = "blueviolet";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 
首先是导出到VisANT

vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)
然后是导出到cytoscape

cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
);
#如果模块包含的基因太多，网络太复杂，还可以进行筛选，比如：

nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]
~~~

~~~
Error1
net = blockwiseModules(datExpr0, power = 5,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,loadTOM =TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
Error in (new("standardGeneric", .Data = function (x, y = NULL, use = "everything",  :  
  参数没有用(weights.x = NULL, weights.y = NULL, cosine = FALSE)
                                                   
 WGCNA包与其他函数冲突导致的；
输入
cor <- WGCNA::cor                                             
~~~