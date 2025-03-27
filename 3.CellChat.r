http://127.0.0.1:8000/notebooks/keti/singlecell_eqtl/code/Untitled2.ipynb?kernel_name=ir

library(CellChat)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(ggalluvial)
library(extrafont)
options(stringsAsFactors = FALSE)

plot_size <- function(x, y) {
    options(repr.plot.width = x, repr.plot.height = y)
}

sce1<-readRDS("/path/singlecell_eqtl/result/RDS/all_subtype_annotated.rds")
sce1

cellchat <- createCellChat(object=sce1,group.by = "Type")
#cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "cell_type")
cellchat


summary(cellchat)
str(cellchat)
levels(cellchat@idents)


#cellchat <- setIdent(cellchat, ident.use = "cell_type")
groupSize <- as.numeric(table(cellchat@idents))  
#查看每个cluster有多少个细胞，后面画图的时候需要用到这个值
groupSize

table(cellchat@idents)

CellChatDB <- CellChatDB.human
#导入小鼠是CellChatDB <- CellChatDB.mouse
str(CellChatDB) #查看数据库信息
#包含interaction、complex、cofactor和geneInfo这4个dataframe
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]


head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new() #下一步不出图的时候运行


showDatabaseCategory(CellChatDB)

#在CellChat中，我们还可以先择特定的信息描述细胞间的相互作用，可以理解为从特定的侧面来刻画细胞间相互作用，比用一个大的配体库又精细了许多。

unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
#选择"Secreted Signaling"进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
CellChatDB.use


# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object

#对表达数据进行预处理，用于细胞间的通信分析。首先在一个细胞组中识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。
#如果配体或受体过表达，则识别过表达配体和受体之间的相互作用。

## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个）
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)#4个进程
cellchat

#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
cellchat

#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human) 
cellchat
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project

#通过为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信。

#推断配体-受体水平细胞通讯网络（结果储存在@net下面，有一个概率值和对应的pval）
#⚠️这一步也是CellChat比CellPhoneDB多的一步
#通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的通信概率。

 options(future.globals.maxSize=4000000000)

#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
cellchat

# 如果某些细胞组中只有少量细胞，则过滤掉细胞间通信
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
head(df.net)

write.csv(df.net, "/path/1net_lr.csv")

#推断信号通路水平的细胞通讯网络（结果储存在@netP下面，有一个概率值和对应的pval）
#我们可以通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
df.netp

a<-subsetCommunication(cellchat)
a

length(unique(df.netp$pathway_name))
max(df.netp$pval)

write.csv(df.netp, "/path/2net_pathway.csv")
#至此，统计分析部分已经完成。#设置slot.name = "netP" 保存显著的结果

#所有细胞群总体观：细胞互作数量与强度统计分析

#统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
cellchat <- aggregateNet(cellchat)
cellchat

plot_size(5,6)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
groupSize

cellchat@net$count

13+	9	+14+	17+	40	+15	+18	+17
13+11+19+15+35+19+15+19

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                     label.edge= F, title.name = "Number of interactions")


netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                     label.edge= F,title.name = "Interaction weights/strength")
# save as TIL/net_number_strength.pdf

#左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。发出箭头的细胞表达配体，箭头指向的细胞表达受体。配体-受体对越多，线越粗。右图：互作的概率/强度值（强度就是概率值相加）

plot_size(18,15)
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                    arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    
}
#每个细胞如何跟别的细胞互作（number of interaction图）

#vertex.weight指定节点权重;weight.scale设置是否对权重进行缩放;arrow.width和arrow.size控制箭头的宽度和大小;
#edge.weight.max设定边的最大权重值;title.name为显示的图标题
# save as TIL/net_number_individual.pdf

## 运行上述代码出现报错 Error in plot.new() : figure margins too large
# par("mar")
## [1] 5.1 4.1 4.1 2.1
# par(mar=c(1,1,1,1))
# 重新运行上面的代码

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as TIL/net_strength_individual.pdf

#⚠️注：层次图、网络图、和弦图、热图只是不同的展示方法，展示的内容和代表的意思一模一样
#比如在前面的功能富集分析或case control的比较中找到了一些信号通路差异，就可以进一步在细胞水平上验证。

cellchat@netP$pathways  #查看都有哪些信号通路
pathways.show <- c("VEGF") 
# 选择其中一个信号通路，比如说TGFb
length(cellchat@netP$pathways )


 pathways.show


plot_size(8,5)
levels(cellchat@idents)    # show all celltype

vertex.receiver = c(1,2) # define a numeric vector （淋系细胞）giving the index of the celltype as targets
#par(mar=c(5.1,4.1,4.1,2.1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout="hierarchy",vertex.weight = groupSize)
# save as TIL/CXCL_hierarchy.pdf



groupSize

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# save as TIL/CXCL_circle.pdf

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# save as TIL/CXCL_chord.pdf

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# 纵轴是发出信号的细胞，横轴是接收信号的细胞，热图颜色深浅代表信号强度。上侧和右侧的柱子是纵轴和横轴强度的累积

#计算配体受体对选定信号通路的贡献值（在这里就是查看哪条信号通路对测试通路VEGF贡献最大）
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.VEGF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE) #提取对VEGF有贡献的所有配体受体 
# save as TIL/CXCL_LR_contribution.pdf
pairLR.VEGF



#提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.VEGF[1,] 
vertex.receiver = c(1,2,4,6) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# save as TIL/CXCL_hierarchy2.pdf


netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# save as TIL/CXCL_circle2.pdf



netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
# save as TIL/CXCL_chord2.pdf



# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Access all the signaling pathways showing significant communications将所有信号通路找出来
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(1,2,4,6) #不画层次图就不需要这一步
dir.create("all_pathways_com_circle") #创建文件夹保存批量画图结果
setwd("all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
      # Visualize communication network associated with both signaling pathway and individual L-R pairs
      netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
                vertex.receiver = vertex.receiver, layout = "circle") #绘制网络图
      # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
      gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
      ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
             plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
    }
    setwd("../")


plot_size(10,12)
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = c(5), 
                     targets.use = c(1,2,3,4,6,7,8), remove.isolate = FALSE)
#ggsave("Mye_Lymph_bubble.pdf", p, width = 8, height = 12) #髓系对淋巴的调节
# save as TIL/Mye_Lymph_bubble.pdf
p


plot_size(10,5)
#比如制定CCL和CXCL这两个信号通路
netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(1,2,3,4,6,7,8), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) based on user's input
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","TGFb"))
netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(1,2,3,4,6,7,8), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)
#CCL，CXCL和TGFb信号通路参与单核细胞和树突状细胞对T细胞的调控作用情况

## Plot the signaling gene expression distribution
p1 = plotGeneExpression(cellchat, signaling = "TGFb")
p1
#ggsave("TGFb_GeneExpression_vln.pdf", p, width = 8, height = 8)
p2 = plotGeneExpression(cellchat, signaling = "TGFb", type = "dot")
p2
#ggsave("TGFb_GeneExpression_dot.pdf", p, width = 8, height = 6)

#通讯网络系统分析使用了三种算法：社会网络分析、NMF分析和流行学习与分类
#⚠️：不同的算法算出来的结果可能会相互矛盾，需要结合生物学知识加以判断

#计算网络中心性权重

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#通过计算每个细胞群的网络中心性指标，识别每类细胞在信号通路中的角色/作用


netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                   width = 15, height = 6, font.size = 10)
# # save as TIL/SNA_CXCL_signalingRole.pdf


plot_size(15,8)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5)
ht1 + ht2
# 横轴是细胞类型，纵轴是pathway。左图是各个细胞类型中各个通路发出信号的强度，右图是各个细胞类型中各个通路接受信号的强度

#计算分解成几个因子(pattern)比较合适（这一步运行比较慢 。在使用NMF对细胞进行亚群细分时，如果不测试的话，最好选择比细胞类型多一点的值）
selectK(cellchat, pattern = "outgoing")
# save as TIL/NMF_outgoing_selectK.pdf



nPatterns = 5 # 挑选曲线中第一个出现下降的点（从5就开始下降了）
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, 
                                          width = 5, height = 9, font.size = 6)
# save as TIL/NMF_outgoing_comPattern_heatmap.pdf
#按selectK算出来的pattern值分为了2个大的pattern，左图纵轴是细胞类型，右图纵轴是信号通路

cellchat

cellchat

netAnalysis_river(cellchat, pattern = "outgoing")
# save as TIL/NMF_outgoing_comPattern_river.pdf  不知道为什么图片没显示

netAnalysis_dot(cellchat, pattern = "outgoing")
# save as TIL/NMF_outgoing_comPattern_dotplot.pdf

selectK(cellchat, pattern = "incoming") 
# save as TIL/NMF_incoming_selectK.pdf

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, 
                                          width = 5, height = 9, font.size = 6)
# save as TIL/NMF_incoming_comPattern_heatmap.pdf

incoming

library(ggalluvial)
netAnalysis_river(cellchat, pattern = "incoming")
# save as TIL/NMF_incoming_comPattern_river.pdf

netAnalysis_dot(cellchat, pattern = "incoming")
# save as TIL/NMF_incoming_comPattern_dotplot.pdf


#把共同起作用的信号通路归纳在一起，分为基于功能的归纳和基于拓扑结构的归纳


cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat
cellchat <- netClustering(cellchat, type = "functional")
cellchat

#Error in do_one(nmeth) : NA/NaN/Inf in foreign function call (arg 1)
p = netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
p
#ggsave("Manifold_functional_cluster.pdf", p, width = 8, height = 6)
#netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
p = netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
p
#ggsave("Manifold_structural_cluster.pdf", p, width = 8, height = 6)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

## The end
saveRDS(cellchat, file = "/path/cellchat_single_sample.rds",compress = F)


