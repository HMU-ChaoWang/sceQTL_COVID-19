monocle2官网：http://cole-trapnell-lab.github.io/monocle-release/docs/

#packageVersion("monocle")

library(Seurat)
library(monocle)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggsci)
library(data.table)
library(tidyverse)
library(harmony)
#library(DoubletFinder)
library(patchwork)
library(ggpubr)
library(reshape2)
library(ggpubr)
#library(progeny)
library(gplots)
library(SingleR)
library(ggpubr)
library(ggsci)
library(data.table)
library(scRNAtoolVis)
library(COSG)
library(RColorBrewer) 
library(ggsci)
library(colorspace)
library(viridis)
library(scales)
library(MatrixEQTL) 

plot_size <- function(x, y) {
    options(repr.plot.width = x, repr.plot.height = y)
}

#导入注释好的seurat对象（已注释）
Monocyte<-readRDS("/path/Monocyte_annotated.rds")

Monocyte@meta.data

table(Monocyte@meta.data$subtype)

Idents(Monocyte) <- Monocyte$subtype
classical <- subset(Monocyte,idents='Monocyte classical')

classical

nonclassical <- subset(Monocyte,idents='Monocyte nonclassical')

nonclassical

Monocyte@meta.data

####筛选部分细胞#####
sample_seob <- function(obj,group.by="subtype",sp.size=NULL,diet="true",sp.total=1000) {
all <- obj
if (diet=="true") {
all <- DietSeurat(all)
}
if (is.null(sp.size)) {
nlen <- length(unique(all@meta.data[,group.by]))
sp.size <- ceiling(sp.total/nlen)
}
seob_list <- list()
i <- 1
for (sc in unique(all@meta.data[,group.by])){
cellist <- colnames(all)[which(all@meta.data[,group.by] == sc)]
ob <- subset(all, cells=cellist)
if (length(colnames(ob)) > sp.size) {
ob <- subset(ob,cells=sample(colnames(ob), sp.size))
}
seob_list[[i]] <- ob
i <- i+1
}
all <- Reduce(merge,seob_list)
return(all)
}

Monocyte<- sample_seob(Monocyte,group.by="subtype",diet="true",sp.total=30000)


Monocyte

table(Monocyte@meta.data$subtype)

Monocyte<-SCTransform(Monocyte)
Monocyte <- FindVariableFeatures(Monocyte)
Monocyte

result_path<-"/path/singlecell_eqtl/result/monocle/monocyte/"
setwd("/path/singlecell_eqtl/result/monocle/monocyte")

##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
data <- as(as.matrix(Monocyte@assays$RNA@counts), 'sparseMatrix')

p_data <- Monocyte@meta.data

##提取表型信息到p_data(phenotype_data)里面
pd <- new('AnnotatedDataFrame', data = p_data)

#p_data$celltype <- pbmc@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加

##提取基因信息 如生物类型、gc含量等
fData <- data.frame(gene_short_name = row.names(Monocyte@assays$RNA), row.names = row.names(Monocyte@assays$RNA))
##data的行数与f_data的行数相同(gene number), data的列数与p_data的行数相同(cell number)

fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
#negbinomial.size()和negbinomial()：输入的表达矩阵为UMI,一般适用于10x的数据；negbinomial()的结果更准确，但是计算更耗时；一般建议采用negbinomial.size()。
#稀疏矩阵用negbinomial.size()
#tobit():适用于输入的表达矩阵为FPKM或者TPM, 构建monocle2的class时会自动进行log化计算
#gaussianff():输入为log化后的FPKM或者TPM。(目前在单细胞数据中，FPKM已不多用，smart-seq2平台数据一般采用TPM)

#monocle_cds <- newCellDataSet(as(as.matrix(data), "sparseMatrix"), #创建一个CDS对象
                       phenoData = pd, 
                       featureData = fd,
                       expressionFamily = tobit())

##读取数据
data <- fread("fpkm.txt",data.table = F,header = T)
pd <-  fread("metadata.txt",data.table = F,header = T)
fd <- fread("gene_annotations.txt",data.table = F,header = T)
##创建
pd <- new("AnnotatedDataFrame", data = pd)   #数据框的另一种形式
fd <- new("AnnotatedDataFrame", data = fd)
HSMM <- newCellDataSet(as.matrix(data),
                       phenoData = pd, featureData = fd,
                       expressionFamily = tobit())
###如果数据量大，建议转化为稀疏矩阵
HSMM <- newCellDataSet(as(as.matrix(data), "sparseMatrix"), #创建一个CDS对象
                       phenoData = pd, 
                       featureData = fd,
                       expressionFamily = tobit())
importCDS(pbmc)



monocle_cds <- estimateSizeFactors(monocle_cds)#size facotr帮助我们标准化细胞之间的mRNA的差异。

monocle_cds <- estimateDispersions(monocle_cds)
#离散度值可以帮助我们进行后续的差异分析 （类似于seurat的数据归一化处理）
#与seurat把标准化后的表达矩阵保存在对象中不同，monocle只保存一些中间结果在对象中，需要用时再用这些中间结果转化。
#经过上面三个函数的计算，mycds对象中多了SizeFactors、Dipersions、num_cells_expressed和num_genes_expressed等信息。

View(pData(monocle_cds))

#因为Seurat已经完成细胞过滤，此步可省略
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.5) #这一操作会在fData(cds)中添加一列num_cells_expressed
print(head(fData(monocle_cds)))#此时有13714个基因
expressed_genes <- row.names(subset(fData(monocle_cds),
    num_cells_expressed >= 5)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。


subset(fData(monocle_cds),
    num_cells_expressed >= 5)

View(fData(monocle_cds))

Monocle官网教程提供了4个分类方法：
Classifying cells by type
Clustering cells without marker genes
Clustering cells using marker genes
Imputing cell type
建议先将细胞注释好再进行Monocle分析，不建议使用monocle做细胞分类。

The ordering workflow
Step 1: choosing genes that define progress
Step 2: reducing the dimensionality of the data
Step 3: ordering the cells in pseudotime
---------------------------------------
Monocle官网教程提供了4个选择方法：
选择发育差异表达基因
选择clusters差异表达基因
选择离散程度高的基因
自定义发育marker基因

理想情况下，我们希望尽可能少地使用正在研究的系统生物学的先验知识。
我们希望从数据中发现重要的排序基因，而不是依赖于文献和教科书，因为这可能会在排序中引入偏见。
我们将从一种更简单的方法开始，但是我们通常推荐一种更复杂的方法，称为“dpFeature”。

disp_table <- dispersionTable(monocle_cds) 
head(disp_table)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
monocle_cds <- setOrderingFilter(monocle_cds, disp.genes)
plot_ordering_genes(monocle_cds)

express_genes <- VariableFeatures(Monocyte)
dim(express_genes)

express_genes <- VariableFeatures(Monocyte)
monocle_cds <- setOrderingFilter(monocle_cds, express_genes)
plot_ordering_genes(monocle_cds)


deg.cluster <- FindAllMarkers(Monocyte)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
monocle_cds <- setOrderingFilter(monocle_cds, express_genes)
plot_ordering_genes(monocle_cds)



#这一步输入的expressed_genes来自于步骤4。
#后续分析使用的是该方法
#也可输入seurat筛选出的高变基因：
express_genes <- VariableFeatures(Monocyte) 
head(express_genes)

monocle_cds <- setOrderingFilter(monocle_cds,express_genes)

monocle_cds <- reduceDimension(monocle_cds,max_components = 2,
               method = 'DRRTree')

monocle_cds <- orderCells(monocle_cds,reverse = T)#reverse = T,颠倒顺序

plot_size(15,8)
epi_trajectory_subtype <- plot_cell_trajectory(monocle_cds,color_by="State", 
                     size=10,show_backbone=TRUE)+
    theme(
        legend.position = 'right',
        legend.title =element_text(size=20),
        legend.text = element_text(size=16)
    )
epi_trajectory_subtype

plot_size(8,8)
plot_cell_trajectory(monocle_cds,color_by="subtype", 
                     size=1,show_backbone=TRUE)

plot_size(8,8)
plot_cell_trajectory(monocle_cds,color_by="Pseudotime", 
                     size=1,show_backbone=TRUE)
#以pseudotime值上色 (Pseudotime是monocle2基于细胞基因表达信息计算的概率，表示时间的先后。

View(pData(monocle_cds))

#saveRDS(Monocyte,'/path/Healthy_sct_reduction.rds',compress = F)

#如果跑monocle的话就不用跑6

Monocyte@meta.data=pData(monocle_cds)
Monocyte@meta.data

Idents(Monocyte) <- Monocyte@meta.data$State
type <- unique(Monocyte$State)
for (i in type) {
    print(i)
    celltype_obj <- subset(Monocyte, idents = i)
    vag_exp <- AverageExpression(celltype_obj, assays = "SCT", group.by = "orig.ident",
        slot = "data",add.ident=NULL)
    write.csv(vag_exp, paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/", i,
        ".csv"), row.names = T, quote = F)
}

#eqtls<-list.files("E:\\陈欣雨\\未完成课题\\sc-eqtl\\data\\need\\type_exp\\",pattern="*.csv")
type<-c("1", "2", "3", "4","5")
setwd("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp")

for(i in 1:length(type))

{
  eqtl<-read.csv(file=paste0(type[i],".csv"),header=TRUE,sep=",")
 write.table(eqtl,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/",
           type[i],".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  print(type[i])

 }


genotype_transpose<-fread("/path/singlecell_eqtl/eqtl/22+X/need/1genotype_transpose.txt",header=F,data.table=F)
colnames(genotype_transpose)<-c(
         
 "ID","SCT.GSM2025728","SCT.GSM2025729","SCT.GSM2025730","SCT.GSM2025731",
"SCT.GSM2025732","SCT.GSM2025733","SCT.GSM2025734","SCT.GSM2025735",
"SCT.GSM2025736","SCT.GSM2025737","SCT.GSM2025738","SCT.GSM2025739",
"SCT.GSM2025740","SCT.GSM2025741","SCT.GSM2025742","SCT.GSM2025743",
"SCT.GSM2025744","SCT.GSM2025745","SCT.GSM2025746","SCT.GSM2025747",
"SCT.GSM2025748","SCT.GSM2025749","SCT.GSM2025750","SCT.GSM2025751",
"SCT.GSM2025752","SCT.GSM2025753","SCT.GSM2025754","SCT.GSM2025755",
"SCT.GSM2025756","SCT.GSM2025757","SCT.GSM2025758","SCT.GSM2025759",
"SCT.GSM2025760","SCT.GSM2025761","SCT.GSM2025762","SCT.GSM2025763",
"SCT.GSM2025764","SCT.GSM2025765","SCT.GSM2025766","SCT.GSM2025767",
"SCT.GSM2025768","SCT.GSM2025769","SCT.GSM2025770","SCT.GSM2025771",
"SCT.GSM2025772","SCT.GSM2025773","SCT.GSM2025774","SCT.GSM2025775",
"SCT.GSM2025776","SCT.GSM2025777","SCT.GSM2025778","SCT.GSM2025779",
"SCT.GSM2025780"
)
rownames(genotype_transpose)<-genotype_transpose[,1]
genotype_transpose<-genotype_transpose[,-1]
#genotype_transpose1<-genotype_transpose[1:20,]
#type<-c("CD4 T","Monocyte","cDC","NK","Platelet","CD8 T","pDC","B")

#subtype<-c("CD4 T memory","Monocyte classical","CD8 T effector","B naive","Monocyte nonclassical","CD4 T naive",
#     "cDC","NK","Platelet","CD8 T naive","pDC","B memory","Plasma")
#subtype_metadata<-read.txt(file=paste0(subtype[i],".csv"),header=TRUE,sep=",")
subtype_metadata<-fread("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/4.txt",header=T,sep="\t")

all_genotype_cellsubtype<-matrix(nrow = dim(genotype_transpose)[1],ncol = dim(subtype_metadata)[2]-1)
  colnames(all_genotype_cellsubtype)<- colnames(subtype_metadata)[2:dim(subtype_metadata)[2]]
  
 rownames(all_genotype_cellsubtype)<-rownames(genotype_transpose)

  for(i in 1:dim(all_genotype_cellsubtype)[2]){  
      loc<-which(colnames(genotype_transpose)%in%colnames(all_genotype_cellsubtype)[i])
     all_genotype_cellsubtype[,i]<-genotype_transpose[,loc]
  }
write.table(all_genotype_cellsubtype,"/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_4.txt",sep=" ",quote=F)
head(all_genotype_cellsubtype,2)

genotype_transpose<-fread("/path/singlecell_eqtl/eqtl/22+X/need/1genotype_transpose.txt",header=F,data.table=F)
colnames(genotype_transpose)<-c(
         
 "ID","SCT.GSM2025728","SCT.GSM2025729","SCT.GSM2025730","SCT.GSM2025731",
"SCT.GSM2025732","SCT.GSM2025733","SCT.GSM2025734","SCT.GSM2025735",
"SCT.GSM2025736","SCT.GSM2025737","SCT.GSM2025738","SCT.GSM2025739",
"SCT.GSM2025740","SCT.GSM2025741","SCT.GSM2025742","SCT.GSM2025743",
"SCT.GSM2025744","SCT.GSM2025745","SCT.GSM2025746","SCT.GSM2025747",
"SCT.GSM2025748","SCT.GSM2025749","SCT.GSM2025750","SCT.GSM2025751",
"SCT.GSM2025752","SCT.GSM2025753","SCT.GSM2025754","SCT.GSM2025755",
"SCT.GSM2025756","SCT.GSM2025757","SCT.GSM2025758","SCT.GSM2025759",
"SCT.GSM2025760","SCT.GSM2025761","SCT.GSM2025762","SCT.GSM2025763",
"SCT.GSM2025764","SCT.GSM2025765","SCT.GSM2025766","SCT.GSM2025767",
"SCT.GSM2025768","SCT.GSM2025769","SCT.GSM2025770","SCT.GSM2025771",
"SCT.GSM2025772","SCT.GSM2025773","SCT.GSM2025774","SCT.GSM2025775",
"SCT.GSM2025776","SCT.GSM2025777","SCT.GSM2025778","SCT.GSM2025779",
"SCT.GSM2025780"
)
rownames(genotype_transpose)<-genotype_transpose[,1]
genotype_transpose<-genotype_transpose[,-1]
#genotype_transpose1<-genotype_transpose[1:20,]
#type<-c("CD4 T","Monocyte","cDC","NK","Platelet","CD8 T","pDC","B")

#subtype<-c("CD4 T memory","Monocyte classical","CD8 T effector","B naive","Monocyte nonclassical","CD4 T naive",
#     "cDC","NK","Platelet","CD8 T naive","pDC","B memory","Plasma")
#subtype_metadata<-read.txt(file=paste0(subtype[i],".csv"),header=TRUE,sep=",")
subtype_metadata<-fread("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/2.txt",header=T,sep="\t")

all_genotype_cellsubtype<-matrix(nrow = dim(genotype_transpose)[1],ncol = dim(subtype_metadata)[2]-1)
  colnames(all_genotype_cellsubtype)<- colnames(subtype_metadata)[2:dim(subtype_metadata)[2]]
  
 rownames(all_genotype_cellsubtype)<-rownames(genotype_transpose)

  for(i in 1:dim(all_genotype_cellsubtype)[2]){  
      loc<-which(colnames(genotype_transpose)%in%colnames(all_genotype_cellsubtype)[i])
     all_genotype_cellsubtype[,i]<-genotype_transpose[,loc]
  }
write.table(all_genotype_cellsubtype,"/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_2.txt",sep=" ",quote=F)
head(all_genotype_cellsubtype,2)

genotype_transpose<-fread("/path/singlecell_eqtl/eqtl/22+X/need/1genotype_transpose.txt",header=F,data.table=F)
colnames(genotype_transpose)<-c(
         
 "ID","SCT.GSM2025728","SCT.GSM2025729","SCT.GSM2025730","SCT.GSM2025731",
"SCT.GSM2025732","SCT.GSM2025733","SCT.GSM2025734","SCT.GSM2025735",
"SCT.GSM2025736","SCT.GSM2025737","SCT.GSM2025738","SCT.GSM2025739",
"SCT.GSM2025740","SCT.GSM2025741","SCT.GSM2025742","SCT.GSM2025743",
"SCT.GSM2025744","SCT.GSM2025745","SCT.GSM2025746","SCT.GSM2025747",
"SCT.GSM2025748","SCT.GSM2025749","SCT.GSM2025750","SCT.GSM2025751",
"SCT.GSM2025752","SCT.GSM2025753","SCT.GSM2025754","SCT.GSM2025755",
"SCT.GSM2025756","SCT.GSM2025757","SCT.GSM2025758","SCT.GSM2025759",
"SCT.GSM2025760","SCT.GSM2025761","SCT.GSM2025762","SCT.GSM2025763",
"SCT.GSM2025764","SCT.GSM2025765","SCT.GSM2025766","SCT.GSM2025767",
"SCT.GSM2025768","SCT.GSM2025769","SCT.GSM2025770","SCT.GSM2025771",
"SCT.GSM2025772","SCT.GSM2025773","SCT.GSM2025774","SCT.GSM2025775",
"SCT.GSM2025776","SCT.GSM2025777","SCT.GSM2025778","SCT.GSM2025779",
"SCT.GSM2025780"
)
rownames(genotype_transpose)<-genotype_transpose[,1]
genotype_transpose<-genotype_transpose[,-1]
#genotype_transpose1<-genotype_transpose[1:20,]
#type<-c("CD4 T","Monocyte","cDC","NK","Platelet","CD8 T","pDC","B")

#subtype<-c("CD4 T memory","Monocyte classical","CD8 T effector","B naive","Monocyte nonclassical","CD4 T naive",
#     "cDC","NK","Platelet","CD8 T naive","pDC","B memory","Plasma")
#subtype_metadata<-read.txt(file=paste0(subtype[i],".csv"),header=TRUE,sep=",")
subtype_metadata<-fread("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/3.txt",header=T,sep="\t")

all_genotype_cellsubtype<-matrix(nrow = dim(genotype_transpose)[1],ncol = dim(subtype_metadata)[2]-1)
  colnames(all_genotype_cellsubtype)<- colnames(subtype_metadata)[2:dim(subtype_metadata)[2]]
  
 rownames(all_genotype_cellsubtype)<-rownames(genotype_transpose)

  for(i in 1:dim(all_genotype_cellsubtype)[2]){  
      loc<-which(colnames(genotype_transpose)%in%colnames(all_genotype_cellsubtype)[i])
     all_genotype_cellsubtype[,i]<-genotype_transpose[,loc]
  }
write.table(all_genotype_cellsubtype,"/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_3.txt",sep=" ",quote=F)
head(all_genotype_cellsubtype,2)

genotype_transpose<-fread("/path/singlecell_eqtl/eqtl/22+X/need/1genotype_transpose.txt",header=F,data.table=F)
colnames(genotype_transpose)<-c(
         
 "ID","SCT.GSM2025728","SCT.GSM2025729","SCT.GSM2025730","SCT.GSM2025731",
"SCT.GSM2025732","SCT.GSM2025733","SCT.GSM2025734","SCT.GSM2025735",
"SCT.GSM2025736","SCT.GSM2025737","SCT.GSM2025738","SCT.GSM2025739",
"SCT.GSM2025740","SCT.GSM2025741","SCT.GSM2025742","SCT.GSM2025743",
"SCT.GSM2025744","SCT.GSM2025745","SCT.GSM2025746","SCT.GSM2025747",
"SCT.GSM2025748","SCT.GSM2025749","SCT.GSM2025750","SCT.GSM2025751",
"SCT.GSM2025752","SCT.GSM2025753","SCT.GSM2025754","SCT.GSM2025755",
"SCT.GSM2025756","SCT.GSM2025757","SCT.GSM2025758","SCT.GSM2025759",
"SCT.GSM2025760","SCT.GSM2025761","SCT.GSM2025762","SCT.GSM2025763",
"SCT.GSM2025764","SCT.GSM2025765","SCT.GSM2025766","SCT.GSM2025767",
"SCT.GSM2025768","SCT.GSM2025769","SCT.GSM2025770","SCT.GSM2025771",
"SCT.GSM2025772","SCT.GSM2025773","SCT.GSM2025774","SCT.GSM2025775",
"SCT.GSM2025776","SCT.GSM2025777","SCT.GSM2025778","SCT.GSM2025779",
"SCT.GSM2025780"
)
rownames(genotype_transpose)<-genotype_transpose[,1]
genotype_transpose<-genotype_transpose[,-1]
#genotype_transpose1<-genotype_transpose[1:20,]
#type<-c("CD4 T","Monocyte","cDC","NK","Platelet","CD8 T","pDC","B")

#subtype<-c("CD4 T memory","Monocyte classical","CD8 T effector","B naive","Monocyte nonclassical","CD4 T naive",
#     "cDC","NK","Platelet","CD8 T naive","pDC","B memory","Plasma")
#subtype_metadata<-read.txt(file=paste0(subtype[i],".csv"),header=TRUE,sep=",")
subtype_metadata<-fread("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/6.txt",header=T,sep="\t")

all_genotype_cellsubtype<-matrix(nrow = dim(genotype_transpose)[1],ncol = dim(subtype_metadata)[2]-1)
  colnames(all_genotype_cellsubtype)<- colnames(subtype_metadata)[2:dim(subtype_metadata)[2]]
  
 rownames(all_genotype_cellsubtype)<-rownames(genotype_transpose)

  for(i in 1:dim(all_genotype_cellsubtype)[2]){  
      loc<-which(colnames(genotype_transpose)%in%colnames(all_genotype_cellsubtype)[i])
     all_genotype_cellsubtype[,i]<-genotype_transpose[,loc]
  }
write.table(all_genotype_cellsubtype,"/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_6.txt",sep=" ",quote=F)
head(all_genotype_cellsubtype,2)

genotype_transpose<-fread("/path/singlecell_eqtl/eqtl/22+X/need/1genotype_transpose.txt",header=F,data.table=F)
colnames(genotype_transpose)<-c(
         
 "ID","SCT.GSM2025728","SCT.GSM2025729","SCT.GSM2025730","SCT.GSM2025731",
"SCT.GSM2025732","SCT.GSM2025733","SCT.GSM2025734","SCT.GSM2025735",
"SCT.GSM2025736","SCT.GSM2025737","SCT.GSM2025738","SCT.GSM2025739",
"SCT.GSM2025740","SCT.GSM2025741","SCT.GSM2025742","SCT.GSM2025743",
"SCT.GSM2025744","SCT.GSM2025745","SCT.GSM2025746","SCT.GSM2025747",
"SCT.GSM2025748","SCT.GSM2025749","SCT.GSM2025750","SCT.GSM2025751",
"SCT.GSM2025752","SCT.GSM2025753","SCT.GSM2025754","SCT.GSM2025755",
"SCT.GSM2025756","SCT.GSM2025757","SCT.GSM2025758","SCT.GSM2025759",
"SCT.GSM2025760","SCT.GSM2025761","SCT.GSM2025762","SCT.GSM2025763",
"SCT.GSM2025764","SCT.GSM2025765","SCT.GSM2025766","SCT.GSM2025767",
"SCT.GSM2025768","SCT.GSM2025769","SCT.GSM2025770","SCT.GSM2025771",
"SCT.GSM2025772","SCT.GSM2025773","SCT.GSM2025774","SCT.GSM2025775",
"SCT.GSM2025776","SCT.GSM2025777","SCT.GSM2025778","SCT.GSM2025779",
"SCT.GSM2025780"
)
rownames(genotype_transpose)<-genotype_transpose[,1]
genotype_transpose<-genotype_transpose[,-1]
#genotype_transpose1<-genotype_transpose[1:20,]
#type<-c("CD4 T","Monocyte","cDC","NK","Platelet","CD8 T","pDC","B")

#subtype<-c("CD4 T memory","Monocyte classical","CD8 T effector","B naive","Monocyte nonclassical","CD4 T naive",
#     "cDC","NK","Platelet","CD8 T naive","pDC","B memory","Plasma")
#subtype_metadata<-read.txt(file=paste0(subtype[i],".csv"),header=TRUE,sep=",")
subtype_metadata<-fread("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/8.txt",header=T,sep="\t")

all_genotype_cellsubtype<-matrix(nrow = dim(genotype_transpose)[1],ncol = dim(subtype_metadata)[2]-1)
  colnames(all_genotype_cellsubtype)<- colnames(subtype_metadata)[2:dim(subtype_metadata)[2]]
  
 rownames(all_genotype_cellsubtype)<-rownames(genotype_transpose)

  for(i in 1:dim(all_genotype_cellsubtype)[2]){  
      loc<-which(colnames(genotype_transpose)%in%colnames(all_genotype_cellsubtype)[i])
     all_genotype_cellsubtype[,i]<-genotype_transpose[,loc]
  }
write.table(all_genotype_cellsubtype,"/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_8.txt",sep=" ",quote=F)
head(all_genotype_cellsubtype,2)

genotype_transpose<-fread("/path/singlecell_eqtl/eqtl/22+X/need/1genotype_transpose.txt",header=F,data.table=F)
colnames(genotype_transpose)<-c(
         
 "ID","SCT.GSM2025728","SCT.GSM2025729","SCT.GSM2025730","SCT.GSM2025731",
"SCT.GSM2025732","SCT.GSM2025733","SCT.GSM2025734","SCT.GSM2025735",
"SCT.GSM2025736","SCT.GSM2025737","SCT.GSM2025738","SCT.GSM2025739",
"SCT.GSM2025740","SCT.GSM2025741","SCT.GSM2025742","SCT.GSM2025743",
"SCT.GSM2025744","SCT.GSM2025745","SCT.GSM2025746","SCT.GSM2025747",
"SCT.GSM2025748","SCT.GSM2025749","SCT.GSM2025750","SCT.GSM2025751",
"SCT.GSM2025752","SCT.GSM2025753","SCT.GSM2025754","SCT.GSM2025755",
"SCT.GSM2025756","SCT.GSM2025757","SCT.GSM2025758","SCT.GSM2025759",
"SCT.GSM2025760","SCT.GSM2025761","SCT.GSM2025762","SCT.GSM2025763",
"SCT.GSM2025764","SCT.GSM2025765","SCT.GSM2025766","SCT.GSM2025767",
"SCT.GSM2025768","SCT.GSM2025769","SCT.GSM2025770","SCT.GSM2025771",
"SCT.GSM2025772","SCT.GSM2025773","SCT.GSM2025774","SCT.GSM2025775",
"SCT.GSM2025776","SCT.GSM2025777","SCT.GSM2025778","SCT.GSM2025779",
"SCT.GSM2025780"
)
rownames(genotype_transpose)<-genotype_transpose[,1]
genotype_transpose<-genotype_transpose[,-1]
#genotype_transpose1<-genotype_transpose[1:20,]
#type<-c("CD4 T","Monocyte","cDC","NK","Platelet","CD8 T","pDC","B")

#subtype<-c("CD4 T memory","Monocyte classical","CD8 T effector","B naive","Monocyte nonclassical","CD4 T naive",
#     "cDC","NK","Platelet","CD8 T naive","pDC","B memory","Plasma")
#subtype_metadata<-read.txt(file=paste0(subtype[i],".csv"),header=TRUE,sep=",")
subtype_metadata<-fread("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/9.txt",header=T,sep="\t")

all_genotype_cellsubtype<-matrix(nrow = dim(genotype_transpose)[1],ncol = dim(subtype_metadata)[2]-1)
  colnames(all_genotype_cellsubtype)<- colnames(subtype_metadata)[2:dim(subtype_metadata)[2]]
  
 rownames(all_genotype_cellsubtype)<-rownames(genotype_transpose)

  for(i in 1:dim(all_genotype_cellsubtype)[2]){  
      loc<-which(colnames(genotype_transpose)%in%colnames(all_genotype_cellsubtype)[i])
     all_genotype_cellsubtype[,i]<-genotype_transpose[,loc]
  }
write.table(all_genotype_cellsubtype,"/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_9.txt",sep=" ",quote=F)
head(all_genotype_cellsubtype,2)



#type<-c("1", "2", "3", "4","5","6",  "7",  "8", "9")
type<-c("2","5") 

for(i in 1:length(type)){
base.dir = find.package("MatrixEQTL") 
useModel = modelLINEAR


SNP_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/1genotype_transpose.txt") 

snps_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-GSE165080-snpsloc.txt")
expression_file_name = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/",type[i],".txt" ,sep="")
                        
gene_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-18-genesloc.txt")
covariates_file_name =character()
    
    
output_file_name_cis = tempfile()
output_file_name_tra = tempfile() 

pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-20

cisDist = 1e6
errorCovariance = numeric()

snps = SlicedData$new() 
snps$fileDelimiter = " "      
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1         
snps$fileSkipColumns = 1      
snps$fileSliceSize = 2000    
snps$LoadFile(SNP_file_name)



gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile( expression_file_name )



snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, sep= "",fill = TRUE,na.strings = "NA")
genepos<-genepos[complete.cases(genepos),]




me = Matrix_eQTL_main(
  snps = snps,# 指定SNP 文件
  gene = gene,# 指定基因表达量文件
  #cvrt = cvrt,
  output_file_name     = output_file_name_tra, # tra指定输出文件
  pvOutputThreshold     = pvOutputThreshold_tra, # tra指定显著性P值
  useModel = useModel, # 指定使用的计算模型
  errorCovariance = errorCovariance, # 指定误差项的协方差矩阵
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis, # cis指定输出文件
  pvOutputThreshold.cis = pvOutputThreshold_cis, # cis指定显著性P值
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name_tra)
unlink(output_file_name_cis)


cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## 绘制近端远端P值得Q-Q图


plot(me)

dev.copy2pdf(file = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i] ,sep="",".pdf"), paper = "a4r")
print(type[i])
#dev.copy2pdf()
#导出文件
write.table(me$trans$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_trans_eqtls.txt",sep=""),sep = "\t",quote = F)
write.table(me$cis$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_cis_eqtls.txt",sep=""),sep = "\t",quote = F)

}

#type<-c("1", "2", "3", "4","5","6",  "7",  "8", "9")
type<-c("4") 

for(i in 1:length(type)){
base.dir = find.package("MatrixEQTL") 
useModel = modelLINEAR


SNP_file_name = ("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_4.txt") 

snps_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-GSE165080-snpsloc.txt")
expression_file_name = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/",type[i],".txt" ,sep="")
                        
gene_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-18-genesloc.txt")
covariates_file_name =character()
    
    
output_file_name_cis = tempfile()
output_file_name_tra = tempfile() 

pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-20

cisDist = 1e6
errorCovariance = numeric()

snps = SlicedData$new() 
snps$fileDelimiter = " "      
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1         
snps$fileSkipColumns = 1      
snps$fileSliceSize = 2000    
snps$LoadFile(SNP_file_name)



gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile( expression_file_name )



snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, sep= "",fill = TRUE,na.strings = "NA")
genepos<-genepos[complete.cases(genepos),]




me = Matrix_eQTL_main(
  snps = snps,# 指定SNP 文件
  gene = gene,# 指定基因表达量文件
  #cvrt = cvrt,
  output_file_name     = output_file_name_tra, # tra指定输出文件
  pvOutputThreshold     = pvOutputThreshold_tra, # tra指定显著性P值
  useModel = useModel, # 指定使用的计算模型
  errorCovariance = errorCovariance, # 指定误差项的协方差矩阵
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis, # cis指定输出文件
  pvOutputThreshold.cis = pvOutputThreshold_cis, # cis指定显著性P值
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name_tra)
unlink(output_file_name_cis)


cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## 绘制近端远端P值得Q-Q图


plot(me)

dev.copy2pdf(file = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i] ,sep="",".pdf"), paper = "a4r")
print(type[i])
#dev.copy2pdf()
#导出文件
write.table(me$trans$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_trans_eqtls.txt",sep=""),sep = "\t",quote = F)
write.table(me$cis$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_cis_eqtls.txt",sep=""),sep = "\t",quote = F)

}

#type<-c("1", "2", "3", "4","5","6",  "7",  "8", "9")
type<-c("2") 

for(i in 1:length(type)){
base.dir = find.package("MatrixEQTL") 
useModel = modelLINEAR


SNP_file_name = ("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_2.txt") 

snps_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-GSE165080-snpsloc.txt")
expression_file_name = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/",type[i],".txt" ,sep="")
                        
gene_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-18-genesloc.txt")
covariates_file_name =character()
    
    
output_file_name_cis = tempfile()
output_file_name_tra = tempfile() 

pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-20

cisDist = 1e6
errorCovariance = numeric()

snps = SlicedData$new() 
snps$fileDelimiter = " "      
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1         
snps$fileSkipColumns = 1      
snps$fileSliceSize = 2000    
snps$LoadFile(SNP_file_name)



gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile( expression_file_name )



snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, sep= "",fill = TRUE,na.strings = "NA")
genepos<-genepos[complete.cases(genepos),]




me = Matrix_eQTL_main(
  snps = snps,# 指定SNP 文件
  gene = gene,# 指定基因表达量文件
  #cvrt = cvrt,
  output_file_name     = output_file_name_tra, # tra指定输出文件
  pvOutputThreshold     = pvOutputThreshold_tra, # tra指定显著性P值
  useModel = useModel, # 指定使用的计算模型
  errorCovariance = errorCovariance, # 指定误差项的协方差矩阵
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis, # cis指定输出文件
  pvOutputThreshold.cis = pvOutputThreshold_cis, # cis指定显著性P值
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name_tra)
unlink(output_file_name_cis)


cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## 绘制近端远端P值得Q-Q图


plot(me)

dev.copy2pdf(file = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i] ,sep="",".pdf"), paper = "a4r")
print(type[i])
#dev.copy2pdf()
#导出文件
write.table(me$trans$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_trans_eqtls.txt",sep=""),sep = "\t",quote = F)
write.table(me$cis$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_cis_eqtls.txt",sep=""),sep = "\t",quote = F)

}

#type<-c("1", "2", "3", "4","5","6",  "7",  "8", "9")
type<-c("3") 

for(i in 1:length(type)){
base.dir = find.package("MatrixEQTL") 
useModel = modelLINEAR



    
SNP_file_name = ("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_3.txt") 

snps_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-GSE165080-snpsloc.txt")
expression_file_name = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/",type[i],".txt" ,sep="")
                        
gene_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-18-genesloc.txt")
    
    
covariates_file_name =character()
    
    
output_file_name_cis = tempfile()
output_file_name_tra = tempfile() 

pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-20

cisDist = 1e6
errorCovariance = numeric()

snps = SlicedData$new() 
snps$fileDelimiter = " "      
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1         
snps$fileSkipColumns = 1      
snps$fileSliceSize = 2000    
snps$LoadFile(SNP_file_name)



gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile( expression_file_name )



snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, sep= "",fill = TRUE,na.strings = "NA")
genepos<-genepos[complete.cases(genepos),]




me = Matrix_eQTL_main(
  snps = snps,# 指定SNP 文件
  gene = gene,# 指定基因表达量文件
  #cvrt = cvrt,
  output_file_name     = output_file_name_tra, # tra指定输出文件
  pvOutputThreshold     = pvOutputThreshold_tra, # tra指定显著性P值
  useModel = useModel, # 指定使用的计算模型
  errorCovariance = errorCovariance, # 指定误差项的协方差矩阵
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis, # cis指定输出文件
  pvOutputThreshold.cis = pvOutputThreshold_cis, # cis指定显著性P值
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name_tra)
unlink(output_file_name_cis)


cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## 绘制近端远端P值得Q-Q图


plot(me)

dev.copy2pdf(file = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i] ,sep="",".pdf"), paper = "a4r")
print(type[i])
#dev.copy2pdf()
#导出文件
write.table(me$trans$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_trans_eqtls.txt",sep=""),sep = "\t",quote = F)
write.table(me$cis$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_cis_eqtls.txt",sep=""),sep = "\t",quote = F)

}

#type<-c("1", "2", "3", "4","5","6",  "7",  "8", "9")
type<-c("6") 

for(i in 1:length(type)){
base.dir = find.package("MatrixEQTL") 
useModel = modelLINEAR


SNP_file_name = ("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_6.txt") 

snps_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-GSE165080-snpsloc.txt")
expression_file_name = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/",type[i],".txt" ,sep="")
                        
gene_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-18-genesloc.txt")
covariates_file_name =character()
    
    
output_file_name_cis = tempfile()
output_file_name_tra = tempfile() 

pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-20

cisDist = 1e6
errorCovariance = numeric()

snps = SlicedData$new() 
snps$fileDelimiter = " "      
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1         
snps$fileSkipColumns = 1      
snps$fileSliceSize = 2000    
snps$LoadFile(SNP_file_name)



gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile( expression_file_name )



snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, sep= "",fill = TRUE,na.strings = "NA")
genepos<-genepos[complete.cases(genepos),]




me = Matrix_eQTL_main(
  snps = snps,# 指定SNP 文件
  gene = gene,# 指定基因表达量文件
  #cvrt = cvrt,
  output_file_name     = output_file_name_tra, # tra指定输出文件
  pvOutputThreshold     = pvOutputThreshold_tra, # tra指定显著性P值
  useModel = useModel, # 指定使用的计算模型
  errorCovariance = errorCovariance, # 指定误差项的协方差矩阵
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis, # cis指定输出文件
  pvOutputThreshold.cis = pvOutputThreshold_cis, # cis指定显著性P值
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name_tra)
unlink(output_file_name_cis)


cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## 绘制近端远端P值得Q-Q图


plot(me)

dev.copy2pdf(file = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i] ,sep="",".pdf"), paper = "a4r")
print(type[i])
#dev.copy2pdf()
#导出文件
write.table(me$trans$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_trans_eqtls.txt",sep=""),sep = "\t",quote = F)
write.table(me$cis$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_cis_eqtls.txt",sep=""),sep = "\t",quote = F)

}

#type<-c("1", "2", "3", "4","5","6",  "7",  "8", "9")
type<-c("8") 

for(i in 1:length(type)){
base.dir = find.package("MatrixEQTL") 
useModel = modelLINEAR


SNP_file_name = ("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_8.txt") 

snps_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-GSE165080-snpsloc.txt")
expression_file_name = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/",type[i],".txt" ,sep="")
                        
gene_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-18-genesloc.txt")
covariates_file_name =character()
    
    
output_file_name_cis = tempfile()
output_file_name_tra = tempfile() 

pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-20

cisDist = 1e6
errorCovariance = numeric()

snps = SlicedData$new() 
snps$fileDelimiter = " "      
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1         
snps$fileSkipColumns = 1      
snps$fileSliceSize = 2000    
snps$LoadFile(SNP_file_name)



gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile( expression_file_name )



snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, sep= "",fill = TRUE,na.strings = "NA")
genepos<-genepos[complete.cases(genepos),]




me = Matrix_eQTL_main(
  snps = snps,# 指定SNP 文件
  gene = gene,# 指定基因表达量文件
  #cvrt = cvrt,
  output_file_name     = output_file_name_tra, # tra指定输出文件
  pvOutputThreshold     = pvOutputThreshold_tra, # tra指定显著性P值
  useModel = useModel, # 指定使用的计算模型
  errorCovariance = errorCovariance, # 指定误差项的协方差矩阵
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis, # cis指定输出文件
  pvOutputThreshold.cis = pvOutputThreshold_cis, # cis指定显著性P值
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name_tra)
unlink(output_file_name_cis)


cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## 绘制近端远端P值得Q-Q图


plot(me)

dev.copy2pdf(file = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i] ,sep="",".pdf"), paper = "a4r")
print(type[i])
#dev.copy2pdf()
#导出文件
write.table(me$trans$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_trans_eqtls.txt",sep=""),sep = "\t",quote = F)
write.table(me$cis$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_cis_eqtls.txt",sep=""),sep = "\t",quote = F)

}

#type<-c("1", "2", "3", "4","5","6",  "7",  "8", "9")
type<-c("9") 

for(i in 1:length(type)){
base.dir = find.package("MatrixEQTL") 
useModel = modelLINEAR


SNP_file_name = ("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/genotype_9.txt") 

snps_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-GSE165080-snpsloc.txt")
expression_file_name = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/",type[i],".txt" ,sep="")
                        
gene_location_file_name = ("/path/singlecell_eqtl/eqtl/22+X/need/COVID-19-18-genesloc.txt")
covariates_file_name =character()
    
    
output_file_name_cis = tempfile()
output_file_name_tra = tempfile() 

pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-20

cisDist = 1e6
errorCovariance = numeric()

snps = SlicedData$new() 
snps$fileDelimiter = " "      
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1         
snps$fileSkipColumns = 1      
snps$fileSliceSize = 2000    
snps$LoadFile(SNP_file_name)



gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile( expression_file_name )



snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, sep= "",fill = TRUE,na.strings = "NA")
genepos<-genepos[complete.cases(genepos),]




me = Matrix_eQTL_main(
  snps = snps,# 指定SNP 文件
  gene = gene,# 指定基因表达量文件
  #cvrt = cvrt,
  output_file_name     = output_file_name_tra, # tra指定输出文件
  pvOutputThreshold     = pvOutputThreshold_tra, # tra指定显著性P值
  useModel = useModel, # 指定使用的计算模型
  errorCovariance = errorCovariance, # 指定误差项的协方差矩阵
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis, # cis指定输出文件
  pvOutputThreshold.cis = pvOutputThreshold_cis, # cis指定显著性P值
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

unlink(output_file_name_tra)
unlink(output_file_name_cis)


cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## 绘制近端远端P值得Q-Q图


plot(me)

dev.copy2pdf(file = paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i] ,sep="",".pdf"), paper = "a4r")
print(type[i])
#dev.copy2pdf()
#导出文件
write.table(me$trans$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_trans_eqtls.txt",sep=""),sep = "\t",quote = F)
write.table(me$cis$eqtls,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte/EQTL/",type[i],"_cis_eqtls.txt",sep=""),sep = "\t",quote = F)

}

#已经跑完
#将eqtl得到的结果提取FDR<0.05的用于下面画图
#rm(list=ls())
setwd("/path/singlecell_eqtl/eqtl/result/result_monocyte1/EQTL/cis/")
eqtls<-list.files("/path/singlecell_eqtl/eqtl/result/result_monocyte1/EQTL/cis/",pattern="*_cis_eqtls.txt")

eqtls


for(i in 1:length(eqtls))
{

  eqtl<-read.table(file=eqtls[i],header=TRUE,sep="\t")
  eqtl_new<-  eqtl[ which(eqtl$FDR<0.05), ]
 write.table(eqtl_new,paste0("/path/singlecell_eqtl/eqtl/result/result_monocyte1/upset_pro/",
       eqtls[i], sep=""),quote=F,row.names=F,col.names=T,sep="\t")

  print(c(eqtls[i],dim(eqtl_new)))
 }

1



















#long long long time
diff <- differentialGeneTest(monocle_cds[express_genes,],fullModelFormulaStr="~subtype",cores=1)  #Monocle做差异分析的主要方法
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)

deg1=diff

##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
#deg <- subset(diff, qval < 0.05) #选出2829个基因
deg=diff
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

##差异基因的结果文件保存
write.table(deg,file="Monocyte_train.monocle.DEG.txt",col.names=T,row.names=F,sep="\t",quote=F)



## 轨迹构建基因可视化
ordergene <- rownames(deg) 
monocle_cds <- setOrderingFilter(monocle_cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
pdf("train.ordergenes.pdf")
plot_ordering_genes(monocle_cds)
dev.off()
plot_ordering_genes(monocle_cds)
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

ordergene <- row.names(deg)[order(deg$qval)][1:2000]
ordergene

monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                                  method = 'DDRTree')
#一旦细胞有序排列，我们就可以在降维空间中可视化轨迹。所以首先选择用于细胞排序的基因，然后使用反向图嵌入(DDRTree)算法对数据进行降维。

# long long time
#num_paths 生物过程中允许的终点细胞状态的数量
#root_cell 要用作排序树根的单元格的名称4
#reverse   是否颠倒所学生物过程的起点和终点
monocle_cds <- orderCells(monocle_cds,reverse = T)    #将细胞排序并完成轨迹构建
#使用root_state参数可以设置拟时间轴的根，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点

plot_size(8,8)
plot_cell_trajectory(monocle_cds,color_by="subtype", 
                     size=1,show_backbone=TRUE,root_states = 4)

plot_size(8,8)
#pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
plot_cell_trajectory(monocle_cds,color_by="State", 
                     size=1,show_backbone=TRUE)
#dev.off()

有时候细胞群不会完美的按照分叉排列。如下，只要一类细胞占某一个分叉细胞的大部分也可以

plot_size(8,8)
plot_cell_trajectory(monocle_cds,color_by="Pseudotime", 
                     size=1,show_backbone=TRUE)
#以pseudotime值上色 (Pseudotime是monocle2基于细胞基因表达信息计算的概率，表示时间的先后。

table(pData(monocle_cds)$State)



View(pData(monocle_cds))

plot_size(8,8)
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters",root_states = 4)

plot_size(8,8)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
#colour=c("#F4C781","#42B1D3","#B54E93","#5FB46F")
fig4A1<-plot_complex_cell_trajectory(monocle_cds,x=1,y=2,color_by = "subtype",
                                     show_backbone=TRUE)+#,root_states = 4
  scale_colour_manual(values =colour)
fig4A1

monocle_cds$Status<-factor(monocle_cds$Status,levels =c("Healthy","asymptomatic case","Moderate","Recover","Severe"))

plot_size(8,8)
colour=c("#DC143C","#9999FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
fig4A2<-plot_complex_cell_trajectory(monocle_cds,x=1,y=2,color_by = "Status",
                                     show_backbone=TRUE)+#root_states = 4)+
  scale_colour_manual(values = colour)
fig4A2

plot_size(10,7)
ggarrange(fig4A1,fig4A2, ncol = 2, nrow = 1)

plot_size(10,8)
#pdf("train.monocle.state.faceted.pdf",width = 10,height = 7)
plot_cell_trajectory(monocle_cds, color_by = "State") + facet_wrap("~State", nrow = 1)
#dev.off()

plot_size(20,8)
p1=plot_cell_trajectory(monocle_cds, color_by = "subtype")  + scale_color_npg() 
p2=plot_cell_trajectory(monocle_cds, color_by = "Status")  + scale_color_nejm()
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
#p3=plot_cell_trajectory(monocle_cds, color_by = "State")  + scale_color_manual(values = colour)
p1|p2

p1 <- plot_cell_trajectory(monocle_cds, x = 1, y = 2, color_by = "subtype") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_manual(values = colour) 
p2 <- plot_complex_cell_trajectory(monocle_cds, x = 1, y = 2,
                                   color_by = "subtype")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 
p1|p2

df <- pData(monocle_cds) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
View(pData(monocle_cds))


plot_size(15,8)
ggplot(df, aes(Pseudotime, colour = subtype, fill=subtype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()



















#ClusterName_color_panel <- c(
  "Naive CD4 T" = "#DC143C", "Memory CD4 T" = "#0000FF", "CD14+ Mono" = "#20B2AA",
  "B" = "#FFA500", "CD8 T" = "#9370DB", "FCGR3A+ Mono" = "#98FB98",
  "NK" = "#F08080", "DC" = "#0000FF", "Platelet" = "#20B2AA"
#)

ClusterName_color_panel <- c(
  "Monocyte classical" = "#DC143C", "Monocyte nonclassical" = "#0000FF"
)


ggplot(df, aes(Pseudotime, colour = subtype, fill=subtype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+ scale_fill_manual(name = "", values = ClusterName_color_panel)+scale_color_manual(name = "", values = ClusterName_color_panel)













#比如对State2的细胞感兴趣
pdata <- Biobase::pData(monocle_cds)
s.cells <- subset(pdata, State=="2") %>% rownames()
save(s.cells, file = "/path/singlecell_eqtl/result/monocle/monocyte/Monocle_state5.rda")

write.csv(pData(monocle_cds), "/path/singlecell_eqtl/result/monocle/monocyte/monocle_cds_pseudotime.csv")
save(monocle_cds, file = "/path/singlecell_eqtl/result/monocle/monocyte/monocle_cds.rda")

keygenes <- head(ordergene,4)
monocle_cds_subset <- monocle_cds[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "subtype")
p3 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
plotc
#ggsave("Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)

plot_size(20,8)
DotPlot(Monocyte,c(     

  
    "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",#Monocytes   
  "CD163" ,  "CD16","LYZ","S100A12"


   
   
    
    
                      
         ),group.by='Type',assay='SCT')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

keygenes <- c("RPS26")
#keygenes <- c("FCGR3A",  "CD68" , "CD14" , "FCN1","LYZ")
monocle_cds_subset <- monocle_cds[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "subtype")
p3 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
plotc
#ggsave("Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)

View(pData(monocle_cds))

table((pData(monocle_cds))$seurat_clusters)
table((pData(monocle_cds))$State)

plot_size(15,8)
keygenes <- c("FCGR3A",  "CD68" , "CD14" , "FCN1","LYZ","RPS26")
monocle_cds_subset <- monocle_cds[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "Status")
p2 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "subtype")
p3 <- plot_genes_in_pseudotime(monocle_cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
plotc

s.genes <- c("FCGR3A",  "CD68" , "CD14" , "FCN1","LYZ","RPS26")
p1 <- plot_genes_jitter(monocle_cds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(monocle_cds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(monocle_cds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
plotc
#ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 8)

colnames(pData(monocle_cds))
pData(monocle_cds)$CD14 = log2(exprs(monocle_cds)['CD14',]+1)
p1=plot_cell_trajectory(monocle_cds, color_by = "CD14")    + scale_color_gsea()
pData(monocle_cds)$FCGR3A = log2( exprs(monocle_cds)['FCGR3A',]+1)
p2=plot_cell_trajectory(monocle_cds, color_by = "FCGR3A")  + scale_color_gsea()
pData(monocle_cds)$RPS26 = log2( exprs(monocle_cds)['RPS26',]+1)
p3=plot_cell_trajectory(monocle_cds, color_by = "RPS26")  + scale_color_gsea()


library(patchwork)
p1+p2+p3

Monocle的主要工作是通过生物过程（如细胞分化）将细胞按顺序排列，而不知道要提前查看哪些基因。一旦这样做了，你就可以分析细胞，找到随着细胞进展而变化的基因。

官方给出的差异分析有三大方法：
1、Basic Differential Analysis
2、Finding Genes that Distinguish Cell Type or State
3、Finding Genes that Change as a Function of Pseudotime
我们重点关注第三个：根据伪时间功能寻找差异基因

sm.ns函数指出Monocle应该通过表达式值拟合自然样条曲线，以帮助它将表达式的变化描述为进程的函数


length(rownames(monocle_cds))

subset <- monocle_cds[which( rownames(monocle_cds)  %in%  ordergene), ]
subset

#long long time
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
#Time_diff <- differentialGeneTest(monocle_cds[ordergene,], cores = 1, 
                                 # fullModelFormulaStr = "~sm.ns(Pseudotime)")
#Time_diff <- Time_diff[,c(5,2,3,4,1,6)] #把gene放前面，也可以不改
#write.csv(Time_diff, "Time_diff_all.csv", row.names = F)



#long long time
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(subset, cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6)] #把gene放前面，也可以不改
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)


dim(Time_diff)
Time_diff

Time_diff1 <- Time_diff[1:50,]
Time_diff1



plot_size(10,10)
Time_genes <- Time_diff1 %>% pull(gene_short_name) %>% as.character()
Time_genes <- unique(Time_genes)
p=plot_pseudotime_heatmap(monocle_cds[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
#cluster_rows 是否对热图的行进行聚类;show_rownames 是否显示表中每一行的名称;return_heatmap 是否将热图对象返回给用户
#num_clusters 分支基因热图的簇数;#hmcols 用于绘制热图的配色方案
#add_annotation_row 要为热图中的每一行显示的其他批注。必须是数据帧，其中 fData 表中cds_subset的每一行对应一行，具有匹配的 ID。
#add_annotation_col 要为热图中的每一列显示的其他注释。必须是数据帧，其中 pData 表中的每个单元格对应一行cds_subset，具有匹配的 ID。
p
#ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

这个图的横轴是拟时间，cluster1的基因是在拟时排序起点高表达的基因，cluster2的基因则是在拟时排序的重点高表达的。
cluster数的多少是由plot_pseudotime_heatmap函数中的num_clusters参数定义的

p$tree_row
clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
#write.csv(clustering, "Time_clustering_all.csv", row.names = F)

Time_genes <- top_n(Time_diff, n = 50, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(monocle_cds[Time_genes,], num_clusters=1, show_rownames=T, return_heatmap=T)
p
#ggsave("Time_heatmapTop100.pdf", p, width = 5, height = 10)

p$tree_row$labels[p$tree_row$order]

View(p$tree_row$labels)

hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
#write.csv(Time_diff_sig, "Time_diff_sig.csv", row.names = F)
hp.genes
Time_diff_sig

marker_genes <- row.names(subset(fData(monocle_cds),
                                 gene_short_name %in% c(    #Monocytes 
                                     "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",
                                     "CD163" ,"CD16","LYZ","S100A12" )))

diff_test_res <- differentialGeneTest(monocle_cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
diff_test_res
sig_gene_names
marker_genes

plot_size(10,8)
plot_pseudotime_heatmap(monocle_cds[sig_gene_names,],
                        num_clusters = 2,
                        cores = 1,
                        show_rownames = T)

单细胞轨迹常常包括分支。这些分支的产生是因为细胞执行不同的基因表达程序。
在发育过程中，当细胞做出命运选择时，分支出现在轨迹中：一个发育谱系沿着一条路径前进，而另一个谱系产生第二条路径。
Monocle包含分析这些分支事件的广泛功能。Monocle提供了一个特殊的统计测试：分支表达式分析建模，或BEAM
#BEAM(Branched expression analysis modeling)是一种统计方法，用于寻找以依赖于分支的方式调控的基因

plot_cell_trajectory(monocle_cds, color_by = "State")

subset <- monocle_cds[which( rownames(monocle_cds)  %in%  ordergene), ]
subset

#long long time
BEAM_res <- BEAM(subset, branch_point = 1, cores = 2) 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)

#write.csv(BEAM_res, "BEAM_res.csv", row.names = F)


a<-BEAM_res[order(-BEAM_res$qval),]
a[1:100,]

p<-plot_genes_branched_heatmap(monocle_cds[row.names(subset(BEAM_res,
                                                  pval < 0.05)),],
                            branch_point = 2, #绘制的是哪个分支
                            num_clusters = 7, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
#use_gene_short_name 是否对每行使用短名称。如果为 FALSE，则使用 fData 表中的行 ID。
p

该热图显示的是同一时间点两个谱系的变化，热图的列是伪时间的点，行是基因。
这张图最上面的条条，灰色的代表分叉前，左边红色代表左边这个cell fate，右边蓝色代表右边这个cell fate，
从热图中间往右读，是伪时间的一个谱系，往左是另一个谱系。基因是被按照等级聚类的，需要结合生物学知识来进行解读。
因为前面设置的branch_point = 1，根据按照State绘制的trajectory图，在1这个分枝上分出的是1和2(好像)这两个status

plot_size(10,25)
BEAM_genes <- top_n(BEAM_res, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(monocle_cds[BEAM_genes,],  branch_point = 2, 
                                 num_clusters = 7, show_rownames = T, return_heatmap = T)
p$ph_res
#ggsave("BEAM_heatmap.pdf", p$ph_res, width = 6.5, height = 10)

plot_size(10,15)
#plot_genes_branched_heatmap(monocle_cds[BEAM_genes,],  branch_point = 1, 
#                                num_clusters = 3, show_rownames = T, return_heatmap = T)
p$ph_res

#显著差异基因(top100)按热图结果排序并保存
##如果要所有的差异基因，就把前面所632个基因的热图存为p
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- BEAM_res[hp.genes, c("gene_short_name", "pval", "qval")]
#write.csv(BEAM_sig, "BEAM_sig.csv", row.names = F)
hp.genes
BEAM_sig



genes <- row.names(subset(fData(monocle_cds),
                          gene_short_name %in% c(  "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R","CD163" ,"CD16","LYZ","S100A12")))

head(genes)

plot_size(10,15)
plot_genes_branched_pseudotime(monocle_cds[genes,],
                               branch_point = 2,
                               color_by = "State",
                               ncol = 1)

注意：不同分组间的细胞尽量不要放在一起做轨迹分析，同一组的生物学重复可以一起分析。明显没有生物学相关性的细胞也不要放在一起做轨迹分析。

可以看到，在pre-branch细胞向这两个分叉分化时，这两个基因都是在state1中高表达的，也就是state1的分支相关基因。

#plot_multiple_branches_pseudotime(cds, branches, branches_name = NULL,min_expr = NULL, cell_size = 0.75, norm_method = c("vstExprs", "log"),nrow = NULL, ncol = 1, panel_order = NULL, color_by = "Branch",
#trend_formula = "~sm.ns(Pseudotime, df=3)", label_by_short_name = TRUE,TPM = FALSE, cores = 1)
#示范命令
#plot_multiple_branches_heatmap(celltrajectory.monocle, branches = c(6,7),
#cluster_rows = TRUE, hclust_method = "ward.D2", num_clusters = 6,
#hmcols = NULL, add_annotation_row = NULL, add_annotation_col = NULL,
#show_rownames = FALSE, use_gene_short_name = TRUE,
#norm_method = c("vstExprs", "log"), scale_max = 3, scale_min = -3,
#trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
#cores = 1)

plot_multiple_branches_heatmap(monocle_cds[BEAM_genes,], branches = c(1,2,3),
cluster_rows = TRUE, hclust_method = "ward.D2", num_clusters = 7,
hmcols = NULL, add_annotation_row = NULL, add_annotation_col = NULL,
show_rownames =T, use_gene_short_name = TRUE,
norm_method = c("vstExprs", "log"), scale_max = 3, scale_min = -3,
trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
cores = 1)

genes <- row.names(subset(fData(monocle_cds),
                          gene_short_name %in% c(  "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R","CD163" ,"CD16","LYZ","S100A12")))



plot_multiple_branches_heatmap(monocle_cds[BEAM_genes,], branches = c(1,2),
cluster_rows = TRUE, hclust_method = "ward.D2", num_clusters = 7,
hmcols = NULL, add_annotation_row = NULL, add_annotation_col = NULL,
show_rownames = T, use_gene_short_name = TRUE,
norm_method = c("vstExprs", "log"), scale_max = 3, scale_min = -3,
trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
cores = 1)

 num_branches(monocle_cds)


pData(monocle_cds)




