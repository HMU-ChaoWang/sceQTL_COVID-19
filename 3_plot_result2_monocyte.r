library(Seurat)
library(ggplot2)
library(harmony)
library(ggsci)
#library(DoubletFinder)
library(dplyr)
library(patchwork)#调整排版设置
library(tidyverse)
library(ggpubr)
library(reshape2)
library(ggpubr)
#library(progeny)
library(gplots)
library(SingleR)
library(scRNAtoolVis)
library(hdf5r)
library(COSG)
library(RColorBrewer) 
library(colorspace)
library(viridis)
library(scales)
library(data.table)
library("UpSetR")
library(grid)
library(vcfR)
library(plyr)
library(PupillometryR)
library(ComplexHeatmap)
library("BuenColors")
library(ggsignif)
library(RImagePalette)
library(VennDiagram)
library(pheatmap)
library(ggrepel)
library(grid)
library(gridExtra)

plot_size <- function(x, y) {
    options(repr.plot.width = x, repr.plot.height = y)
}

path = "/path/"
files <- list.files(path)
#这一步配合下面的 Read10X(dir[i]) 能给barcode添加前缀
dir = paste0(path,'/',files)
names(dir) <- files

#显示色板名字和颜色数量
#maxcolors，调色板中的颜色数
#category，调色板分类，有：div, qual, seq 三种
#colorblind，对色盲是否友好
#连续型（sequential）：单渐变色，一种颜色由浅到深。
#离散型（divergent）：双渐变色，一种颜色到另外一种颜色的渐变。
#定性型（qualitative）：区分色，几种区分度很高的颜色组合。
rownames(brewer.pal.info)
brewer.pal.info
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=FALSE)
#display.brewer.pal(8,"Set1")
brewer.pal(9,"Set1")





plot_size (5,5)
setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3")
scales::show_col(setbarcolor)

setwd("/path/")
fileFDR0.05_path <- "/path/"
eqtl_files<-list.files(fileFDR0.05_path)
eqtl_files

eGene<-list()
eSNP<-list()
eQTL<-list()
 
    for(i in 1:length(eqtl_files)){
        file<-read.table(paste0(fileFDR0.05_path,eqtl_files[i]),header = T,sep="\t")
        eGene[[i]]<-file$gene
        eSNP[[i]]<-file$snps
        eQTL[[i]]<-paste0(file$snps,"_",file$gene)
        print(i)
    }

    listInput_eGene_Tissue_T<-list(  
                            #"Macrophage_NLRP3_N"=eGene[[9]],"Macrophage_NLRP3_PBMC"=eGene[[10]],
                            "1"=eGene[[1]],
                            #"TAM_C1QA_N"=eGene[[12]],"TAM_C1QA_PBMC"=eGene[[13]],
                            "2"=eGene[[2]],
                            #"TAM_MKI67_N"=eGene[[15]],
                            "3"=eGene[[3]],
                            #"TAM_UBB_N"=eGene[[17]],"TAM_UBB_PBMC"=eGene[[18]],
                            "4"=eGene[[4]],
"5"=eGene[[5]])


    listInput_eSNP_Tissue_T<-list(  
                            #"Macrophage_NLRP3_N"=eSNP[[9]],"Macrophage_NLRP3_PBMC"=eSNP[[10]],
                            "1"=eSNP[[1]],
                            #"TAM_C1QA_N"=eSNP[[12]],"TAM_C1QA_PBMC"=eSNP[[13]],
                            "2"=eSNP[[2]],
                            #"TAM_MKI67_N"=eSNP[[15]],
                            "3"=eSNP[[3]],
                            #"TAM_UBB_N"=eSNP[[17]],"TAM_UBB_PBMC"=eSNP[[18]],
                            "4"=eSNP[[4]],
"5"=eSNP[[5]])

    listInput_eQTL_Tissue_T<-list(  
                            #"Macrophage_NLRP3_N"=eQTL[[9]],"Macrophage_NLRP3_PBMC"=eQTL[[10]],
                            "1"=eQTL[[1]],
                            #"TAM_C1QA_N"=eQTL[[12]],"TAM_C1QA_PBMC"=eQTL[[13]],
                            "2"=eQTL[[2]],
                            #"TAM_MKI67_N"=eQTL[[15]],
                            "3"=eQTL[[3]],
                            #"TAM_UBB_N"=eQTL[[17]],"TAM_UBB_PBMC"=eQTL[[18]],
                            "4"=eQTL[[4]],
"5"=eQTL[[5]])

options(repr.plot.width = 16, repr.plot.height = 12)

#setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3","#BC8F8F","#20B2AA","#08519c")
setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3")
#setbarcolor<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")#"#FFFF33","#A65628","#F781BF","#999999"

upset(fromList(listInput_eGene_Tissue_T),
      sets.bar.color = setbarcolor, #柱子颜色
      nsets = 100,     # 绘制的最大集合个数
     # nsets =dim(fromList(listInput_eGene_Tissue_T))[2],#即可统计数据集的个数
      nintersects = NA, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据集合的大小排序。
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE
      mb.ratio = c(0.7,0.3),   # 左侧和上方条形图的比例关系
      text.scale = 1.5, # 文字标签的大小
      number.angles = 0, #柱子角度
      point.size = 3, #点的大小
      line.size = 1,
      main.bar.color = "#C5C9C6", #交集柱子颜色
      matrix.color = "#292D7C", #矩阵点的颜色
      mainbar.y.label = "eGene Intersections",
      sets.x.label = "cell types Intersections",

      queries = list(
                    list(query = intersects,params = list("1"),
                          color="#DB1E25",active = T),
          
                    list(query = intersects,params = list("2"),
                          color="#E76124",active = T),
          
                   
                    list(query = intersects,params = list("3"),
                          color="#F5C040",active = T),
      
                  
                    list(query = intersects,params = list("4"),
                          color="#F2932A",active = T),
          list(query = intersects,params = list("5"),
                          color="#FF34B3",active = T))

     )
dev.copy2pdf( file="/path//monocyte_eGene.pdf",paper = "a4r")

inter_eGene <- get.venn.partitions(listInput_eGene_Tissue_T)
inter_eGene



for (i in 1:nrow(inter_eGene)) inter_eGene[i,'values'] <- paste(inter_eGene[[i,'..values..']], collapse = ', ')
write.table(inter_eGene[-c(6, 7)],paste0("/path//","Monocyte_eGene_inter.txt",sep=""), row.names = FALSE, sep = '\t', quote = FALSE)
inter_eGene

options(repr.plot.width = 16, repr.plot.height = 12)

#setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3","#BC8F8F","#20B2AA","#08519c")
setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3")
#setbarcolor<-c( '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999')

upset(fromList(listInput_eSNP_Tissue_T),
      sets.bar.color = setbarcolor, #柱子颜色
      nsets = 100,     # 绘制的最大集合个数
      nintersects = NA, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据集合的大小排序。
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE
      mb.ratio = c(0.7,0.3),   # 左侧和上方条形图的比例关系
      text.scale = 1.5, # 文字标签的大小
      number.angles = 0, #柱子角度
      point.size = 3, #点的大小
      line.size = 1,
      main.bar.color = "#C5C9C6", #交集柱子颜色
      matrix.color = "#292D7C", #矩阵点的颜色
      mainbar.y.label = "eGene Intersections",
      sets.x.label = "cell types Intersections",

      queries = list(
                     list(query = intersects,params = list("1"),
                          color="#F2932A",active = T),
          
                    list(query = intersects,params = list("2"),
                          color="#E76124",active = T),
          
                   
                    list(query = intersects,params = list("3"),
                          color="#F5C040",active = T),
      
                  
                    list(query = intersects,params = list("4"),
                          color="#DB1E25",active = T),
          list(query = intersects,params = list("5"),
                          color="#FF34B3",active = T))

     )
dev.copy2pdf( file="/path//monocyte_eSNP.pdf",paper = "a4r")

inter_eSNP <- get.venn.partitions(listInput_eSNP_Tissue_T)
inter_eSNP



for (i in 1:nrow(inter_eSNP)) inter_eSNP[i,'values'] <- paste(inter_eSNP[[i,'..values..']], collapse = ', ')
write.table(inter_eSNP[-c(6, 7)],paste0("/path//","Monocyte_eSNP_inter.txt",sep=""), row.names = FALSE, sep = '\t', quote = FALSE)
inter_eSNP

options(repr.plot.width = 16, repr.plot.height = 12)

#setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3","#BC8F8F","#20B2AA","#08519c")
setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3")
#setbarcolor<-c( '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999')

upset(fromList(listInput_eQTL_Tissue_T),
      sets.bar.color = setbarcolor, #柱子颜色
      nsets = 100,     # 绘制的最大集合个数
      nintersects = NA, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据集合的大小排序。
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE
      mb.ratio = c(0.7,0.3),   # 左侧和上方条形图的比例关系
      text.scale = 1.5, # 文字标签的大小
      number.angles = 0, #柱子角度
      point.size = 3, #点的大小
      line.size = 1,
      main.bar.color = "#C5C9C6", #交集柱子颜色
      matrix.color = "#292D7C", #矩阵点的颜色
      mainbar.y.label = "eGene Intersections",
      sets.x.label = "cell types Intersections",

      queries = list(
                    
                    list(query = intersects,params = list("1"),
                          color="#F5C040",active = T),
          
                    list(query = intersects,params = list("2"),
                          color="#E76124",active = T),
          
                   
                    list(query = intersects,params = list("3"),
                          color="#F2932A",active = T),
      
                  
                    list(query = intersects,params = list("4"),
                          color="#DB1E25",active = T),
          list(query = intersects,params = list("5"),
                          color="#FF34B3",active = T))

     )
dev.copy2pdf( file="/path//monocyte_eQTL.pdf",paper = "a4r")

inter_eQTL <- get.venn.partitions(listInput_eQTL_Tissue_T)
inter_eQTL

for (i in 1:nrow(inter_eQTL)) inter_eQTL[i,'values'] <- paste(inter_eQTL[[i,'..values..']], collapse = ', ')
write.table(inter_eQTL[-c(6, 7)],paste0("/path//","Monocyte_eQTL_inter.txt",sep=""), row.names = FALSE, sep = '\t', quote = FALSE)
inter_eQTL



Monocyte_exp_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/result_monocyte/exp/txt/"
subtype_exp_ststus<-list.files(Monocyte_exp_path)

state_1 <- read.table( paste0(Monocyte_exp_path,subtype_exp_ststus[1]),header=T,row.names = 1)
state_2      <- read.table( paste0(Monocyte_exp_path,subtype_exp_ststus[2]),header=T,row.names = 1)
state_3     <- read.table( paste0(Monocyte_exp_path,subtype_exp_ststus[3]),header=T,row.names = 1)
state_4      <- read.table( paste0(Monocyte_exp_path,subtype_exp_ststus[4]),header=T,row.names = 1)
state_5       <- read.table( paste0(Monocyte_exp_path,subtype_exp_ststus[5]),header=T,row.names = 1)

subtype_exp_ststus
state_5

cellsubtype_genotype_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/result_monocyte/EQTL_genotype/"
cellsubtype_all <- list.files(cellsubtype_genotype_path)

genotype_transpose_1<-fread(paste0(cellsubtype_genotype_path,cellsubtype_all[1],sep=""),
                                                 sep=" ",data.table=F)
#asymptomatic_genotype_transpose<-column_to_rownames(asymptomatic_genotype_transpose,var="V1")
colnames(genotype_transpose_1)[1]<-"ID"

genotype_transpose_2<-fread(paste0(cellsubtype_genotype_path,cellsubtype_all[1],sep=""),
                                                 sep=" ",data.table=F)
#Healthy_genotype_transpose<-column_to_rownames(Healthy_genotype_transpose,var="V1")
colnames(genotype_transpose_2)[1]<-"ID"


genotype_transpose_3<-fread(paste0(cellsubtype_genotype_path,cellsubtype_all[2],sep=""),
                                                 sep=" ",data.table=F)
#Moderate_genotype_transpose<-column_to_rownames(Moderate_genotype_transpose,var="V1")
colnames(genotype_transpose_3)[1]<-"ID"

genotype_transpose_4<-fread(paste0(cellsubtype_genotype_path,cellsubtype_all[3],sep=""),
                                                 sep=" ",data.table=F)
#Recover_genotype_transpose<-column_to_rownames(Recover_genotype_transpose,var="V1")
colnames(genotype_transpose_4)[1]<-"ID"

genotype_transpose_5<-fread(paste0(cellsubtype_genotype_path,cellsubtype_all[1],sep=""),
                                                 sep=" ",data.table=F)
#Severe_genotype_transpose<-column_to_rownames(Severe_genotype_transpose,var="V1")
colnames(genotype_transpose_5)[1]<-"ID"

genotype_transpose_1

eQTL_expression<-function(Gene,SNP){

    ###基因型数据###
    SNP_genotype_1 <- genotype_transpose_1[which(genotype_transpose_1$ID == SNP),] #定位SNP
    SNP_genotype_1<- SNP_genotype_1[,2:dim(SNP_genotype_1)[2]] #提取需要样本的行数
    SNP_genotype_1 <- t(as.data.frame(SNP_genotype_1))
    
    SNP_genotype_2 <- genotype_transpose_2[which(genotype_transpose_2$ID == SNP),] #定位SNP
    SNP_genotype_2<- SNP_genotype_2[,2:dim(SNP_genotype_2)[2]] #提取需要样本的行数
    SNP_genotype_2 <- t(as.data.frame(SNP_genotype_2))   
    
    SNP_genotype_3 <- genotype_transpose_3[which(genotype_transpose_3$ID == SNP),] #定位SNP
    SNP_genotype_3<- SNP_genotype_3[,2:dim(SNP_genotype_3)[2]] #提取需要样本的行数
    SNP_genotype_3 <- t(as.data.frame(SNP_genotype_3))

    SNP_genotype_4 <- genotype_transpose_4[which(genotype_transpose_4$ID == SNP),] #定位SNP
    SNP_genotype_4<- SNP_genotype_4[,2:dim(SNP_genotype_4)[2]] #提取需要样本的行数
    SNP_genotype_4 <- t(as.data.frame(SNP_genotype_4))
    
    SNP_genotype_5 <- genotype_transpose_5[which(genotype_transpose_5$ID == SNP),] #定位SNP
    SNP_genotype_5<- SNP_genotype_5[,2:dim(SNP_genotype_5)[2]] #提取需要样本的行数
    SNP_genotype_5 <- t(as.data.frame(SNP_genotype_5))
    
     
    ###表达数据###  
    Monocyte_gene_exp_1 <- state_1[which(rownames(state_1) == Gene),] #定位基因表达信息
    Monocyte_gene_exp_1 <- rbind(Monocyte_gene_exp_1,"Asymptomatic")
    Monocyte_gene_exp_1 <- t(Monocyte_gene_exp_1)
    Monocyte_gene_exp_1 <- as.data.frame(Monocyte_gene_exp_1)

    Monocyte_gene_exp_2 <- state_2[which(rownames(state_2) == Gene),] #定位基因表达信息
    Monocyte_gene_exp_2 <- rbind(Monocyte_gene_exp_2,"Healthy")
    Monocyte_gene_exp_2 <- t(Monocyte_gene_exp_2)
    Monocyte_gene_exp_2 <- as.data.frame(Monocyte_gene_exp_2)
    
    Monocyte_gene_exp_3 <- state_3[which(rownames(state_3) == Gene),] #定位基因表达信息
    Monocyte_gene_exp_3 <- rbind(Monocyte_gene_exp_3,"Moderate")
    Monocyte_gene_exp_3 <- t(Monocyte_gene_exp_3)
    Monocyte_gene_exp_3 <- as.data.frame(Monocyte_gene_exp_3)

    Monocyte_gene_exp_4 <- state_4[which(rownames(state_4) == Gene),] #定位基因表达信息
    Monocyte_gene_exp_4 <- rbind(Monocyte_gene_exp_4,"Recover")
    Monocyte_gene_exp_4 <- t(Monocyte_gene_exp_4)
    Monocyte_gene_exp_4 <- as.data.frame(Monocyte_gene_exp_4)
    
    Monocyte_gene_exp_5 <- state_5[which(rownames(state_5) == Gene),] #定位基因表达信息
    Monocyte_gene_exp_5 <- rbind(Monocyte_gene_exp_5,"Severe")
    Monocyte_gene_exp_5 <- t(Monocyte_gene_exp_5)
    Monocyte_gene_exp_5 <- as.data.frame(Monocyte_gene_exp_5)
    
    ###构建eqtl表达数据集###
    Monocyte_eqtl_expression_data_1 <- cbind(Monocyte_gene_exp_1,SNP_genotype_1)
    Monocyte_eqtl_expression_data_2 <- cbind(Monocyte_gene_exp_2,SNP_genotype_2)
    Monocyte_eqtl_expression_data_3 <- cbind(Monocyte_gene_exp_3,SNP_genotype_3)
    Monocyte_eqtl_expression_data_4 <- cbind(Monocyte_gene_exp_4,SNP_genotype_4)
    Monocyte_eqtl_expression_data_5 <- cbind(Monocyte_gene_exp_5,SNP_genotype_5)
      
    Monocyte_eqtl_expression_data_1 <- cbind(Monocyte_gene_exp_1,SNP_genotype_1)
    colnames(Monocyte_eqtl_expression_data_1)<-c("expression","celltype","genotype")
    
    Monocyte_eqtl_expression_data_2 <- cbind(Monocyte_gene_exp_2,SNP_genotype_2)
    colnames(Monocyte_eqtl_expression_data_2)<-c("expression","celltype","genotype")
    
    Monocyte_eqtl_expression_data_3 <- cbind(Monocyte_gene_exp_3,SNP_genotype_3)
    colnames(Monocyte_eqtl_expression_data_3)<-c("expression","celltype","genotype")
    
    Monocyte_eqtl_expression_data_4 <- cbind(Monocyte_gene_exp_4,SNP_genotype_4)
    colnames(Monocyte_eqtl_expression_data_4)<-c("expression","celltype","genotype")    
    
    Monocyte_eqtl_expression_data_5 <- cbind(Monocyte_gene_exp_5,SNP_genotype_5)
    colnames(Monocyte_eqtl_expression_data_5)<-c("expression","celltype","genotype") 
    
    allcellsubtype_expression_data<-rbind(Monocyte_eqtl_expression_data_1,
                                          Monocyte_eqtl_expression_data_2,
                                          Monocyte_eqtl_expression_data_3,
                                          Monocyte_eqtl_expression_data_4,
                                          Monocyte_eqtl_expression_data_5 )

    allcellsubtype_expression_data$expression<-as.numeric(allcellsubtype_expression_data$expression)
    allcellsubtype_expression_data$genotype<-as.character(allcellsubtype_expression_data$genotype)
    allcellsubtype_expression_data<-as.data.frame(allcellsubtype_expression_data)
    return(allcellsubtype_expression_data)
}  

rs705704_RPS26, rs6464100_TMEM176B, rs9614637_FAM118A, rs773111_RPS26, rs2074038_ACCS, rs310625_PPDPF!!!, rs10824071_ADK!!!, rs767809_ADK, rs738171_FAM118A, rs1139147_PPDPF!!!, rs2350630_FAM118A, rs226492_FAM118A, rs10528596;rs757062209;rs79050799_FAM118A, rs2549796_ERAP2, rs1056893_ERAP2, rs1230381_ERAP2, rs6464101_TMEM176B, rs1128966_NT5C3B!!!, rs1428835003;rs397804906_ADK, rs3829655_SIGLEC14!!!, rs2548525_ERAP2!!!, rs2549798_ERAP2, rs2549799_ERAP2, rs2113191_ERAP2, rs2548524_ERAP2, rs1019503_ERAP2, rs2548526_ERAP2!!!, rs3794623_CYBA, rs773109_RPS26, rs10546363;rs397746675_ERAP2

Gene<-c("STAT3")
SNP<-c("rs1128966")
expression_data<-eQTL_expression(Gene = Gene,SNP = SNP)
dim(expression_data)
expression_data

options(repr.plot.width = 12, repr.plot.height = 5)

expression_data$celltype<-factor(expression_data$celltype,
                                     levels=c("Asymptomatic","Healthy",
                                              "Moderate","Recover",
                                             "Severe"))



    ggplot(data = expression_data, aes(x = genotype, y = expression)) +
    geom_violin(aes(colour = celltype),trim=FALSE,width=0.8,size=1)+ 
    geom_boxplot(aes(colour = celltype),width=0.1,size = 0.5)+   
    geom_jitter(aes(colour = celltype ),width=0.1,shape=21,size=0.8) +
    scale_colour_manual(values = c("#F5C040","#F2932A",
                                   "#E76124","#DB1E25","#999999"))+
    theme_test() +
    labs(title=c(paste0(Gene," / ",SNP)),
        x = "Genotype",
        y = c(paste0("Expression of ",Gene)))+
    facet_grid(. ~ celltype)+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=14),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####去除图例###
        legend.position = 'none',
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill="#E8DAEF", size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1) #图形线框

     )

Gene<-c("STAT3")
SNP<-c("rs1128966")
expression_data<-eQTL_expression(Gene = Gene,SNP = SNP)
dim(expression_data)
expression_data

options(repr.plot.width = 12, repr.plot.height = 5)

expression_data$celltype<-factor(expression_data$celltype,
                                     levels=c("Asymptomatic","Health",
                              "Moderate","Recover","Severe"))


    ggplot(data = expression_data, aes(x = genotype, y = expression)) +
    geom_violin(aes(colour = celltype),trim=FALSE,width=0.8,size=1)+ 
    geom_boxplot(aes(colour = celltype),width=0.1,size = 0.5)+   
    geom_jitter(aes(colour = celltype ),width=0.1,shape=21,size=0.8) +
    scale_colour_manual(values = c("#F5C040","#F2932A",
                                   "#E76124","#DB1E25","#999999"))+
    theme_test() +
    labs(title=c(paste0(Gene," / ",SNP)),
        x = "Genotype",
        y = c(paste0("Expression of ",Gene)))+
    facet_grid(. ~ celltype)+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=14),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####去除图例###
        legend.position = 'none',
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill="#E8DAEF", size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1) #图形线框

     )

Monocyte_exp_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/4status_upset/Monocyte/type/"
ststus_exp_file<-list.files(Monocyte_exp_path)

asymptomatic_Monocyte <- paste0(Monocyte_exp_path,ststus_exp_file[1])
Healthy_Monocyte      <- paste0(Monocyte_exp_path,ststus_exp_file[2])
Moderate_Monocyte     <- paste0(Monocyte_exp_path,ststus_exp_file[3])
Recover_Monocyte      <- paste0(Monocyte_exp_path,ststus_exp_file[4])
Severe_Monocyte       <- paste0(Monocyte_exp_path,ststus_exp_file[5])

allfiles <-c(asymptomatic_Monocyte,Healthy_Monocyte,Moderate_Monocyte,Recover_Monocyte,Severe_Monocyte)
allfiles
subtype_exp_ststus
Monocyte_exp_path

Monocyte_eGene_num_status<-c()

for(i in 1:length(allfiles)){
    eqtl_file<-read.table(allfiles[i],header=T,sep="\t")
    eGene_num<-length(unique(eqtl_file$gene))
    status<-strsplit(strsplit(allfiles[i],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
    status_eGene_num<-cbind(status,eGene_num)
    Monocyte_eGene_num_status<-rbind(Monocyte_eGene_num_status,status_eGene_num)
    #All_eGene_num_type<-as.data.frame(cbind(c(celltype),All_eGene_num_type))
    print(i)
}

#colnames(All_eGene_num_type)<-c("type","subtype","eGene_num")
colnames(Monocyte_eGene_num_status)<-c("status","eGene_num")
Monocyte_eGene_num_status<-as.data.frame(Monocyte_eGene_num_status)
Monocyte_eGene_num_status

#All_result_eGene_num$sample_Freq<-53
#All_result_eGene_num

Monocyte_status_freq<-matrix(ncol = 3,nrow = 5)
colnames(Monocyte_status_freq)<-c("status","type","sample_Freq")
Monocyte_status_freq<-as.data.frame(Monocyte_status_freq)

Monocyte_status_freq$type<-"Monocyte"
Monocyte_status_freq$status<-c("asymptomatic_Monocyte","Healthy_Monocyte",
                              "Moderate_Monocyte","Recover_Monocyte","Severe_Monocyte")
Monocyte_status_freq$sample_Freq<-c(5,11,13,12,12)

Monocyte_status_freq

Monocyte_status_eGene_num<-cbind(Monocyte_status_freq,Monocyte_eGene_num_status$eGene_num)
colnames(Monocyte_status_eGene_num)<-c("status","type","sample_Freq","eGene_num")
Monocyte_status_eGene_num<-as.data.frame(Monocyte_status_eGene_num)
Monocyte_status_eGene_num$sample_Freq<-as.numeric(Monocyte_status_eGene_num$sample_Freq)
Monocyte_status_eGene_num$eGene_num<-as.numeric(Monocyte_status_eGene_num$eGene_num)
Monocyte_status_eGene_num

options(repr.plot.width = 10, repr.plot.height = 8)
ggplot(data=Monocyte_status_eGene_num,aes(x=sample_Freq,y=eGene_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=sample_Freq,y=eGene_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=14),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 20),
        legend.text=element_text(face="bold",size = 15)

     )

Monocyte_exp_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/4status_upset/Monocyte/type/"
ststus_exp_file<-list.files(Monocyte_exp_path)

asymptomatic_Monocyte <- paste0(Monocyte_exp_path,ststus_exp_file[1])
Healthy_Monocyte      <- paste0(Monocyte_exp_path,ststus_exp_file[2])
Moderate_Monocyte     <- paste0(Monocyte_exp_path,ststus_exp_file[3])
Recover_Monocyte      <- paste0(Monocyte_exp_path,ststus_exp_file[4])
Severe_Monocyte       <- paste0(Monocyte_exp_path,ststus_exp_file[5])

allfiles <-c(asymptomatic_Monocyte,Healthy_Monocyte,Moderate_Monocyte,Recover_Monocyte,Severe_Monocyte)
allfiles
subtype_exp_ststus
Monocyte_exp_path

Monocyte_eGene_num_status<-c()
i=1

    eqtl_file<-read.table(allfiles[i],header=T,sep="\t")
eqtl_file

Monocyte_esnp_num_status<-c()

for(i in 1:length(allfiles)){
    eqtl_file<-read.table(allfiles[i],header=T,sep="\t")
    esnp_num<-length(unique(eqtl_file$snps))
    status<-strsplit(strsplit(allfiles[i],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
    status_esnp_num<-cbind(status,esnp_num)
    Monocyte_esnp_num_status<-rbind(Monocyte_esnp_num_status,status_esnp_num)
    #All_eGene_num_type<-as.data.frame(cbind(c(celltype),All_eGene_num_type))
    print(i)
}

#colnames(All_eGene_num_type)<-c("type","subtype","eGene_num")
colnames(Monocyte_esnp_num_status)<-c("status","esnp_num")
Monocyte_esnp_num_status<-as.data.frame(Monocyte_esnp_num_status)
Monocyte_esnp_num_status

#All_result_eGene_num$sample_Freq<-53
#All_result_eGene_num

Monocyte_status_freq<-matrix(ncol = 3,nrow = 5)
colnames(Monocyte_status_freq)<-c("status","type","sample_Freq")
Monocyte_status_freq<-as.data.frame(Monocyte_status_freq)

Monocyte_status_freq$type<-"Monocyte"
Monocyte_status_freq$status<-c("asymptomatic_Monocyte","Healthy_Monocyte",
                              "Moderate_Monocyte","Recover_Monocyte","Severe_Monocyte")
Monocyte_status_freq$sample_Freq<-c(5,11,13,12,12)

Monocyte_status_freq

Monocyte_esnp_num_status<-cbind(Monocyte_status_freq,Monocyte_esnp_num_status$esnp_num)
colnames(Monocyte_esnp_num_status)<-c("status","type","sample_Freq","esnp_num")
Monocyte_esnp_num_status<-as.data.frame(Monocyte_esnp_num_status)
Monocyte_esnp_num_status$sample_Freq<-as.numeric(Monocyte_esnp_num_status$sample_Freq)
Monocyte_esnp_num_status$esnp_num<-as.numeric(Monocyte_esnp_num_status$esnp_num)
Monocyte_esnp_num_status

options(repr.plot.width = 10, repr.plot.height = 8)
ggplot(data=Monocyte_esnp_num_status,aes(x=sample_Freq,y=esnp_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=sample_Freq,y=esnp_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=14),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 20),
        legend.text=element_text(face="bold",size = 15)

     )

sce<-readRDS("/data_alluser/CXY/keti/singlecell_eqtl/result/RDS/all_subtype_annotated.rds")
sce

sce@meta.data

sce_Monocyte<-sce@meta.data[which(sce@meta.data$Type=="Monocyte"),]
sce_Monocyte

Monocyte_exp_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/4status_upset/Monocyte/type/"
ststus_exp_file<-list.files(Monocyte_exp_path)

asymptomatic_Monocyte <- paste0(Monocyte_exp_path,ststus_exp_file[1])
Healthy_Monocyte      <- paste0(Monocyte_exp_path,ststus_exp_file[2])
Moderate_Monocyte     <- paste0(Monocyte_exp_path,ststus_exp_file[3])
Recover_Monocyte      <- paste0(Monocyte_exp_path,ststus_exp_file[4])
Severe_Monocyte       <- paste0(Monocyte_exp_path,ststus_exp_file[5])

allfiles <-c(asymptomatic_Monocyte,Healthy_Monocyte,Moderate_Monocyte,Recover_Monocyte,Severe_Monocyte)
allfiles
subtype_exp_ststus
Monocyte_exp_path

Monocyte_eqtl_num_status<-c()

for(i in 1:length(allfiles)){
    eqtl_file<-read.table(allfiles[i],header=T,sep="\t")
   # esnp_num<-length(unique(eqtl_file$snps))
    eqtl_num<-dim(eqtl_file)[1]
    status<-strsplit(strsplit(allfiles[i],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
    status_eqtl_num<-cbind(status,eqtl_num)
    Monocyte_eqtl_num_status<-rbind(Monocyte_eqtl_num_status,status_eqtl_num)
    #All_eGene_num_type<-as.data.frame(cbind(c(celltype),All_eGene_num_type))
    print(i)
}

#colnames(All_eGene_num_type)<-c("type","subtype","eGene_num")
colnames(Monocyte_eqtl_num_status)<-c("status","eqtl_num")
Monocyte_eqtl_num_status<-as.data.frame(Monocyte_eqtl_num_status)
Monocyte_eqtl_num_status

#All_result_eGene_num$sample_Freq<-53
#All_result_eGene_num

Monocyte_status_freq<-matrix(ncol = 3,nrow = 5)
colnames(Monocyte_status_freq)<-c("status","type","sample_Freq")
Monocyte_status_freq<-as.data.frame(Monocyte_status_freq)

Monocyte_status_freq$type<-"Monocyte"
Monocyte_status_freq$status<-c("asymptomatic_Monocyte","Healthy_Monocyte",
                              "Moderate_Monocyte","Recover_Monocyte","Severe_Monocyte")
Monocyte_status_freq$sample_Freq<-c(5,11,13,12,12)

Monocyte_status_freq

#要换成总的EQTL数,没有FDR<0.05(手动打开文件看的)
Monocyte_eqtl_num_status$eqtl_num<-c(30303,35696,38324,38514,48454)
Monocyte_eqtl_num_status

Monocyte_eqtl_num_status<-cbind(Monocyte_status_freq,Monocyte_eqtl_num_status$eqtl_num)
colnames(Monocyte_eqtl_num_status)<-c("status","type","sample_Freq","eqtl_num")
Monocyte_eqtl_num_status<-as.data.frame(Monocyte_eqtl_num_status)
Monocyte_eqtl_num_status$sample_Freq<-as.numeric(Monocyte_eqtl_num_status$sample_Freq)
Monocyte_eqtl_num_status$eqtl_num<-as.numeric(Monocyte_eqtl_num_status$eqtl_num)
Monocyte_eqtl_num_status

Monocyte_eqtl_num_status$status_cellnum<-data.frame(table(sce_Monocyte$Status))[,2]
Monocyte_eqtl_num_status

options(repr.plot.width = 10, repr.plot.height = 8)
ggplot(data=Monocyte_eqtl_num_status,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=14),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 20),
        legend.text=element_text(face="bold",size = 15)

     )



























options(repr.plot.width = 12, repr.plot.height = 5)
ggplot(data=Monocyte_status_eGene_num,aes(x=status,y=eGene_num,fill=status))+
geom_bar(stat="identity",width=0.5)+facet_grid(. ~ type)+
theme_bw()+
scale_fill_manual(values=c("#C5CC62","#388246","#117A82","#1B66B0","#C57FCC","#892E73",
                           
                           "#F5C040","#F2932A","#E76124","#DB1E25"))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_blank(),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill=c("#EA9390","#682342"), size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1), #图形线框
        
        ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 12)
     )




Monocyte_exp_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/4status_upset/Monocyte/type/"
ststus_exp_file<-list.files(Monocyte_exp_path)

asymptomatic_Monocyte <- paste0(Monocyte_exp_path,ststus_exp_file[1])
Healthy_Monocyte      <- paste0(Monocyte_exp_path,ststus_exp_file[2])
Moderate_Monocyte     <- paste0(Monocyte_exp_path,ststus_exp_file[3])
Recover_Monocyte      <- paste0(Monocyte_exp_path,ststus_exp_file[4])
Severe_Monocyte       <- paste0(Monocyte_exp_path,ststus_exp_file[5])

allfiles <-c(asymptomatic_Monocyte,Healthy_Monocyte,Moderate_Monocyte,Recover_Monocyte,Severe_Monocyte)
allfiles
subtype_exp_ststus
Monocyte_exp_path

















options(repr.plot.width = 12, repr.plot.height = 5)
ggplot(data=Monocyte_esnp_num_status,aes(x=status,y=esnp_num,fill=status))+
geom_bar(stat="identity",width=0.5)+facet_grid(. ~ type)+
theme_bw()+
scale_fill_manual(values=c("#C5CC62","#388246","#117A82","#1B66B0","#C57FCC","#892E73",
                           
                           "#F5C040","#F2932A","#E76124","#DB1E25"))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_blank(),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill=c("#EA9390","#682342"), size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1), #图形线框
        
        ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 12)
     )


All_eQTL_num_ststus_AllCell<-fread("/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/2EQTL_result/All_eQTL_num_ststus_AllCell.txt")

All_eQTL_num_ststus_AllCell

 options(repr.plot.width = 10, repr.plot.height = 8)
ggplot(data=All_eQTL_num_ststus_AllCell,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=type),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )

options(repr.plot.width = 8, repr.plot.height = 5)

Monocyte
B
CD4
CD8
cDC
NK
pDC
Platelet

Monocyte<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$type=="Monocyte"),]
Monocyte

 options(repr.plot.width = 10, repr.plot.height = 8)
Monocyte<-ggplot(data=Monocyte,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
Monocyte

B<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$type=="B"),]
B

 options(repr.plot.width = 10, repr.plot.height = 8)
B<-ggplot(data=B,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
B

CD4<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$type=="CD4"),]
CD4

 options(repr.plot.width = 10, repr.plot.height = 8)
CD4<-ggplot(data=CD4,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
CD4

CD8<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$type=="CD8"),]
CD8

 options(repr.plot.width = 10, repr.plot.height = 8)
CD8<-ggplot(data=CD8,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
CD8

cDC<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$type=="cDC"),]
cDC

 options(repr.plot.width = 10, repr.plot.height = 8)
cDC<-ggplot(data=cDC,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
cDC

NK<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$type=="NK"),]
NK

 options(repr.plot.width = 10, repr.plot.height = 8)
NK<-ggplot(data=NK,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
NK

pDC<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$type=="pDC"),]
pDC

 options(repr.plot.width = 10, repr.plot.height = 8)
pDC<-ggplot(data=pDC,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
pDC

Platelet<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$type=="Platelet"),]
Platelet

 options(repr.plot.width = 10, repr.plot.height = 8)
Platelet<-ggplot(data=Platelet,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
Platelet

asymptomatic
Healthy
Moderate
Recover
Severe

asymptomatic<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$Status=="asymptomatic"),]
asymptomatic

 options(repr.plot.width = 10, repr.plot.height = 8)
asymptomatic<-ggplot(data=asymptomatic,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = type))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
asymptomatic

Healthy<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$Status=="Healthy"),]
Healthy

 options(repr.plot.width = 10, repr.plot.height = 8)
Healthy<-ggplot(data=Healthy,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = type))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
Healthy

Moderate<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$Status=="Moderate"),]
Moderate

 options(repr.plot.width = 10, repr.plot.height = 8)
Moderate<-ggplot(data=Moderate,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = type))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
Moderate

Recover<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$Status=="Recover"),]
Recover

 options(repr.plot.width = 10, repr.plot.height = 8)
Recover<-ggplot(data=Recover,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = type))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
Recover

Severe<-All_eQTL_num_ststus_AllCell[which(All_eQTL_num_ststus_AllCell$Status=="Severe"),]
Severe

 options(repr.plot.width = 10, repr.plot.height = 8)
Severe<-ggplot(data=Severe,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=status),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = type))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'#'#999999'

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )
Severe

setwd("/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_sub/")
fileFDR0.05_path <- "/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_sub/"
eqtl_files<-list.files(fileFDR0.05_path)
eqtl_files

eGene<-list()
eSNP<-list()
eQTL<-list()
 
    for(i in 1:length(eqtl_files)){
        file<-read.table(paste0(fileFDR0.05_path,eqtl_files[i]),header = T,sep="\t")
        eGene[[i]]<-file$gene
        eSNP[[i]]<-file$snps
        eQTL[[i]]<-paste0(file$snps,"_",file$gene)
        
    }

    listInput_eGene_Tissue_T<-list(  
                            #"Macrophage_NLRP3_N"=eGene[[9]],"Macrophage_NLRP3_PBMC"=eGene[[10]],
                            "B memory"=eGene[[1]],
                            #"TAM_C1QA_N"=eGene[[12]],"TAM_C1QA_PBMC"=eGene[[13]],
                            "B naive"=eGene[[2]],
                            "CD4 T memory"=eGene[[3]],
                            #"TAM_UBB_N"=eGene[[17]],"TAM_UBB_PBMC"=eGene[[18]],
                            "CD4 T naive"=eGene[[4]],
                            "CD8 T effector"=eGene[[5]],
                            "CD8 T naive"=eGene[[6]],
                            "cDC"=eGene[[7]],
                            "Monocyte classical"=eGene[[8]],
                            "Monocyte nonclassical"=eGene[[9]],      
                            "NK"=eGene[[10]],
                            "pDC"=eGene[[11]],
                            "Plasma"=eGene[[12]],
                            "Platelet"=eGene[[13]])


    listInput_eSNP_Tissue_T<-list(  
                            "B memory"=eSNP[[1]],
                            "B naive"=eSNP[[2]],
                            "CD4 T memory"=eSNP[[3]],
                            "CD4 T naive"=eSNP[[4]],
                            "CD8 T effector"=eSNP[[5]],
                            "CD8 T naive"=eSNP[[6]],
                            "cDC"=eSNP[[7]],
                            "Monocyte classical"=eSNP[[8]],
                            "Monocyte nonclassical"=eSNP[[9]],      
                            "NK"=eSNP[[10]],
                            "pDC"=eSNP[[11]],
                            "Plasma"=eSNP[[12]],
                            "Platelet"=eSNP[[13]])





                           
    listInput_eQTL_Tissue_T<-list(  
                            "B memory"=eQTL[[1]],
                            "B naive"=eQTL[[2]],
                            "CD4 T memory"=eQTL[[3]],
                            "CD4 T naive"=eQTL[[4]],
                            "CD8 T effector"=eQTL[[5]],
                            "CD8 T naive"=eQTL[[6]],
                            "cDC"=eQTL[[7]],
                            "Monocyte classical"=eQTL[[8]],
                            "Monocyte nonclassical"=eQTL[[9]],      
                            "NK"=eQTL[[10]],
                            "pDC"=eQTL[[11]],
                            "Plasma"=eQTL[[12]],
                            "Platelet"=eQTL[[13]])





                    

options(repr.plot.width = 16, repr.plot.height = 12)

#setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3","#BC8F8F","#20B2AA","#08519c")
   
setbarcolor<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c")

upset(fromList(listInput_eGene_Tissue_T),
      sets.bar.color = setbarcolor, #柱子颜色
      nsets = 100,     # 绘制的最大集合个数
      nintersects = 50, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据集合的大小排序。
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE
      mb.ratio = c(0.7,0.3),   # 左侧和上方条形图的比例关系
      text.scale = 1.5, # 文字标签的大小
      number.angles = 0, #柱子角度
      point.size = 3, #点的大小
      line.size = 1,
      main.bar.color = "#C5C9C6", #交集柱子颜色
      matrix.color = "#292D7C", #矩阵点的颜色
      mainbar.y.label = "eGene Intersections",
      sets.x.label = "cell types Intersections",

      queries = list(
      
  
      
                    list(query = intersects,params = list("B memory"),
                          color="#E41A1C",active = T),
                    list(query = intersects,params = list("B naive"),
                          color="#377EB8",active = T),
                    list(query = intersects,params = list("CD4 T memory"),
                          color="#4DAF4A",active = T),
          
                    list(query = intersects,params = list("CD4 T naive"),
                          color="#984EA3",active = T),
                    list(query = intersects,params = list("CD8 T effector"),
                          color="#FF7F00",active = T),
                    list(query = intersects,params = list("CD8 T naive"),
                          color="#FFFF33",active = T),
        
          
                    list(query = intersects,params = list("cDC"),
                          color="#A65628",active = T),
                    list(query = intersects,params = list("Monocyte classical"),
                          color="#F781BF",active = T),
      
                    list(query = intersects,params = list("Monocyte nonclassical"),
                          color="#999999",active = T),
                    list(query = intersects,params = list("NK"),
                          color="#FF34B3",active = T),
                    list(query = intersects,params = list("pDC"),
                          color="#BC8F8F",active = T),
                    list(query = intersects,params = list("Plasma"),
                          color="#20B2AA",active = T),
                    list(query = intersects,params = list("Platelet"),
                          color="#08519c",active = T))

     )
#dev.copy2pdf( file="/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_pdf/eGene.pdf",paper = "a4r")

inter_eGene <- get.venn.partitions(listInput_eGene_Tissue_T)
inter_eGene


for (i in 1:nrow(inter_eGene)) inter_eGene[i,'values'] <- paste(inter_eGene[[i,'..values..']], collapse = ', ')
write.table(inter_eGene[-c(14,15)],paste0("/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_pdf/","sub_All_eGene_inter.txt",sep=""), row.names = FALSE, sep = '\t', quote = FALSE)
inter_eGene

options(repr.plot.width = 16, repr.plot.height = 12)

#setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3","#BC8F8F","#20B2AA","#08519c")
   
setbarcolor<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c")

upset(fromList(listInput_eSNP_Tissue_T),
      sets.bar.color = setbarcolor, #柱子颜色
      nsets = 100,     # 绘制的最大集合个数
      nintersects = 50, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据集合的大小排序。
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE
      mb.ratio = c(0.7,0.3),   # 左侧和上方条形图的比例关系
      text.scale = 1.5, # 文字标签的大小
      number.angles = 0, #柱子角度
      point.size = 3, #点的大小
      line.size = 1,
      main.bar.color = "#C5C9C6", #交集柱子颜色
      matrix.color = "#292D7C", #矩阵点的颜色
      mainbar.y.label = "eGene Intersections",
      sets.x.label = "cell types Intersections",

      queries = list(
      
  
      
                    list(query = intersects,params = list("B memory"),
                          color="#E41A1C",active = T),
                    list(query = intersects,params = list("B naive"),
                          color="#377EB8",active = T),
                    list(query = intersects,params = list("CD4 T memory"),
                          color="#4DAF4A",active = T),
          
                    list(query = intersects,params = list("CD4 T naive"),
                          color="#984EA3",active = T),
                    list(query = intersects,params = list("CD8 T effector"),
                          color="#FF7F00",active = T),
                    list(query = intersects,params = list("CD8 T naive"),
                          color="#FFFF33",active = T),
        
          
                    list(query = intersects,params = list("cDC"),
                          color="#A65628",active = T),
                    list(query = intersects,params = list("Monocyte classical"),
                          color="#F781BF",active = T),
      
                    list(query = intersects,params = list("Monocyte nonclassical"),
                          color="#999999",active = T),
                    list(query = intersects,params = list("NK"),
                          color="#FF34B3",active = T),
                    list(query = intersects,params = list("pDC"),
                          color="#BC8F8F",active = T),
                    list(query = intersects,params = list("Plasma"),
                          color="#20B2AA",active = T),
                    list(query = intersects,params = list("Platelet"),
                          color="#08519c",active = T))

     )
#dev.copy2pdf( file="/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_pdf/eSNP.pdf",paper = "a4r")

inter_eSNP <- get.venn.partitions(listInput_eSNP_Tissue_T)
inter_eSNP

for (i in 1:nrow(inter_eSNP)) inter_eSNP[i,'values'] <- paste(inter_eSNP[[i,'..values..']], collapse = ', ')
write.table(inter_eSNP[-c(14,15)],paste0("/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_pdf/","sub_All_eSNP_inter.txt",sep=""), row.names = FALSE, sep = '\t', quote = FALSE)
inter_eSNP

options(repr.plot.width = 16, repr.plot.height = 12)

#setbarcolor<-c("#F2932A","#F5C040","#DB1E25","#E76124","#FF34B3","#BC8F8F","#20B2AA","#08519c")
   
setbarcolor<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c")

upset(fromList(listInput_eQTL_Tissue_T),
      sets.bar.color = setbarcolor, #柱子颜色
      nsets = 100,     # 绘制的最大集合个数
      nintersects = 50, #绘制的最大交集个数，NA则全部绘制
      order.by = "freq", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据集合的大小排序。
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE
      mb.ratio = c(0.7,0.3),   # 左侧和上方条形图的比例关系
      text.scale = 1.5, # 文字标签的大小
      number.angles = 0, #柱子角度
      point.size = 3, #点的大小
      line.size = 1,
      main.bar.color = "#C5C9C6", #交集柱子颜色
      matrix.color = "#292D7C", #矩阵点的颜色
      mainbar.y.label = "eGene Intersections",
      sets.x.label = "cell types Intersections",

      queries = list(
      
  
      
                    list(query = intersects,params = list("B memory"),
                          color="#E41A1C",active = T),
                    list(query = intersects,params = list("B naive"),
                          color="#377EB8",active = T),
                    list(query = intersects,params = list("CD4 T memory"),
                          color="#4DAF4A",active = T),
          
                    list(query = intersects,params = list("CD4 T naive"),
                          color="#984EA3",active = T),
                    list(query = intersects,params = list("CD8 T effector"),
                          color="#FF7F00",active = T),
                    list(query = intersects,params = list("CD8 T naive"),
                          color="#FFFF33",active = T),
        
          
                    list(query = intersects,params = list("cDC"),
                          color="#A65628",active = T),
                    list(query = intersects,params = list("Monocyte classical"),
                          color="#F781BF",active = T),
      
                    list(query = intersects,params = list("Monocyte nonclassical"),
                          color="#999999",active = T),
                    list(query = intersects,params = list("NK"),
                          color="#FF34B3",active = T),
                    list(query = intersects,params = list("pDC"),
                          color="#BC8F8F",active = T),
                    list(query = intersects,params = list("Plasma"),
                          color="#20B2AA",active = T),
                    list(query = intersects,params = list("Platelet"),
                          color="#08519c",active = T))

     )
#dev.copy2pdf( file="/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_pdf/eQTL.pdf",paper = "a4r")

inter_eQTL <- get.venn.partitions(listInput_eQTL_Tissue_T)
inter_eQTL


for (i in 1:nrow(inter_eQTL))inter_eQTL[i,'values'] <- paste(inter_eQTL[[i,'..values..']], collapse = ', ')
write.table(inter_eQTL[-c(14,15)],paste0("/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_pdf/","sub_All_eQTL_inter.txt",sep=""), row.names = FALSE, sep = '\t', quote = FALSE)
inter_eQTL

All_exp_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/1EQTL_exp/subtype_exp/All_subtype_exp/txt/"
subtype_exp_ststus<-list.files(All_exp_path)
All_exp_path
subtype_exp_ststus

All_exp_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/1EQTL_exp/subtype_exp/All_subtype_exp/0_txt/"
subtype_exp_ststus<-list.files(All_exp_path)

B_memory              <- read.table( paste0(All_exp_path,subtype_exp_ststus[1]),header=T,row.names = 1)
B_naive               <- read.table( paste0(All_exp_path,subtype_exp_ststus[2]),header=T,row.names = 1)
CD4_T_memory          <- read.table( paste0(All_exp_path,subtype_exp_ststus[3]),header=T,row.names = 1)
CD4_T_naive           <- read.table( paste0(All_exp_path,subtype_exp_ststus[4]),header=T,row.names = 1)
CD8_T_effector        <- read.table( paste0(All_exp_path,subtype_exp_ststus[5]),header=T,row.names = 1)
CD8_T_naive           <- read.table( paste0(All_exp_path,subtype_exp_ststus[6]),header=T,row.names = 1)
cDC                   <- read.table( paste0(All_exp_path,subtype_exp_ststus[7]),header=T,row.names = 1)
Monocyte_classical    <- read.table( paste0(All_exp_path,subtype_exp_ststus[8]),header=T,row.names = 1)
Monocyte_nonclassical <- read.table( paste0(All_exp_path,subtype_exp_ststus[9]),header=T,row.names = 1)
NK                    <- read.table( paste0(All_exp_path,subtype_exp_ststus[10]),header=T,row.names = 1)
pDC                   <- read.table( paste0(All_exp_path,subtype_exp_ststus[11]),header=T,row.names = 1)
Plasma                <- read.table( paste0(All_exp_path,subtype_exp_ststus[12]),header=T,row.names = 1)
Platelet              <- read.table( paste0(All_exp_path,subtype_exp_ststus[13]),header=T,row.names = 1)


B_memory              <- read.table( paste0(All_exp_path,subtype_exp_ststus[1]),header=T,row.names = 1)
B_memory

All_subtype_genotype_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/22+X/need/"
subtype_all <- list.files(All_subtype_genotype_path)
All_subtype_genotype_path
subtype_all

All_subtype_genotype_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/22+X/need/"
subtype_all <- list.files(All_subtype_genotype_path)

All_genotype_transpose<-fread(paste0(All_subtype_genotype_path,subtype_all[2],sep=""),
                                                 sep=" ",data.table=F)
#asymptomatic_genotype_transpose<-column_to_rownames(asymptomatic_genotype_transpose,var="V1")
#colnames(asymptomatic_genotype_transpose)[1]<-"ID"

All_genotype_transpose

B_memory             
B_naive               
CD4_T_memory          
CD4_T_naive           
CD8_T_effector        
CD8_T_naive           
cDC                  
Monocyte_classical    
Monocyte_nonclassical 
NK                  
pDC                  
Plasma               
Platelet              


eQTL_expression<-function(Gene,SNP){

    ###基因型数据###
    SNP_genotype <- All_genotype_transpose[which(All_genotype_transpose$ID == SNP),] #定位SNP
    SNP_genotype<- SNP_genotype[,2:dim(SNP_genotype)[2]] #提取需要样本的行数
    SNP_genotype <- t(as.data.frame(SNP_genotype))
    
    
    
    ###表达数据###  
    B_memory_gene_exp <- B_memory[which(rownames(B_memory) == Gene),] #定位基因表达信息
    B_memory_gene_exp <- rbind(B_memory_gene_exp,"B_memory")
    B_memory_gene_exp <- t(B_memory_gene_exp)
    B_memory_gene_exp <- as.data.frame(B_memory_gene_exp)
    
    B_naive_gene_exp <- B_naive[which(rownames(B_naive) == Gene),] #定位基因表达信息
    B_naive_gene_exp <- rbind(B_naive_gene_exp,"B_naive")
    B_naive_gene_exp <- t(B_naive_gene_exp)
    B_naive_gene_exp <- as.data.frame(B_naive_gene_exp)    

    CD4_T_memory_gene_exp <- CD4_T_memory[which(rownames(CD4_T_memory) == Gene),] #定位基因表达信息
    CD4_T_memory_gene_exp <- rbind(CD4_T_memory_gene_exp,"CD4_T_memory")
    CD4_T_memory_gene_exp <- t(CD4_T_memory_gene_exp)
    CD4_T_memory_gene_exp <- as.data.frame(CD4_T_memory_gene_exp)
    
    CD4_T_naive_gene_exp <- CD4_T_naive[which(rownames(CD4_T_naive) == Gene),] #定位基因表达信息
    CD4_T_naive_gene_exp <- rbind(CD4_T_naive_gene_exp,"CD4_T_naive")
    CD4_T_naive_gene_exp <- t(CD4_T_naive_gene_exp)
    CD4_T_naive_gene_exp <- as.data.frame(CD4_T_naive_gene_exp)
      
    CD8_T_effector_gene_exp <- CD8_T_effector[which(rownames(CD8_T_effector) == Gene),] #定位基因表达信息
    CD8_T_effector_gene_exp <- rbind(CD8_T_effector_gene_exp,"CD8_T_effector")
    CD8_T_effector_gene_exp <- t(CD8_T_effector_gene_exp)
    CD8_T_effector_gene_exp <- as.data.frame(CD8_T_effector_gene_exp)
    
    CD8_T_naive_gene_exp <- CD8_T_naive[which(rownames(CD8_T_naive) == Gene),] #定位基因表达信息
    CD8_T_naive_gene_exp <- rbind(CD8_T_naive_gene_exp,"CD8_T_naive")
    CD8_T_naive_gene_exp <- t(CD8_T_naive_gene_exp)
    CD8_T_naive_gene_exp <- as.data.frame(CD8_T_naive_gene_exp)
     

    cDC_gene_exp <- cDC[which(rownames(cDC) == Gene),] #定位基因表达信息
    cDC_gene_exp <- rbind(cDC_gene_exp,"cDC")
    cDC_gene_exp <- t(cDC_gene_exp)
    cDC_gene_exp <- as.data.frame(cDC_gene_exp)
    
    
    
    Monocyte_classical_gene_exp <- Monocyte_classical[which(rownames(Monocyte_classical) == Gene),] #定位基因表达信息
    Monocyte_classical_gene_exp <- rbind(Monocyte_classical_gene_exp,"Monocyte_classical")
    Monocyte_classical_gene_exp <- t(Monocyte_classical_gene_exp)
    Monocyte_classical_gene_exp <- as.data.frame(Monocyte_classical_gene_exp)
    
    Monocyte_nonclassical_gene_exp <- Monocyte_nonclassical[which(rownames(Monocyte_nonclassical) == Gene),] #定位基因表达信息
    Monocyte_nonclassical_gene_exp <- rbind(Monocyte_nonclassical_gene_exp,"Monocyte_nonclassical")
    Monocyte_nonclassical_gene_exp <- t(Monocyte_nonclassical_gene_exp)
    Monocyte_nonclassical_gene_exp <- as.data.frame(Monocyte_nonclassical_gene_exp)
       
    
    NK_gene_exp <- NK[which(rownames(NK) == Gene),] #定位基因表达信息
    NK_gene_exp <- rbind(NK_gene_exp,"NK")
    NK_gene_exp <- t(NK_gene_exp)
    NK_gene_exp <- as.data.frame(NK_gene_exp)
    
    pDC_gene_exp <- pDC[which(rownames(pDC) == Gene),] #定位基因表达信息
    pDC_gene_exp <- rbind(pDC_gene_exp,"pDC")
    pDC_gene_exp <- t(pDC_gene_exp)
    pDC_gene_exp <- as.data.frame(pDC_gene_exp)
    
    
    Plasma_gene_exp <- Plasma[which(rownames(Plasma) == Gene),] #定位基因表达信息
    Plasma_gene_exp <- rbind(Plasma_gene_exp,"Plasma")
    Plasma_gene_exp <- t(Plasma_gene_exp)
    Plasma_gene_exp <- as.data.frame(Plasma_gene_exp)
    
    
    Platelet_gene_exp <- Platelet[which(rownames(Platelet) == Gene),] #定位基因表达信息
    Platelet_gene_exp <- rbind(Platelet_gene_exp,"Platelet")
    Platelet_gene_exp <- t(Platelet_gene_exp)
    Platelet_gene_exp <- as.data.frame(Platelet_gene_exp)
        
     


    ###构建eqtl表达数据集###
        
    B_memory_eqtl_expression_data        <- cbind(B_memory_gene_exp,SNP_genotype)
    colnames(B_memory_eqtl_expression_data)<-c("expression","subtype","genotype")
    
    B_naive_eqtl_expression_data        <- cbind(B_naive_gene_exp,SNP_genotype)
    colnames(B_naive_eqtl_expression_data)<-c("expression","subtype","genotype")
    
    
    CD4_T_memory_eqtl_expression_data      <- cbind(CD4_T_memory_gene_exp,SNP_genotype)
    colnames(CD4_T_memory_eqtl_expression_data)<-c("expression","subtype","genotype")
    
    
    CD4_T_naive_eqtl_expression_data      <- cbind(CD4_T_naive_gene_exp,SNP_genotype)
    colnames(CD4_T_naive_eqtl_expression_data)<-c("expression","subtype","genotype")
    
    
    CD8_T_effector_eqtl_expression_data      <- cbind(CD8_T_effector_gene_exp,SNP_genotype)
    colnames(CD8_T_effector_eqtl_expression_data)<-c("expression","subtype","genotype")
    
    
    CD8_T_naive_eqtl_expression_data      <- cbind(CD8_T_naive_gene_exp,SNP_genotype)
    colnames(CD8_T_naive_eqtl_expression_data)<-c("expression","subtype","genotype")
    
    
    cDC_eqtl_expression_data      <- cbind(cDC_gene_exp,SNP_genotype)
    colnames(cDC_eqtl_expression_data)<-c("expression","subtype","genotype")    
    
    Monocyte_classical_eqtl_expression_data <- cbind(Monocyte_classical_gene_exp,SNP_genotype)
    colnames(Monocyte_classical_eqtl_expression_data)<-c("expression","subtype","genotype") 
    
    Monocyte_nonclassical_eqtl_expression_data <- cbind(Monocyte_nonclassical_gene_exp,SNP_genotype)
    colnames(Monocyte_nonclassical_eqtl_expression_data)<-c("expression","subtype","genotype") 
    
    
    
    NK_eqtl_expression_data       <- cbind(NK_gene_exp,SNP_genotype)
    colnames(NK_eqtl_expression_data)<-c("expression","subtype","genotype") 
    
    pDC_eqtl_expression_data      <- cbind(pDC_gene_exp,SNP_genotype)
    colnames(pDC_eqtl_expression_data)<-c("expression","subtype","genotype") 
    
    Plasma_eqtl_expression_data      <- cbind(Plasma_gene_exp,SNP_genotype)
    colnames(Plasma_eqtl_expression_data)<-c("expression","subtype","genotype") 
    
    
    Platelet_eqtl_expression_data <- cbind(Platelet_gene_exp,SNP_genotype)
    colnames(Platelet_eqtl_expression_data)<-c("expression","subtype","genotype") 
    

    
    allsubtype_expression_data<-rbind(B_memory_eqtl_expression_data,
                                       B_naive_eqtl_expression_data,
                                       CD4_T_memory_eqtl_expression_data,
                                       CD4_T_naive_eqtl_expression_data,
                                       CD8_T_effector_eqtl_expression_data,
                                       CD8_T_naive_eqtl_expression_data,
                                       cDC_eqtl_expression_data,
                                       Monocyte_classical_eqtl_expression_data ,
                                       Monocyte_nonclassical_eqtl_expression_data,
                                       NK_eqtl_expression_data,
                                       pDC_eqtl_expression_data,
                                       Plasma_eqtl_expression_data,
                                       Platelet_eqtl_expression_data )
 
    allsubtype_expression_data$expression<-as.numeric(allsubtype_expression_data$expression)
    allsubtype_expression_data$genotype<-as.character(allsubtype_expression_data$genotype)
    allsubtype_expression_data<-as.data.frame(allsubtype_expression_data)
    return(allsubtype_expression_data)
}  

Gene<-c("RPS26")
SNP<-c("rs705704")
expression_data<-eQTL_expression(Gene = Gene,SNP = SNP)
expression_data

options(repr.plot.width = 20, repr.plot.height = 5)

expression_data$subtype<-factor(expression_data$subtype,
                                     levels=c(
                                             "B_memory",
                                             "B_naive",               
                                             "CD4_T_memory",          
                                             "CD4_T_naive" ,          
                                             "CD8_T_effector",        
                                             "CD8_T_naive",           
                                             "cDC",                  
                                             "Monocyte_classical",    
                                             "Monocyte_nonclassical", 
                                             "NK",                  
                                             "pDC",                  
                                             "Plasma" ,              
                                             "Platelet"  
                                             ))


    ggplot(data = expression_data, aes(x = genotype, y = expression)) +
    geom_violin(aes(colour = subtype),trim=FALSE,width=0.8,size=1)+ 
    geom_boxplot(aes(colour = subtype),width=0.1,size = 0.5)+   
    geom_jitter(aes(colour = subtype ),width=0.1,shape=21,size=0.8) +
    scale_colour_manual(values = c(
    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF',
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c"
))+
    theme_test() +
    labs(title=c(paste0(Gene," / ",SNP)),
        x = "Genotype",
        y = c(paste0("Expression of ",Gene)))+
    facet_grid(. ~ subtype)+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=14),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####去除图例###
        legend.position = 'none',
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill="#E8DAEF", size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1) #图形线框

     )

Gene<-c("SNHG5")
SNP<-c("rs1059307")
expression_data<-eQTL_expression(Gene = Gene,SNP = SNP)
expression_data

options(repr.plot.width = 20, repr.plot.height = 5)

expression_data$subtype<-factor(expression_data$subtype,
                                     levels=c(
                                             "B_memory",
                                             "B_naive",               
                                             "CD4_T_memory",          
                                             "CD4_T_naive" ,          
                                             "CD8_T_effector",        
                                             "CD8_T_naive",           
                                             "cDC",                  
                                             "Monocyte_classical",    
                                             "Monocyte_nonclassical", 
                                             "NK",                  
                                             "pDC",                  
                                             "Plasma" ,              
                                             "Platelet"  
                                             ))


    ggplot(data = expression_data, aes(x = genotype, y = expression)) +
    geom_violin(aes(colour = subtype),trim=FALSE,width=0.8,size=1)+ 
    geom_boxplot(aes(colour = subtype),width=0.1,size = 0.5)+   
    geom_jitter(aes(colour = subtype ),width=0.1,shape=21,size=0.8) +
    scale_colour_manual(values = c(
    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF',
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c"
))+
    theme_test() +
    labs(title=c(paste0(Gene," / ",SNP)),
        x = "Genotype",
        y = c(paste0("Expression of ",Gene)))+
    facet_grid(. ~ subtype)+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=14),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####去除图例###
        legend.position = 'none',
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill="#E8DAEF", size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1) #图形线框

     )

Gene<-c("SNHG5")
SNP<-c("rs3087978")
expression_data<-eQTL_expression(Gene = Gene,SNP = SNP)
expression_data

options(repr.plot.width = 20, repr.plot.height = 5)

expression_data$subtype<-factor(expression_data$subtype,
                                     levels=c(
                                             "B_memory",
                                             "B_naive",               
                                             "CD4_T_memory",          
                                             "CD4_T_naive" ,          
                                             "CD8_T_effector",        
                                             "CD8_T_naive",           
                                             "cDC",                  
                                             "Monocyte_classical",    
                                             "Monocyte_nonclassical", 
                                             "NK",                  
                                             "pDC",                  
                                             "Plasma" ,              
                                             "Platelet"  
                                             ))


    ggplot(data = expression_data, aes(x = genotype, y = expression)) +
    geom_violin(aes(colour = subtype),trim=FALSE,width=0.8,size=1)+ 
    geom_boxplot(aes(colour = subtype),width=0.1,size = 0.5)+   
    geom_jitter(aes(colour = subtype ),width=0.1,shape=21,size=0.8) +
    scale_colour_manual(values = c(
    '#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF',
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c"
))+
    theme_test() +
    labs(title=c(paste0(Gene," / ",SNP)),
        x = "Genotype",
        y = c(paste0("Expression of ",Gene)))+
    facet_grid(. ~ subtype)+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=14),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####去除图例###
        legend.position = 'none',
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill="#E8DAEF", size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1) #图形线框

     )

allfiles_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_sub/"
files<-list.files(allfiles_path)

B_memory_all              <-paste0(allfiles_path,files[1])
B_naive_all               <-paste0(allfiles_path,files[2])
CD4_T_memory_all          <-paste0(allfiles_path,files[3])
CD4_T_naive_all           <-paste0(allfiles_path,files[4])
CD8_T_effector_all        <-paste0(allfiles_path,files[5])
CD8_T_naive_all           <-paste0(allfiles_path,files[6])
cDC_all                   <-paste0(allfiles_path,files[7])
Monocyte_classical_all    <-paste0(allfiles_path,files[8])
Monocyte_nonclassical_all <-paste0(allfiles_path,files[9])
NK_all                    <-paste0(allfiles_path,files[10])
pDC_all                   <-paste0(allfiles_path,files[11])
Plasma_all                <-paste0(allfiles_path,files[12])
Platelet_all              <-paste0(allfiles_path,files[13])


allfiles <-c(B_memory_all,B_naive_all,CD4_T_memory_all,CD4_T_naive_all,CD8_T_effector_all,CD8_T_naive_all,cDC_all,
             Monocyte_classical_all,Monocyte_nonclassical_all,NK_all,pDC_all,Plasma_all,Platelet_all)
allfiles
allfiles_path
files


             

All_eGene_num_subtype<-c()
for(i in 1:length(allfiles)){
    eqtl_file<-read.table(allfiles[i],header=T,sep="\t")
    eGene_num<-length(unique(eqtl_file$gene))
    subtype<-strsplit(strsplit(allfiles[i],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
    subtype_eGene_num<-cbind(subtype,eGene_num)
    All_eGene_num_subtype<-rbind(All_eGene_num_subtype,subtype_eGene_num)
    #All_eGene_num_type<-as.data.frame(cbind(c(celltype),All_eGene_num_type))
    print(i)
}

#colnames(All_eGene_num_type)<-c("type","subtype","eGene_num")
colnames(All_eGene_num_subtype)<-c("subtype","eGene_num")
All_result_eGene_num<-as.data.frame(All_eGene_num_subtype)
All_result_eGene_num

#All_result_eGene_num$sample_Freq<-53
#All_result_eGene_num

All_subtype_freq<-matrix(ncol = 2,nrow = 13)
colnames(All_subtype_freq)<-c("subtype","sample_Freq")
All_subtype_freq<-as.data.frame(All_subtype_freq)

All_subtype_freq$subtype<-All_result_eGene_num$subtype
#All_subtype_freq$sample_Freq<-53
All_subtype_freq$sample_Freq[1:13]<-c(47,53,53,53,53,53,53,53,53,53,53,44,53)
All_subtype_freq







Cell_eGene_num<-cbind(All_subtype_freq,All_result_eGene_num$eGene_num)
colnames(Cell_eGene_num)<-c("subtype","sample_Freq","eGene_num")
Cell_eGene_num<-as.data.frame(Cell_eGene_num)
Cell_eGene_num$sample_Freq<-as.numeric(Cell_eGene_num$sample_Freq)
Cell_eGene_num$eGene_num<-as.numeric(Cell_eGene_num$eGene_num)
Cell_eGene_num

options(repr.plot.width = 10, repr.plot.height = 8)
ggplot(data=Cell_eGene_num,aes(x=sample_Freq,y=eGene_num))+geom_point(aes(colour=subtype),size=5)+
geom_text_repel(aes(x=sample_Freq,y=eGene_num, label = subtype))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c"

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )

allfiles_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_sub/"
files<-list.files(allfiles_path)

B_memory_all              <-paste0(allfiles_path,files[1])
B_naive_all               <-paste0(allfiles_path,files[2])
CD4_T_memory_all          <-paste0(allfiles_path,files[3])
CD4_T_naive_all           <-paste0(allfiles_path,files[4])
CD8_T_effector_all        <-paste0(allfiles_path,files[5])
CD8_T_naive_all           <-paste0(allfiles_path,files[6])
cDC_all                   <-paste0(allfiles_path,files[7])
Monocyte_classical_all    <-paste0(allfiles_path,files[8])
Monocyte_nonclassical_all <-paste0(allfiles_path,files[9])
NK_all                    <-paste0(allfiles_path,files[10])
pDC_all                   <-paste0(allfiles_path,files[11])
Plasma_all                <-paste0(allfiles_path,files[12])
Platelet_all              <-paste0(allfiles_path,files[13])


allfiles <-c(B_memory_all,B_naive_all,CD4_T_memory_all,CD4_T_naive_all,CD8_T_effector_all,CD8_T_naive_all,cDC_all,
             Monocyte_classical_all,Monocyte_nonclassical_all,NK_all,pDC_all,Plasma_all,Platelet_all)
allfiles
allfiles_path
files


All_esnp_num_subtype<-c()
for(i in 1:length(allfiles)){
    eqtl_file<-read.table(allfiles[i],header=T,sep="\t")
    esnp_num<-length(unique(eqtl_file$snps))
    subtype<-strsplit(strsplit(allfiles[i],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
    subtype_esnp_num<-cbind(subtype,esnp_num)
    All_esnp_num_subtype<-rbind(All_esnp_num_subtype,subtype_esnp_num)
    #All_eGene_num_type<-as.data.frame(cbind(c(celltype),All_eGene_num_type))
    print(i)
}

#colnames(All_eGene_num_type)<-c("type","subtype","eGene_num")
colnames(All_esnp_num_subtype)<-c("subtype","esnp_num")
All_esnp_num_subtype<-as.data.frame(All_esnp_num_subtype)
All_esnp_num_subtype


All_subtype_freq<-matrix(ncol = 2,nrow = 13)
colnames(All_subtype_freq)<-c("subtype","sample_Freq")
All_subtype_freq<-as.data.frame(All_subtype_freq)

All_subtype_freq$subtype<-All_result_eGene_num$subtype
#All_subtype_freq$sample_Freq<-53
All_subtype_freq$sample_Freq[1:13]<-c(47,53,53,53,53,53,53,53,53,53,53,44,53)
All_subtype_freq

Cell_esnp_num<-cbind(All_subtype_freq,All_result_esnp_num$esnp_num)
colnames(Cell_esnp_num)<-c("subtype","sample_Freq","esnp_num")
Cell_esnp_num<-as.data.frame(Cell_esnp_num)
Cell_esnp_num$sample_Freq<-as.numeric(Cell_esnp_num$sample_Freq)
Cell_esnp_num$esnp_num<-as.numeric(Cell_esnp_num$esnp_num)
Cell_esnp_num

 options(repr.plot.width = 10, repr.plot.height = 8)
ggplot(data=Cell_esnp_num,aes(x=sample_Freq,y=esnp_num))+geom_point(aes(colour=subtype),size=5)+
geom_text_repel(aes(x=sample_Freq,y=esnp_num, label = subtype))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

      "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c"

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )

#sub细胞num
all<-readRDS("/data_alluser/CXY/keti/singlecell_eqtl/result/RDS/all_subtype_annotated.rds")
all@meta.data
z<-all@meta.data[,c(5,10,11)]
z
table(c(z[which(z$subtype=="CD4 T memory"),])$Status)

sub_All_eQTL_num_ststus_AllCell<-fread("/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/2EQTL_result/sub_All_eQTL_num_ststus_AllCell.txt")
sub_All_eQTL_num_ststus_AllCell

options(repr.plot.width = 10, repr.plot.height = 8)
ggplot(data=sub_All_eQTL_num_ststus_AllCell,aes(x=status_cellnum,y=eqtl_num))+geom_point(aes(colour=sub),size=5)+
geom_text_repel(aes(x=status_cellnum,y=eqtl_num, label = status))+
theme_classic(base_size = 16)+
scale_colour_manual(values=c(

"#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
               "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c"

))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=15,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 15,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 15,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_text(face="bold",color="black",size=11),
         axis.text.y = element_text(face="bold",color="black",size=11),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.3),
        axis.ticks.y = element_line(color="black",size=0.3),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
                ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 10)

     )









options(repr.plot.width = 12, repr.plot.height = 5)
ggplot(data=Cell_eGene_num,aes(x=subtype,y=eGene_num,fill=subtype))+
geom_bar(stat="identity",width=0.5)+facet_grid(. ~ sample_Freq)+
theme_bw()+
scale_fill_manual(values=c("#C5CC62","#388246","#117A82","#1B66B0","#C57FCC","#892E73",
                           
                           "#F5C040","#F2932A","#E76124","#DB1E25",
                            "#E41A1C","#377EB8","#4DAF4A"
        #"#984EA3","#FF7F00","#FFFF33","#A65628", "#F781BF","#999999","#FF34B3","#BC8F8F","#20B2AA","#08519c"
                          ))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_blank(),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill=c("#EA9390","#682342"), size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1), #图形线框
        
        ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 12)
     )










options(repr.plot.width = 12, repr.plot.height = 5)
ggplot(data=Cell_esnp_num,aes(x=subtype,y=esnp_num,fill=subtype))+
geom_bar(stat="identity",width=0.5)+facet_grid(. ~ sample_Freq)+
theme_bw()+
scale_fill_manual(values=c("#C5CC62","#388246","#117A82","#1B66B0","#C57FCC","#892E73",
                           
                           "#F5C040","#F2932A","#E76124","#DB1E25",
                            "#E41A1C","#377EB8","#4DAF4A"
                          ))+
    theme(
        ####坐标轴标签####
        plot.title=element_text(face="bold.italic", #字体
                                color="black",      #颜色
                                size=24,            #大小
                                hjust=0.5,            #位置
                                vjust=0.5),
         axis.title.x = element_text(face="bold", 
                                size = 24,
                                vjust= 0),
         axis.title.y = element_text(face="bold",
                                size = 24,
                                vjust=1),
        ####坐标轴刻度####
         axis.text.x = element_blank(),
         axis.text.y = element_text(face="bold",color="black",size=14),
       
        ####坐标轴刻度线段####
        axis.ticks.x = element_line(color="black",size=0.5),
        axis.ticks.y = element_line(color="black",size=0.5),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        
        ####分面标签、背景####
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background.x = element_rect(color="black", fill=c("#EA9390","#682342"), size=1, linetype="solid"),
        panel.border = element_rect(color = 'black',size = 1), #图形线框
        
        ####图例####
        legend.title=element_text(face="bold",size = 15),
        legend.text=element_text(face="bold",size = 12)
     )


allfiles_path<-"/data_alluser/CXY/keti/singlecell_eqtl/eqtl/result/3EQTL_upset/All/upset_sub/"
files<-list.files(allfiles_path)


B_memory_all              <-paste0(allfiles_path,files[1])
B_naive_all               <-paste0(allfiles_path,files[2])
CD4_T_memory_all          <-paste0(allfiles_path,files[3])
CD4_T_naive_all           <-paste0(allfiles_path,files[4])
CD8_T_effector_all        <-paste0(allfiles_path,files[5])
CD8_T_naive_all           <-paste0(allfiles_path,files[6])
cDC_all                   <-paste0(allfiles_path,files[7])
Monocyte_classical_all    <-paste0(allfiles_path,files[8])
Monocyte_nonclassical_all <-paste0(allfiles_path,files[9])
NK_all                    <-paste0(allfiles_path,files[10])
pDC_all                   <-paste0(allfiles_path,files[11])
Plasma_all                <-paste0(allfiles_path,files[12])
Platelet_all              <-paste0(allfiles_path,files[13])

allfiles <-c(B_memory_all,B_naive_all,CD4_T_memory_all,CD4_T_naive_all,CD8_T_effector_all,CD8_T_naive_all,cDC_all,
             Monocyte_classical_all,Monocyte_nonclassical_all,NK_all,pDC_all,Plasma_all,Platelet_all)
allfiles
allfiles_path
files


 result_all_eGene<-c()
for(i in 1:length(allfiles)){
     subtype<-strsplit(strsplit(allfiles[i],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
    file<-read.table(allfiles[i],header=T,sep="\t")
    egene_file<-unique(file$gene)

    overlapped_result<-c()
    
    for(j in 1:length(allfiles)){
        subtype_2<-strsplit(strsplit(allfiles[j],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
        file_2<-read.table(allfiles[j],header=T,sep="\t")
        egene_file_2<-unique(file_2$gene)
        intersect_egene<-intersect(egene_file, egene_file_2)
        overlapped_pt<-length(intersect_egene)/length(egene_file)
        overlapped_result<-rbind(overlapped_result,overlapped_pt)
        rownames(overlapped_result)[j]<-subtype_2
    }
    
    result_all_eGene<-cbind(result_all_eGene,overlapped_result)
    
    colnames(result_all_eGene)[i]<-subtype
    
}
result_all_eGene

#函数使用：
#pheatmap(
  mat = test, # 热图的数据源，要保证为数值矩阵，类型为numeric
  scale = "none", # 数值标准化，可以规定按行（row）或按列（column）
  cluster_rows = TRUE, #  是否按行聚类
  cluster_cols = TRUE, #  是否按列聚类
  legend = TRUE, #  图例是否显示
  legend_breaks = NA, # 图例是否断点标注
  legend_labels = NA, # 图例断点标注的标题
  show_rownames = TRUE, # 是否显示行名
  show_colnames = TRUE, # 是否显示列名
  main = NA, #  热图标题
  fontsize = 10, #  热图字体大小
  fontsize_row = 10, #  热图行名字体大小，默认为fontsize
  fontsize_col = 10, #  热图列名字体大小，默认为fontsize
  angle_col = 0, #  热图列名角度
  display_numbers = FALSE, #  矩阵的数值是否显示在热图上
  number_format = "%.2f", # 单元格中数值的显示方式
  fontsize_number = 8, #  显示在热图上的矩阵数值的大小，默认为0.8*fontsize
  gaps_row = NULL, #  在行未聚类时使用，确定将间隙放置在热图中的位置。
  gaps_col = NULL, #  在列未聚类时使用，确定将间隙放置在热图中的位置。
  labels_row = NULL, #  使用行标签代替行名
  labels_col = NULL, #  使用列标签代替列名
  filename = NA, #  保存路径和文件名
  width = NA, # 保存图片的宽度
  height = NA, #  保存图片的高度
  na_col = "#DDDDDD" #  对“NA”值对应的单元格填充颜色
)

colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)

plot_size(10,8)
pheatmap(result_all_eGene,cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c( "#FFFFFF", "#0070B0"))(50),fontsize_row = 10,
         display_numbers = T, #  矩阵的数值是否显示在热图上
  number_format = "%.2f", # 单元格中数值的显示方式
  fontsize_number = 10, #  显示在热图上的矩阵数值的大小，默认为0.8*fontsize
          gaps_row = NULL, #  在行未聚类时使用，确定将间隙放置在热图中的位置。
  gaps_col = NULL, #  在列未聚类时使用，确定将间隙放置在热图中的位置。
         filename = NA, #  保存路径和文件名
  width = NA, # 保存图片的宽度
  height = NA, #  保存图片的高度
         cellwidth=40,#热图单位元素（cell）的宽度，NA表示依据窗口自动调整
        
         border_color = "#B1B4B6"#热图的单位元素的描边颜色，NA表示不描边默认：“gray60”
)

result_all_SNP<-c()
for(i in 1:length(allfiles)){
     subtype<-strsplit(strsplit(allfiles[i],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
    file<-read.table(allfiles[i],header=T,sep="\t")
    snps_file<-unique(file$snps)

    overlapped_result<-c()
    
    for(j in 1:length(allfiles)){
        subtype_2<-strsplit(strsplit(allfiles[j],split="_cis_eqtls")[[1]][1],split="/")[[1]][11]
        file_2<-read.table(allfiles[j],header=T,sep="\t")
        snps_file_2<-unique(file_2$snps)
        intersect_snps<-intersect(snps_file, snps_file_2)
        overlapped_pt<-length(intersect_snps)/length(snps_file)
        overlapped_result<-rbind(overlapped_result,overlapped_pt)
        rownames(overlapped_result)[j]<-subtype_2
    }
    
    result_all_SNP<-cbind(result_all_SNP,overlapped_result)
    
    colnames(result_all_SNP)[i]<-subtype
    
}
result_all_SNP

plot_size(10,8)
pheatmap(result_all_SNP,cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c( "#FFFFFF", "#0070B0"))(50),fontsize_row = 10,
          display_numbers = T, #  矩阵的数值是否显示在热图上
           number_format = "%.2f", # 单元格中数值的显示方式
  fontsize_number = 10, #  显示在热图上的矩阵数值的大小，默认为0.8*fontsize
          gaps_row = NULL, #  在行未聚类时使用，确定将间隙放置在热图中的位置。
  gaps_col = NULL, #  在列未聚类时使用，确定将间隙放置在热图中的位置。
         filename = NA, #  保存路径和文件名
  width = NA, # 保存图片的宽度
  height = NA, #  保存图片的高度
         cellwidth=40,#热图单位元素（cell）的宽度，NA表示依据窗口自动调整
        
         border_color = "#B1B4B6")














