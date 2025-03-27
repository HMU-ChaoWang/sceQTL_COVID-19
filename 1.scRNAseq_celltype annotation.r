library(Seurat)
library(ggplot2)
library(harmony)
library(ggsci)
#library(DoubletFinder)
library(dplyr)
library(patchwork)
library(tidyverse)
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
library(hdf5r)
library(COSG)
library(RColorBrewer) 
library(ggsci)
library(colorspace)
library(viridis)
library(scales)

plot_size <- function(x, y) {
    options(repr.plot.width = x, repr.plot.height = y)
}

path = "path"
files <- list.files(path)
dir = paste0(path,'/',files)
names(dir) <- files

files

sample_list <- list()
for(i in files){
    spleen.data <- Matrix::readMM(file = paste0(path,'/',i,'/matrix.mtx.gz'))
    
    barcode <- fread(paste0(path,'/',i,'/barcodes.tsv.gz'),header=F)
    barcode$V1 <- paste0(i,'_',barcode$V1)
                     
    features <- fread(paste0(path,'/',i,'/features.tsv.gz'),header=F)
    
    colnames(spleen.data) = barcode$V1
    rownames(spleen.data) = features$V2
    print(i)
    sample_list[[i]]<- (spleen.data,
                                          project=i,
                                          min.cells=3,
                                          min.features=200)
    
}

 spleen.data <- Matrix::readMM('/path/matrix.mtx.gz')
    spleen.data 






sample_list <- list()
for(i in files){
    spleen.data <- Matrix::readMM('/path/matrix.mtx.gz')
    
    barcode <- fread('/path/barcodes.tsv.gz',header=F)
   
                     
    features <- fread('/path/features.tsv.gz',header=F)
    
    colnames(spleen.data) = barcode$V1
    rownames(spleen.data) = features$V2
    
    sample_list<- CreateSeuratObject(spleen.data,
                                          project="GSE165080",
                                          min.cells=3,
                                          min.features=200)
    
}

sample_all <- merge(sample_list[[1]],sample_list[2:length(sample_list)])

sample_all

plot_size(15,7)
sample_all[['percent.MT']] <- PercentageFeatureSet(sample_all,pattern='^MT')
feats <- c('nFeature_RNA','nCount_RNA',"percent.MT")

sample_all@meta.data

VlnPlot(sample_all,feature=feats,group.by='orig.ident',pt.size=0.1,ncol=3)+
   NoLegend()

sample_all <- subset(sample_all,subset= nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 20000 &
                      percent.MT < 30)

VlnPlot(sample_all, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()

sample_all

plot_size(15,8)
FeatureScatter(sample_all,'nFeature_RNA','nCount_RNA')+
    geom_vline(aes(xintercept=200),colour='red',linetype='dashed')+
   facet_wrap(~sample_all@meta.data$orig.ident)

meta_annation = c(
        'GSM2025728'='asymptomatic case',
        'GSM2025729'='asymptomatic case',
        'GSM2025730'='asymptomatic case',
        'GSM2025731'='asymptomatic case',
        'GSM2025732'='asymptomatic case',
        'GSM2025733'='Moderate',
        'GSM2025734'='Moderate',
        'GSM2025735'='Moderate',
        'GSM2025736'='Moderate',
        'GSM2025737'='Moderate',
        'GSM2025738'='Moderate',
        'GSM2025739'='Moderate',
        'GSM2025740'='Moderate',
        'GSM2025741'='Moderate',
        'GSM2025742'='Moderate',
        'GSM2025743'='Moderate',
        'GSM2025744'='Moderate',
        'GSM2025745'='Moderate',
        'GSM2025746'='Severe',
        'GSM2025747'='Severe',
        'GSM2025748'='Severe',
        'GSM2025749'='Severe',
        'GSM2025750'='Severe',
        'GSM2025751'='Severe',
        'GSM2025752'='Severe',
        'GSM2025753'='Severe',
        'GSM2025754'='Severe',
        'GSM2025755'='Severe',
        'GSM2025756'='Severe',
        'GSM2025757'='Severe',
        'GSM2025758'='Recover',
        'GSM2025759'='Recover',
        'GSM2025760'='Recover',
        'GSM2025761'='Recover',
        'GSM2025762'='Recover',
        'GSM2025763'='Recover',
        'GSM2025764'='Recover',
        'GSM2025765'='Recover',
        'GSM2025766'='Recover',
        'GSM2025767'='Recover',
        'GSM2025768'='Recover',
        'GSM2025769'='Recover',
        'GSM2025770'='Healthy',
        'GSM2025771'='Healthy',
        'GSM2025772'='Healthy',
        'GSM2025773'='Healthy',
        'GSM2025774'='Healthy',
        'GSM2025775'='Healthy',
        'GSM2025776'='Healthy',
        'GSM2025777'='Healthy',
        'GSM2025778'='Healthy',
        'GSM2025779'='Healthy',
        'GSM2025780'='Healthy'
       
       )
sample_all[['Status']] <- unname(meta_annation[sample_all@meta.data$orig.ident])

saveRDS(sample_all,'/data_alluser/CXY/keti/singlecell_eqtl/EQTL/result/after_QC.rds')

sample_all<-readRDS('path/after_QC.rds')

sample_all@meta.data

sample_all_sct <- SCTransform(sample_all, vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>% 
                  RunPCA(npcs = 50, verbose = F) %>%
                  RunHarmony(group.by.vars = "orig.ident",assay.use = "SCT",max.iter.harmony = 20)

sample_all_sct

plot_size(10,5)
ElbowPlot(sample_all_sct,ndims=50)


dim = 25
sample_all_sct <- FindNeighbors(sample_all_sct,dims=1:dim,assay='SCT',reduction = 'harmony') %>%
    FindClusters(resolution = 0.1) %>%
    RunTSNE(dims = 1:dim,reduction = 'harmony')%>%
    RunUMAP(dims = 1:dim,reduction = 'harmony')

plot_size(20,8)
tsne <- DimPlot(sample_all_sct, reduction = "tsne", group.by = "seurat_clusters",
    pt.size = 1, label = T ,label.size = 4,label.box = F)
umap <- DimPlot(sample_all_sct, reduction = "umap", group.by = "seurat_clusters",
    pt.size = 1, label = T ,label.size = 4,label.box = F)
p1 <- tsne | umap
p1

sample_all_sct

saveRDS(sample_all_sct,'path/after_reduction.rds')

sample_all_sct<-readRDS('path/after_reduction.rds')



clusters_markers <- cosg(sample_all_sct,groups="all", slot = "counts", assay = "RNA",mu=1,n_genes_user=10)

clusters_markers 

write.csv(clusters_markers,
         '/data_alluser/CXY/keti/singlecell_eqtl/result/marker/clusters_markers_top10.csv',
         row.names=F)
#write.csv(clusters_markers,
        # '/data_alluser/CXY/cxy/DataSet/linbajie/csv/AllCells/clusters_markers.csv',
        # row.names=F)

 #lymphoid
"CD38","SLAMF7","POU2AF1",               #Plasma cell
"CD19","MS4A1","CD79A","CD79B","BLK",               #B cell 
"LILRA4","PLD4","LAMP5","IRF4","SMPD3","SPIB","PACSIN1","ITM2C","PTCRA","TCL1A","CD123","IL3RA",#pDC  "DERL3","GZMB",
"NCAM1","NKG7","ITGAL","ITGB2","KLRD1","KLRF1","CD247","GNLY","PRF1",        #NK cell
"AC092580.4","RPLP1","IL32","CCL5","RPS12","RPL41","RPS27","TMSB4X","TPT1","RPS29","KLRB1","CORO1A","ZNF683",#NKT  
"CD3D","CD3E","CD3G","TRAC","CD4","IL7R",                #CD4+ T cell
"CD8A","CD8B","GZMA","FYN",             #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",
"CD2","CTLA4","ICOS","CD28","SIRPG","GIMAP5","LTB","IL2RA","PTPRCAP","SPOCK2",#Treg  "CD3G","CD3D",

#Myeloid
"CD1C","CD1E","CCL17","RNASE6","FCER1A","CLEC10A","HLA-DQB2","WFDC21P","CD1A","CSF2RA","CCL22","S100B", "HLA-DPB1","BIRC3","CLEC9A","THBD","CCL19",#cDC 
"VCAN","CFP","CD300E","MPEG1","FGL2","S100A12","LILRA5","LGALS2","THBS1","FCAR","CD16","FCN1","APOBEC3A","KLRC4-KLRK1",#Monocyte  "RNASE6",
"FUT4","MPO","CEACAM8","ELANE","CXCR1","CXCR2",  #Neutrophil
"CD68","CD163","C1QA","MRC1","MS4A6A","MSR1","MERTK","SLC11A1",             #Macrophage
"LILRA1","LILRB2","NAPSB","ADGRE2","PSAP","LRRC25","GPBAR1",       #CD16+ monocyte  "CFP","LILRA5",
"KIT","CPA3","AREG","SAMSN1",                  #Mast cell

"GNLY", "PRF1" ,"SPON2" ,"GZMB",  "IL2RB" , "MYOM2", "ARL4C" ,"CD247" ,"NCAM1", #NK
"MS4A1" ,"BANK1" ,"CD79A", "CD74", "PAX5" ,"RALGPS2" ,"CD22", "IGHD" ,"CD79B", "HLA-DRA" ,"FCRL1" ,#B
"CD83", "HLA-DQB1" , "HLA-DQA1", 'HLA-DRB1', #B  

"IGHA1" ,"IGJ" ,"IGLC3" ,"IGKC", "IGHM" ,"IGHG1" ,"IGHG2", "IGHA2", "IGLC2" ,"IGHG3" ,"IGHG4" ,#IgA Plasma Cells
"MZB1",  #IgM Plasma Cells
"HSP90B1" ,      #IgG Plasma Cells
    
"CD163",  "HLA-DPA1", "THBS1", "EREG" ,"EMP1", "IL1R2", "HLA-DRB5",#CDC
"IRF8", "TCF4", "UGCG" ,"SERPINF1", "CCDC50", "ITM2C" ,"LILRA4", "PTPRS" ,"SPIB","MAP1A", #pDCs
    
"CDK6", "GATA2", "SOX4" ,"LAPTM4B" ,"MYB" ,"ANKRD28", "MSRB3", "H1F0", "PTRF", "FHL1" ,#HSCs
"HBB" ,"HBA2" ,"HBA1", "HBD" ,"ALAS2" ,"SLC25A37" ,"SNCA", "DCAF12", "CA1" ,"AHSP", #RBCs
"LCN2", "S100P", "BASP1" ,"CSF3R" , "GCA"  ,"MSRB1","EEF2", "NCF4" ,#Neutrophils
"FAM129C" ,"FCRLA",  #Myeloid/B-cell progenitor
"MKI67", "STMN1" ,"TOP2A", "CENPF" ,"TUBB" ,"TUBA1B" ,"MCM4", "PCNA", "SMC4" ,"RRM2" ,#Proliferating Lymphocytes
"TUBB1", "PPBP" ,"PF4" ,"SDPR" ,"ITGA2B", "PRKAR2B" ,"F13A1", "NRGN", "HIST1H2AC", "RGS18" #Platelets


plot_size(35,15)
DotPlot(sample_all_sct,c(    
"S100A9" ,"S100A8", "LYZ" ,"FCN1", "VCAN" ,"IL8" ,#CD14
"IFI27" , "IFITM3" ,"IFI6" ,  "PLBD1" ,"MNDA" ,"TYMP" ,#Inflammatory CD14
"IFI30" ,"S100A12" , "C19orf59", "FPR1" ,#Inflammatory CD14
"FCGR3A" ,"MS4A7" ,"CSF1R" ,"CDKN1C", "LST1" ,"C1QA", "WARS", "LILRB2", "PSAP" ,"AIF1",  #CD16
    
"TCF7" ,"IL7R" ,"LEF1", "DGKA", "RPSAP58" ,"LTB" ,"CCR7", "RPS3A" ,"RPL18A" ,"ABLIM1", #CD4
"CD8A" ,"GZMH", "CCL5" ,"TRGC2", "SYNE1", "TRGC1" ,"CD3G" ,"SYNE2", "NKG7", "FGFBP2", #CD8
"TNFAIP3" ,"CXCR4", "GZMA" ,"GZMK" ,"ZFP36L2" ,"SUN2" ,"CD3D",
"RPL30"  ,"XIST", "ITGA6" ,"ETS1", "IL6ST" , "PIK3IP1" ,"LDHB" ,#Activated CD4
    
"GNLY", "PRF1" ,"SPON2" ,"GZMB",  "IL2RB" , "MYOM2", "ARL4C" ,"CD247" ,"NCAM1", #NK
"MS4A1" ,"BANK1" ,"CD79A", "CD74", "PAX5" ,"RALGPS2" ,"CD22", "IGHD" ,"CD79B", "HLA-DRA" ,"FCRL1" ,#B
"CD83", "HLA-DQB1" , "HLA-DQA1", 'HLA-DRB1', #B  

"IGHA1" ,"IGJ" ,"IGLC3" ,"IGKC", "IGHM" ,"IGHG1" ,"IGHG2", "IGHA2", "IGLC2" ,"IGHG3" ,"IGHG4" ,#IgA Plasma Cells
"MZB1",  #IgM Plasma Cells
"HSP90B1" ,      #IgG Plasma Cells
    
"CD163",  "HLA-DPA1", "THBS1", "EREG" ,"EMP1", "IL1R2", "HLA-DRB5",#CDC
"IRF8", "TCF4", "UGCG" ,"SERPINF1", "CCDC50", "ITM2C" ,"LILRA4", "PTPRS" ,"SPIB","MAP1A", #pDCs
    
"CDK6", "GATA2", "SOX4" ,"LAPTM4B" ,"MYB" ,"ANKRD28", "MSRB3", "H1F0", "PTRF", "FHL1" ,#HSCs
"HBB" ,"HBA2" ,"HBA1", "HBD" ,"ALAS2" ,"SLC25A37" ,"SNCA", "DCAF12", "CA1" ,"AHSP", #RBCs
"LCN2", "S100P", "BASP1" ,"CSF3R" , "GCA"  ,"MSRB1","EEF2", "NCF4" ,#Neutrophils
"FAM129C" ,"FCRLA",  #Myeloid/B-cell progenitor
#"MKI67", "STMN1" ,"TOP2A", "CENPF" ,"TUBB" ,"TUBA1B" ,"MCM4", "PCNA", "SMC4" ,"RRM2" ,#Proliferating Lymphocytes
"TUBB1", "PPBP" ,"PF4" ,"SDPR" ,"ITGA2B", "PRKAR2B" ,"F13A1", "NRGN", "HIST1H2AC", "RGS18" #Platelets
),group.by='seurat_clusters',assay='SCT')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(40,20)
DotPlot(sample_all_sct,
        c(  #lymphoid
"CD38","SLAMF7","POU2AF1",               #Plasma cell
"CD19","MS4A1","CD79A","CD79B","BLK",               #B cell 
"LILRA4","PLD4","LAMP5","IRF4","SMPD3","SPIB","PACSIN1","ITM2C","PTCRA","TCL1A","CD123","IL3RA",#pDC  "DERL3","GZMB",
"NCAM1","NKG7","ITGAL","ITGB2","KLRD1","KLRF1","CD247","GNLY","PRF1",        #NK cell
"AC092580.4","RPLP1","IL32","CCL5","RPS12","RPL41","RPS27","TMSB4X","TPT1","RPS29","KLRB1","CORO1A","ZNF683","LAG3",""#NKT  
"CD3D","CD3E","CD3G","TRAC","CD4","IL7R",                #CD4+ T cell
"CD8A","CD8B","GZMA","FYN",             #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",
"CD2","CTLA4","ICOS","CD28","SIRPG","GIMAP5","LTB","IL2RA","PTPRCAP","SPOCK2",#Treg  "CD3G","CD3D",

#Myeloid
"CD1C","CD1E","CCL17","RNASE6","FCER1A","CLEC10A","HLA-DQB2","WFDC21P","CD1A","CSF2RA","CCL22","S100B", "HLA-DPB1","BIRC3","CLEC9A","THBD","CCL19",#cDC 
"VCAN","CFP","CD300E","MPEG1","FGL2","S100A12","LILRA5","LGALS2","THBS1","FCAR","CD16","FCN1","APOBEC3A","KLRC4-KLRK1",#Monocyte  "RNASE6",
"FUT4","MPO","CEACAM8","ELANE","CXCR1","CXCR2","FCGR3B", #Neutrophil
"CD68","CD163","C1QA","MRC1","MS4A6A","MSR1","MERTK","SLC11A1",             #Macrophage
"LILRA1","LILRB2","NAPSB","ADGRE2","PSAP","LRRC25","GPBAR1",       #CD16+ monocyte  "CFP","LILRA5",
"KIT","CPA3","AREG","SAMSN1",                  #Mast cell

 'PPBP','GP9',"PF4","CST3",#Platelets
   'CD34',"CD133","ABCG2","STRO-1","NCAM"     #Stem Cells  
          
         ),group.by='seurat_clusters',assay='SCT')+
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=25),axis.text.y = element_text(size=30))

 plot_size(20,8)
DotPlot(sample_all_sct,
        c(  
            "NCAM1","CD160",#NK
            "LAG3", #NKT  "CD160","CD8A","CD8B",
            "CD68","CSF1R","IL1B","CD14","CD33","HLA-DMA", #classical monocytes
            #Macrophages                        
"CTSD", "KCNMA1", "CTSB", "MRC1", "APOE" ,"CTSL", #Macrophages
            "CLEC10A","S100A9", "IRF8", "TCF4", "UGCG" ,"SERPINF1", "CCDC50", "ITM2C" ,"LILRA4", "PTPRS" ,"SPIB","MAP1A", #pDCs
            
            "CD19","MS4A1","CD79A","CD79B","BLK", #B cell 
             'PPBP','GP9',"PF4","CST3",#Platelets
 "CD1C","CD1E","CCL17","RNASE6","FCER1A","HLA-DQB2","WFDC21P","CD1A","CSF2RA","CCL22","S100B", "HLA-DPB1","BIRC3","CLEC9A","THBD","CCL19",#cDC 
 "CD3D","CD3E","CD3G","TRAC","CD4","IL7R", "FOXP3","IL2RA",#Treg "CD4","IL7R",               #CD4+ T cell
"CD8A","CD8B","GZMA","FYN"            #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",         
          
          
         ),group.by='seurat_clusters',assay='SCT')+
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))


plot_size(25,50)
FeaturePlot(sample_all_sct,c(
 "NCAM1","CD160",#NK
            "LAG3", #NKT  "CD160","CD8A","CD8B",
            "CD68","CSF1R","IL1B","CD14","CD33","HLA-DMA", #classical monocytes
            "CLEC10A","S100A9",#pDC
            "CTSD", "KCNMA1", "CTSB", "MRC1", "APOE" ,"CTSL", #Macrophages
            "CD19","MS4A1","CD79A","CD79B","BLK", #B cell 
             'PPBP','GP9',"PF4","CST3",#Platelets
 "CD1C","CD1E","CCL17","RNASE6","FCER1A","HLA-DQB2","WFDC21P","CD1A","CSF2RA","CCL22","S100B", "HLA-DPB1","BIRC3","CLEC9A","THBD","CCL19",#cDC 
 "CD3D","CD3E","CD3G","TRAC","CD4","IL7R", "FOXP3","IL2RA",#Treg "CD4","IL7R",               #CD4+ T cell
"CD8A","CD8B","GZMA","FYN"             #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",

                         ),reduction='umap',pt.size=0.3,ncol=4)

plot_size(25,20)

FeaturePlot(sample_all_sct,c(
 "CD3D","CD3E","CD3G","TRAC","CD4","IL7R",                #CD4+ T cell
"CD8A","CD8B","GZMA","FYN",             #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",
"CD2","CTLA4","ICOS","CD28","SIRPG","GIMAP5","LTB","IL2RA","PTPRCAP","SPOCK2"#Treg  "CD3G","CD3D"
   
    

                         ),reduction='umap',pt.size=0.3,ncol=4)



plot_size(25,100)

FeaturePlot(sample_all_sct,c(
 #lymphoid
"CD38","SLAMF7","POU2AF1",               #Plasma cell
"CD19","MS4A1","CD79A","CD79B","BLK",               #B cell 
"LILRA4","PLD4","LAMP5","IRF4","SMPD3","SPIB","PACSIN1","ITM2C","PTCRA","TCL1A","CD123","IL3RA",#pDC  "DERL3","GZMB",
"NCAM1","NKG7","ITGAL","ITGB2","KLRD1","KLRF1","CD247","GNLY","PRF1",        #NK cell
"AC092580.4","RPLP1","IL32","CCL5","RPS12","RPL41","RPS27","TMSB4X","TPT1","RPS29","KLRB1","CORO1A","ZNF683",#NKT  
"CD3D","CD3E","CD3G","TRAC","CD4","IL7R",                #CD4+ T cell
"CD8A","CD8B","GZMA","FYN",             #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",
"CD2","CTLA4","ICOS","CD28","SIRPG","GIMAP5","LTB","IL2RA","PTPRCAP","SPOCK2",#Treg  "CD3G","CD3D",

#Myeloid
"CD1C","CD1E","CCL17","RNASE6","FCER1A","CLEC10A","HLA-DQB2","WFDC21P","CD1A","CSF2RA","CCL22","S100B", "HLA-DPB1","BIRC3","CLEC9A","THBD","CCL19",#cDC 
"VCAN","CFP","CD300E","MPEG1","FGL2","S100A12","LILRA5","LGALS2","THBS1","FCAR","CD16","FCN1","APOBEC3A","KLRC4-KLRK1",#Monocyte  "RNASE6",
"FUT4","MPO","CEACAM8","ELANE","CXCR1","CXCR2","FCGR3B", #Neutrophil
"CD68","CD163","C1QA","MRC1","MS4A6A","MSR1","MERTK","SLC11A1",             #Macrophage
"LILRA1","LILRB2","NAPSB","ADGRE2","PSAP","LRRC25","GPBAR1",       #CD16+ monocyte  "CFP","LILRA5",
"KIT","CPA3","AREG","SAMSN1",                  #Mast cell

 'PPBP','GP9',"PF4","CST3",#Platelets
   'CD34',"CD133","ABCG2","STRO-1","NCAM"     #Stem Cells  
   
    

                         ),reduction='umap',pt.size=0.3,ncol=4)


cluster_annatation=c(
                                    '0'='CD4 T',
                                    '1'='Monocyte',
                                    '2'='CD8 T',
                                    '3'='B',
                                    '4'='Monocyte',
                                    '5'='CD8 T',
                                    '6'='CD8 T',
                                    '7'='NK',
                                    '8'='Platelet',
                                    '9'='CD4 T',
                                    '10'='cDC',
                                    '11'='CD8 T',
                                    '12'='B',
                                    '13'='pDC',
                                    '14'='B',
                                    '15'='B',
                                    '16'='CD8 T',
                                    '17'='Monocyte',
                                    '18'='CD4 T'
)

sample_all_sct[['Type']] <- unname(cluster_annatation[sample_all_sct@meta.data$seurat_clusters])



sample_all_sct<-readRDS('path/after_reduction.rds')


saveRDS(sample_all_sct,'path/after_annotation_Type.rds',compress=F)



unique(sample_all_sct@meta.data)

plot_size(20,8)

tsne <- DimPlot(sample_all_sct, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label = T ,label.size = 4,label.box = F)

umap <- DimPlot(sample_all_sct, reduction = "umap", group.by = "Type",

    pt.size = 0.4, label = T ,label.size = 4,label.box = F)




p1 <- tsne | umap 
p1

ggsave("/path/UAMP_Type.pdf", plot = Type, width = 10, height = 6)
ggsave("/path/UAMP_Type.png", plot = Type, width = 10, height = 6)

sample_all_sct<-readRDS("path/after_annotation_Type.rds")

sample_all_sct@meta.data

table(sample_all_sct@meta.data$Type=="CD8 T")



Idents(sample_all_sct) <- sample_all_sct$Type

B<- subset(sample_all_sct,idents='B') %>%
       SCTransform(vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>%
       RunPCA(npcs = 50,verbose = F)



plot_size(8,5)
ElbowPlot(B,ndims=50)

dim=25
B <- FindNeighbors(B, dims = 1:dim,assay = "SCT",reduction = 'harmony') %>%
    FindClusters(resolution = 0.2,assay = "SCT") %>%
    RunTSNE(reduction = "harmony",dims = 1:dim) %>%
    RunUMAP(reduction = "harmony",dims = 1:dim) 

plot_size(12,5)
tsne <- DimPlot(B,reduction = "tsne", group.by = "seurat_clusters",pt.size=1,label=T,label.size = 3,label.box = T)
umap <- DimPlot(B,reduction = "umap",group.by  = "seurat_clusters",pt.size=1,label=T,label.size = 3,label.box = T)
tsne|umap


plot_size(20,6)
DotPlot(B,
        feature=c( 
                  "TCL1A","FCER2","IL4R",#B_IN
            "TACI","CD27",#B_MEM
            "IGJ","BCMA"  #PLASMA
                 ),assay='SCT')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(20,6)
DotPlot(B,
        feature=c( 
                  'JCHAIN', "IGHG1","MZB1","SDC1",#Plasma  "CD79A"
            
                
            
  
          "CD19", "CD79A", "MS4A1", "FGR", "FCRL5","CXCR3", "CXCR4", "CXCR5", "CXCR6",#Memory "CD20",
                 'IGHD',"CD20","CD34","CD38","CD45R","TCL1A" ,"BTG1","SIGLEC6"#Naive B
                 ),assay='SCT')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(20,6)
DotPlot(B,
        feature=c(#"CD37", 'CD79A','CD79B',#B
                           'TNFRSF17','CD3E', 
                  'JCHAIN', "IGHG1","MZB1","SDC1",#Plasma  "CD79A"
                  "IgA", "IgG", "IgE",  "CD27", "CD40", "CD80", "PDL-2","CXCR3", "CXCR4", "CXCR5", "CXCR6",'MS4A1',#Memory "CD20",
                  "CD1", "CD21", "Notch2" , #Marginal zone B cells  "CD27",
                  "IgD",  "CD5",  "CD24", "TLR4","IL-10",  #Regulatory B cells  "CD21","CD1",
                 "CD22", "CD23",  #Follicular B cells    "CD21", "IgD", 
        # "LILRA4","PLD4","LAMP5","IRF4","SMPD3","SPIB","PACSIN1","ITM2C","PTCRA","TCL1A","CD123","IL3RA","DERL3","GZMB",#pDC  
                 'IGHD',"IgCD19","CD20","CD34","CD38","CD45R" #Naive B
                 ),assay='SCT')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(20,9)
VlnPlot(B,c('CD19','CD79A','CD79B','MS4A1','IGHD', #B
                            'JCHAIN','TNFRSF17','CD38','CD3E', 
                  "IGHG1","MZB1","SDC1",#Plasma  "CD79A"
                  "IgA", "IgG", "IgE", "CD20", "CD27", "CD40", "CD80", "PDL-2","CXCR3", "CXCR4", "CXCR5", "CXCR6",#Memory
                  "CD1", "CD21", "Notch2" , #Marginal zone B cells  "CD27",
                  "IgD",  "CD5",  "CD24", "TLR4","IL-10",  #Regulatory B cells  "CD21","CD1",
                 "CD22", "CD23",  #Follicular B cells    "CD21", "IgD", 
         "LILRA4","PLD4","LAMP5","IRF4","SMPD3","SPIB","PACSIN1","ITM2C","PTCRA","TCL1A","CD123","IL3RA","DERL3","GZMB"#pDC  
                 
             )
                         ,pt.size = 0,ncol=6)+
    NoLegend()+
scale_x_discrete("")

plot_size(20,35)
FeaturePlot(B,feature=c(
                        "CD37", 'CD79A','CD79B',#B
                           'TNFRSF17','CD3E', 
                  'JCHAIN', "IGHG1","MZB1","SDC1",#Plasma  "CD79A"
                  "IgA", "IgG", "IgE",  "CD27", "CD40", "CD80", "PDL-2","CXCR3", "CXCR4", "CXCR5", "CXCR6",'MS4A1',#Memory "CD20",
                  "CD1", "CD21", "Notch2" , #Marginal zone B cells  "CD27",
                  "IgD",  "CD5",  "CD24", "TLR4","IL-10",  #Regulatory B cells  "CD21","CD1",
                 "CD22", "CD23",  #Follicular B cells    "CD21", "IgD", 
        # "LILRA4","PLD4","LAMP5","IRF4","SMPD3","SPIB","PACSIN1","ITM2C","PTCRA","TCL1A","CD123","IL3RA","DERL3","GZMB",#pDC  
                 'IGHD',"CD19","CD20","CD34","CD38","CD45R" #Naive B
                        
                 ),reduction = 'tsne',ncol=4)

clusters_markers_B <- FindAllMarkers(B,only.pos = T,test.use = 'wilcox',assay='RNA',slot='counts',
               logfc.threshold = 0.5,min.pct = 0.25)


clusters_markers_B_top10 <- clusters_markers_B %>% group_by(cluster) %>% top_n(n=50,wt=avg_log2FC)
subset(clusters_markers_B_top10,cluster=='9')

B_annotation=c(
    
    '0'='B naive',
    '1'='B naive',
    '2'='B naive',
    '3'='B naive',
    '4'='B naive',
    '5'='B naive',
    '6'='B naive',
    '7'='B naive',
    '8'='B naive',
    '9'='B memory',
    '10'='Plasma',
    '11'='B naive',
    '12'='B naive',
    '13'='B naive',
    '14'='B naive',
    '15'='B naive'
   

)
B[['subtype']] <- unname(B_annotation[B@meta.data$seurat_clusters])


saveRDS(B,'path/B_annotated.rds',compress = F)



    

CD4<- subset(sample_all_sct,idents='CD4 T') %>%
       SCTransform(vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>%
       RunPCA(npcs = 50,verbose = F)

CD41<- subset(sample_all_sct,idents='CD4 T')
CD41

plot_size(8,5)
ElbowPlot(CD4,ndims=50)

dim=25
CD4 <- FindNeighbors(CD4, dims = 1:dim,assay = "SCT",reduction = 'harmony') %>%
    FindClusters(resolution = 0.1,assay = "SCT") %>%
    RunTSNE(reduction = "harmony",dims = 1:dim) %>%
    RunUMAP(reduction = "harmony",dims = 1:dim)

plot_size(12,5)
tsne <- DimPlot(CD4,reduction = "tsne", group.by = "seurat_clusters",pt.size=1,label.size = 5,label.box = T)
umap <- DimPlot(CD4,reduction = "umap",group.by  = "seurat_clusters",pt.size=1,label.size = 5,label.box = T)
tsne|umap

plot_size(20,6)
DotPlot(CD4,
        feature=c( 

"SOX4", "ID2",
"KLRB1","GZMK","TNFSF13B",#CD4+ effector T cell  "CCR7","SELL",
"CCR7","SELL","LRRN3","GPR183","IL7R" #CD4+ naive T cell

              
                       
                 ),assay='RNA')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))




plot_size(20,6)
DotPlot(CD4,
        feature=c( 

"CD3E", "CD4", #CD4+ T cell：
"CCR7","SELL","CD5",#CD4+ naive T cell：
"CD44","IFNG","S100A4","GPR183","IL7R", #CD4+ memory T cell：
"FASLG", "FAS", #CD4+ effector T cell："IFNG",
  'TCF7','CD69',  #Tex/耗竭T细胞
"CXCR5","ICOS","IL21","BCL6" ,#CD4+ Tfh：
   "CCR10","CD25","CD52","CMTM7","STAT5","IL10","IL12"  #Regulatory T   
              
                       
                 ),assay='SCT')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))


plot_size(20,9)
VlnPlot(CD4,c("CD3E", "CD4", #CD4+ T cell：
"CCR7","SELL","CD5",#CD4+ naive T cell：
"CD44","IFNG","S100A4","GPR183", "IL7R",#CD4+ memory T cell：
"FASLG", "FAS", #CD4+ effector T cell："IFNG",

"FOXP3","IL2RA"# CD4+ regulatory T (Treg) cell ：
             )
                         ,pt.size = 0,ncol=6)+
    NoLegend()+
scale_x_discrete("")

plot_size(20,25)
FeaturePlot(CD4,feature=c( "CD3E", "CD4", #CD4+ T cell：
"CCR7","SELL","CD5",#CD4+ naive T cell：
"CD44","IFNG","S100A4","GPR183","IL7R", #CD4+ memory T cell：
"FASLG", "FAS", #CD4+ effector T cell："IFNG",

"FOXP3","IL2RA"# CD4+ regulatory T (Treg) cell ：
                 ),reduction = 'umap',ncol=4)

plot_size(20,25)
FeaturePlot(CD4,feature=c( 'CD3D','CD3E','CD3G',
                        'CD40LG','PDGFRA',#Th
                 'CD8A','CD8B',#Tc  "CD3D","CD3E",
                  "TRDV2","TRGV9" , #gdT 
                     "S100A4",  #cd4 naive 
                        "IL32","ISG15", #cd4 ifn激活
                  "CD27","CCR7",#Naive T   "CD8A","CD8B"
                  "CCR10","CD25","IL2RA","CD52","CMTM7","FOXP3","STAT5","IL10","IL12", #Regulatory T
                  'CD4',"IL7R",  #CD4+ memory "CD27","CCR7",
                 "FASLG", "IFNG","CD44hi","FAS",   #CD4+ effector T
                  "ZNF683",'GNLY', #NKT cells  "CD8A","CD8B",
                  "GZMK", #CD8+ T cells "CD8A","CD8B",
                  "SELL" ,    #CD8 naive    "CCR7",
                 "FGFBP2","FCGR3A","CX3CR1","NKG7","TYOBP","PRF1","CD160","CD247","CCL3","GZMB",#NK  "GNLY",
                  'TCF7','CD69',  #Tex/耗竭T细胞
                        'IFITM3','IFI27'
                 ),reduction = 'tsne',ncol=4)

#clusters_markers_Mes <- FindAllMarkers(Mes,only.pos = T,test.use = 'wilcox',assay='RNA',slot='counts',
               logfc.threshold = 0.5,min.pct = 0.25)
#clusters_markers_Mes_top10 <- clusters_markers_Mes %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
subset(clusters_markers_Mes_top10,cluster=='3')

CD4_annotation=c(
    '0'='CD4 T naive',
    '1'='CD4 T memory',
    '2'='CD4 T memory',
    '3'='CD4 T memory',
    '4'='CD4 T naive'
    
)
CD4[['subtype']] <- unname(CD4_annotation[CD4@meta.data$seurat_clusters])


plot_size(20,8)

tsne <- DimPlot(CD4, reduction = "tsne", group.by = "subtype",

    pt.size = 0.4, label = T ,label.size = 4,label.box = F)

umap <- DimPlot(CD4, reduction = "umap", group.by = "subtype",

    pt.size = 0.4, label = T ,label.size = 4,label.box = F)

p1 <- tsne | umap

p1

saveRDS(CD4,'path/CD4_annotated.rds')

CD4<-readRDS("path/CD4_annotated.rds")

CD4







Monocyte<-readRDS("path/Monocyte_annotated.rds")

Monocyte<- subset(sample_all_sct,idents='Monocyte') %>%
       SCTransform(vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>%
       RunPCA(npcs = 50,verbose = F)

plot_size(8,5)
ElbowPlot(Monocyte,ndims=50)

dim=25
Monocyte <- FindNeighbors(Monocyte, dims = 1:dim,assay = "SCT",reduction = 'harmony') %>%
    FindClusters(resolution = 0.1,assay = "SCT") %>%
    RunTSNE(reduction = "harmony",dims = 1:dim) %>%
    RunUMAP(reduction = "harmony",dims = 1:dim) 

plot_size(12,5)
tsne <- DimPlot(Monocyte,reduction = "tsne", group.by = "seurat_clusters",pt.size=1,label.size = 5,label.box = T)
umap <- DimPlot(Monocyte,reduction = "umap",group.by  = "seurat_clusters",pt.size=1,label.size = 5,label.box = T)
tsne|umap

plot_size(20,6)
DotPlot(Monocyte,
        feature=c("CD14","CD16","LYZ","CCR2","CX3CR1"
                 ),assay='RNA')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(20,6)
DotPlot(Monocyte,
        feature=c("S100A9" ,"S100A8", "LYZ" ,"FCN1", "VCAN" ,"IL8" ,#CD14
"IFI27" , "IFITM3" ,"IFI6" ,  "PLBD1" ,"MNDA" ,"TYMP" ,#Inflammatory CD14
"IFI30" ,"S100A12" , "C19orf59", "FPR1" ,#Inflammatory CD14
"FCGR3A" ,"MS4A7" ,"CSF1R" ,"CDKN1C", "LST1" ,"C1QA", "WARS", "LILRB2", "PSAP" ,"AIF1"  ,"FCGR3"#CD16
                 ),assay='RNA')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(20,25)
FeaturePlot(Monocyte,feature=c( "S100A9" ,"S100A8", "LYZ" ,"FCN1", "VCAN" ,"IL8" ,#CD14
"IFI27" , "IFITM3" ,"IFI6" ,  "PLBD1" ,"MNDA" ,"TYMP" ,#Inflammatory CD14
"IFI30" ,"S100A12" , "C19orf59", "FPR1" ,#Inflammatory CD14
"FCGR3A" ,"MS4A7" ,"CSF1R" ,"CDKN1C", "LST1" ,"C1QA", "WARS", "LILRB2", "PSAP" ,"AIF1"  #CD16
                 ),reduction = 'tsne',ncol=4)

clusters_markers_Mes <- FindAllMarkers(Mes,only.pos = T,test.use = 'wilcox',assay='RNA',slot='counts',
               logfc.threshold = 0.5,min.pct = 0.25)
clusters_markers_Mes_top10 <- clusters_markers_Mes %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
subset(clusters_markers_Mes_top10,cluster=='3')

Monocyte_annotation=c(
    '0'='Monocyte classical',
    '1'='Monocyte classical',
    '2'='Monocyte nonclassical',
    '3'='Monocyte classical',
    '4'='Monocyte classical',
    '5'='Monocyte classical',
    '6'='Monocyte classical',
    '7'='Monocyte classical',
    '8'='Monocyte classical'
   
)
Monocyte[['subtype']] <- unname(Monocyte_annotation[Monocyte@meta.data$seurat_clusters])

saveRDS(Monocyte,'path/Monocyte_annotated.rds',compress = F)

CD8<- subset(sample_all_sct,idents='CD8 T') %>%
       SCTransform(vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>%
       RunPCA(npcs = 50,verbose = F)

plot_size(8,5)
ElbowPlot(CD8,ndims=50)

dim=25
CD8 <- FindNeighbors(CD8, dims = 1:dim,assay = "SCT",reduction = 'harmony') %>%
    FindClusters(resolution = 0.1,assay = "SCT") %>%
    RunTSNE(reduction = "harmony",dims = 1:dim) %>%
    RunUMAP(reduction = "harmony",dims = 1:dim) 

plot_size(12,5)
tsne <- DimPlot(CD8,reduction = "tsne", group.by = "seurat_clusters",pt.size=1,label.size = 5,label.box = T)
umap <- DimPlot(CD8,reduction = "umap",group.by  = "seurat_clusters",pt.size=1,label.size = 5,label.box = T)
tsne|umap

plot_size(20,6)
DotPlot(CD8,
        feature=c("CD3E","CD8A", #CD8+ T cell：
"CCR7", "SELL", "LTB","PASK", #CD8+ naïve：
"FASLG","CD44","GNLY","NKG7",# CD8+ effector ：
                  'TCF7','CD69'  #Tex/耗竭T细胞
                 ),assay='RNA')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))


plot_size(20,6)
DotPlot(CD8,
        feature=c("CD3E","CD8A", #CD8+ T cell：
"CCR7", "SELL", "LTB","PASK", #CD8+ naïve：
"FASLG","CD44","GNLY","NKG7",# FAS CD8+ effector ：
                  'TCF7','CD69'  #Tex/耗竭T细胞
                 ),assay='RNA')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(20,6)
DotPlot(CD8,
        feature=c( 'CD3D','CD3E','CD3G',
                        'CD40LG','PDGFRA',#Th
                 'CD8A','CD8B',#Tc  "CD3D","CD3E",
                  "TRDV2","TRGV9" , #gdT 
                  "CD27","CCR7",#Naive T   "CD8A","CD8B"
                  "CCR10","CD25","IL2RA","CD52","CMTM7","FOXP3","STAT5","IL10","IL12","CTLA4", #Regulatory T
                  'CD4',"IL7R",  #CD4+ memory "CD27","CCR7",
                 "FASLG", "IFNG","CD44hi","FAS",   #CD4+ effector T
                  "ZNF683",'GNLY', #NKT cells  "CD8A","CD8B",
                  "GZMK","CXCR3" ,#CD8+ T cells "CD8A","CD8B",
                  "SELL" ,    #CD8 naive    "CCR7",
                # "FGFBP2","FCGR3A","CX3CR1","NKG7","TYOBP","PRF1","CD160","CD247","CCL3","GZMB",#NK  "GNLY",
                  "MK167","TYMS",#Proliferating T 
                         'TCF7','CD69','IFITM3','ISG15','IFI27'
                 ),assay='RNA')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(20,15)
FeaturePlot(CD8,feature=c( "CD3E","CD8A", #CD8+ T cell：
"CCR7", "SELL", "LTB","PASK", #CD8+ naive：
"FASLG","CD44","GNLY","NKG7",# FAS CD8+ effector ：

                  'TCF7','CD69'  #Tex/耗竭T细胞
                 ),reduction = 'tsne',ncol=4,label=T)

clusters_markers_Mes <- FindAllMarkers(Mes,only.pos = T,test.use = 'wilcox',assay='RNA',slot='counts',
               logfc.threshold = 0.5,min.pct = 0.25)
clusters_markers_Mes_top10 <- clusters_markers_Mes %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
subset(clusters_markers_Mes_top10,cluster=='3')

CD8_annotation=c(
             '0'='CD8 T effector',
             '1'='CD8 T effector',
             '2'='CD8 T naive',
             '3'='CD8 T naive',
             '4'='CD8 T naive',
             '5'='CD8 T effector'
)
CD8[['subtype']] <- unname(CD8_annotation[CD8@meta.data$seurat_clusters])

saveRDS(CD8,'path/CD8_annotated.rds',compress = F)

CD8<-readRDS("/data_alluser/CXY/keti/singlecell_eqtl/resBult/RDS/CD8_annotated.rds")

CD8

NK<- subset(sample_all_sct,idents='NK') %>%
       SCTransform(vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>%
       RunPCA(npcs = 50,verbose = F)

plot_size(8,5)
ElbowPlot(NK,ndims=50)

dim=25
NK <- FindNeighbors(NK, dims = 1:dim,assay = "SCT",reduction = 'harmony') %>%
    FindClusters(resolution = 0.1,assay = "SCT") %>%
    RunTSNE(reduction = "harmony",dims = 1:dim) %>%
    RunUMAP(reduction = "harmony",dims = 1:dim) 

plot_size(12,5)
tsne <- DimPlot(NK,reduction = "tsne", group.by = "seurat_clusters",pt.size=1,label.size = 5,label.box = T)
umap <- DimPlot(NK,reduction = "umap",group.by  = "seurat_clusters",pt.size=1,label.size = 5,label.box = T)
tsne|umap

plot_size(20,6)
DotPlot(NK,
        feature=c("GZMB","GZMa","XCL1","XCL2","GZMK"
                 ),assay='RNA')+
theme(axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(20,15)
FeaturePlot(NK,feature=c( "CD3E","CD8A", #CD8+ T cell：
"CCR7", "SELL", "LTB","PASK", #CD8+ naïve：
"FASLG","CD44","GNLY","NKG7",# FAS CD8+ effector ：
"IL7R",
                  'TCF7','CD69'  #Tex/耗竭T细胞
                 ),reduction = 'tsne',ncol=4)

clusters_markers_Mes <- FindAllMarkers(Mes,only.pos = T,test.use = 'wilcox',assay='RNA',slot='counts',
               logfc.threshold = 0.5,min.pct = 0.25)
clusters_markers_Mes_top10 <- clusters_markers_Mes %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
subset(clusters_markers_Mes_top10,cluster=='3')

CD8_annotation=c(
    '0'='NK',
    '1'='NK',
    '2'='NK',
    '3'='NK',
    '4'='NK',
    '5'='NK'
)
CD8[['subtype']] <- unname(CD8_annotation[CD8@meta.data$seurat_clusters])

saveRDS(CD8,'/data_alluser/CXY/keti/singlecell_eqtl/EQTL/result/CD8_annotated.rds')

NK <- subset(sample_all_sct,idents='NK')

NK[['subtype']] <- 'NK'

saveRDS(NK,'path/NK_annotated.rds',compress = F)



Platelet <- subset(sample_all_sct,idents='Platelet')

Platelet[['subtype']] <- 'Platelet'

saveRDS(Platelet,'path/Platelet_annotated.rds',compress = F)

cDC <- subset(sample_all_sct,idents='cDC')

cDC[['subtype']] <- 'cDC'

saveRDS(cDC,'path/cDC_annotated.rds',compress = F)

pDC <- subset(sample_all_sct,idents='pDC')

pDC[['subtype']] <- 'pDC'

saveRDS(pDC,'path/pDC_annotated.rds',compress = F)

CD4<-readRDS("path/CD4_annotated.rds")
CD4

CD8<-readRDS("path/CD8_annotated.rds")
CD8

Monocyte<-readRDS("path/Monocyte_annotated.rds")
Monocyte

Platelet<-readRDS("path/Platelet_annotated.rds")
Platelet

B<-readRDS("path/B_annotated.rds")
B

cDC<-readRDS("path/cDC_annotated.rds")
cDC

pDC<-readRDS("path/pDC_annotated.rds")
pDC

NK<-readRDS("path/NK_annotated.rds")
NK

sample_all_sct1<-readRDS("path/all_subtype_annotated.rds")

sample_all_sct1@meta.data

table(sample_all_sct1@meta.data$Type=="CD8 T")
table(sample_all_sct1@meta.data$subtype=="Platelet")
z<-sample_all_sct1@meta.data[,c(5,10,11)]
table(c(z[which(z$subtype=="CD4 T naive"),])$Status)





info <- rbind(CD4[[c('Type','subtype')]],Monocyte[[c('Type','subtype')]],
             B[[c('Type','subtype')]],CD8[[c('Type','subtype')]],
             Platelet[[c('Type','subtype')]],cDC[[c('Type','subtype')]],
             pDC[[c('Type','subtype')]], NK[[c('Type','subtype')]])

merged <- AddMetaData(sample_all_sct,metadata=info)

merged@meta.data

unique(merged@meta.data$subtype)

plot_size(20,8)

tsne <- DimPlot(merged, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label.size = 4,label.box = F)

umap <- DimPlot(merged, reduction = "umap", group.by = "subtype",

    pt.size = 0.4, label.size = 4,label.box = F)

p1 <- tsne | umap

p1

saveRDS(merged,'path/all_subtype_annotated.rds',compress = F)

#merged@meta.data$subtype[which(merged@meta.data$Type =='Platelet')] <- 'Platelet'

#merged@meta.data$subtype[which(merged@meta.data$Type =='cDC')] <- 'cDC'

#merged@meta.data$subtype[which(merged@meta.data$Type =='pDC')] <- 'pDC'

All_t_b<-merged

All_t_b@meta.data$Status1<-"NA"

All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025728")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025729")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025730")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025731")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025732")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025733")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025734")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025735")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025736")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025737")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025738")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025739")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025740")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025741")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025742")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025743")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025744")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025745")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025746")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025747")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025748")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025749")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025750")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025751")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025752")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025753")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025754")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025755")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025756")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025757")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025758")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025759")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025760")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025761")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025762")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025763")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025764")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025765")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025766")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025767")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025768")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025769")]<-"COVID-19"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025770")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025771")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025772")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025773")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025774")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025775")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025776")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025777")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025778")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025779")]<-"Healthy"
All_t_b@meta.data$Status1[which(All_t_b@meta.data$orig.ident=="GSM2025780")]<-"Healthy"

All_t_b@meta.data$orig.ident1<-"NA"

All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025728")]<-"asymptomatic case1"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025729")]<-"asymptomatic case2"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025730")]<-"asymptomatic case3"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025731")]<-"asymptomatic case4"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025732")]<-"asymptomatic case5"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025733")]<-"Moderate1"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025734")]<-"Moderate2"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025735")]<-"Moderate3"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025736")]<-"Moderate4"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025737")]<-"Moderate5"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025738")]<-"Moderate6"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025739")]<-"Moderate7"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025740")]<-"Moderate8"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025741")]<-"Moderate9"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025742")]<-"Moderate10"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025743")]<-"Moderate11"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025744")]<-"Moderate12"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025745")]<-"Moderate13"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025746")]<-"Severe1"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025747")]<-"Severe2"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025748")]<-"Severe3"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025749")]<-"Severe4"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025750")]<-"Severe5"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025751")]<-"Severe6"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025752")]<-"Severe7"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025753")]<-"Severe8"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025754")]<-"Severe9"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025755")]<-"Severe10"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025756")]<-"Severe11"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025757")]<-"Severe12"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025758")]<-"Recover1"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025759")]<-"Recover2"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025760")]<-"Recover3"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025761")]<-"Recover4"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025762")]<-"Recover5"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025763")]<-"Recover6"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025764")]<-"Recover7"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025765")]<-"Recover8"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025766")]<-"Recover9"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025767")]<-"Recover10"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025768")]<-"Recover11"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025769")]<-"Recover12"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025770")]<-"Healthy1"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025771")]<-"Healthy2"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025772")]<-"Healthy3"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025773")]<-"Healthy4"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025774")]<-"Healthy5"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025775")]<-"Healthy6"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025776")]<-"Healthy7"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025777")]<-"Healthy8"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025778")]<-"Healthy9"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025779")]<-"Healthy10"
All_t_b@meta.data$orig.ident1[which(All_t_b@meta.data$orig.ident=="GSM2025780")]<-"Healthy11"

All_t_b@meta.data


plot_size(20,8)

tsne <- DimPlot(All_t_b, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label.size = 4,label.box = F)

umap <- DimPlot(All_t_b, reduction = "umap", group.by = "Type",label = T,

    pt.size = 0.4, label.size = 4,label.box = F)

p1 <- tsne | umap

p1


saveRDS(All_t_b,'path/all_subtype_annotated.rds',compress = F)