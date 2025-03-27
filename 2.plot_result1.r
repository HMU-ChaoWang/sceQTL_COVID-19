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
library(paletteer)
library(scico)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(showtext)

plot_size <- function(x, y) {
    options(repr.plot.width = x, repr.plot.height = y)
}

path = "/path/result"
files <- list.files(path)
dir = paste0(path,'/',files)
names(dir) <- files


rownames(brewer.pal.info)
brewer.pal.info
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=FALSE)
#display.brewer.pal(8,"Set1")
brewer.pal(9,"Set1")

sce<-readRDS("/path/result/RDS/all_subtype_annotated.rds")

sce@meta.data

palettes_d_names[order(-palettes_d_names $length),][290:350,]

palettes_d_names[order(-palettes_d_names $length),][2210:2225,]


names <-palettes_d_names[2200:2225,]
for (i in 1:dim(names)[1]) {
  palette <-  paletteer_d(paste(names[i,1], names[i,2], sep = "::"))
scales::show_col(palette)
    
}


scales::show_col(colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(15))

paletteer_d("ggsci::default_nejm",n=8)
paletteer_d("basetheme::clean",n=8)
paletteer_d("RColorBrewer::Paired",n=8)
paletteer_c("scico::roma", n = 8)
paletteer_dynamic("cartography::green.pal", 8)
paletteer_c("scico::roma", n = 8)
paletteer_d("ggsci::default_nejm",n=8)
scico_palette_show()

colorRampPalette(c("blue", "red"))(100)
display.brewer.pal(12,"Paired")
scales::show_col(colorRampPalette(brewer.pal(12,"Paired"))(16))
scales::show_col(scico(13, begin = 0, end = 1, direction = -1, palette = "roma"))
x <- colorRampPalette(c("royalblue","firebrick3"))(16)
scales::show_col(x)



scales::show_col(brewer.pal(9,"Set1"))

scales::show_col(colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))
scales::show_col(brewer.pal(12,"Paired"))

color<-colorRampPalette(paletteer_d("ggsci::default_nejm",direction = -1))(13)
scales::show_col(color)

cloor(1,4,7,8,11,23,15,18,12,20,22,24,21)


#COLOR<-c("#F0C7DA","#C30E6F","#F3AD5E","#78AFD5","#EB736C","#BCBCBC","#249089","#82BF83","#B75C26","#EEEC96",
#        "#623C89","#C1B0CA","#EC7A1F","#5E3700","#E0091B","#EA8C91","#30983C","#A0CB8A","#0871B4" ,"#9CC8E0" )                


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ,"#9CC8E0")  

COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ) 
0871B4 82BF83 EC7A1F
scales::show_col(COLOR)

scales::show_col("#E41A1C")

plot_size(22,8)

color1<-c("#B6B6B7","#BDD7E7","#6BAED6","#3182BD","#08519C")

tsne_Status <- DimPlot(sce, reduction = "tsne", group.by = "Status",

   pt.size = 0.4, label.size = 4,label.box = F,cols=color1)

#tsne_Status <- DimPlot(sce, reduction = "tsne", group.by = "Status",

#    pt.size = 0.4, label.size = 4,label.box = F,cols=brewer.pal(12,"Paired"))



sample <- DimPlot(sce, reduction = "tsne", group.by = "orig.ident1",

    pt.size = 0.4, label.size = 4,label.box = F)



p <- sample| tsne_Status

p
ggsave("/path/result/PLOT/result_plot/blue_tsne_sample_status.pdf", plot = p, width = 15, height = 5)
ggsave("/path/result/PLOT/result_plot/blue_tsne_sample_status.png", plot = p, width = 15, height = 5)

color1<-c("#B6B6B7","#BDD7E7","#6BAED6","#3182BD","#08519C")
plot_size(22,8)

tsne_Status <- DimPlot(sce, reduction = "umap", group.by = "Status",

    pt.size = 0.4, label.size = 4,label.box = F,cols=color1)



sample <- DimPlot(sce, reduction = "umap", group.by = "orig.ident1",

    pt.size = 0.4, label.size = 4,label.box = F)



p <- sample| tsne_Status
p
ggsave("/path/result/PLOT/result_plot/blue_umap_sample_status.pdf", plot = p, width = 15, height = 5)
ggsave("/path/result/PLOT/result_plot/blue_umap_sample_status.png", plot = p, width = 15, height = 5)

plot_size(22,8)

tsne_Status <- DimPlot(sce, reduction = "umap", group.by = "Status",

    pt.size = 0.4, label.size = 4,label.box = F,cols=brewer.pal(12,"Paired"))



sample <- DimPlot(sce, reduction = "umap", group.by = "orig.ident1",

    pt.size = 0.4, label.size = 4,label.box = F)



p <- sample| tsne_Status
p
ggsave("/path/result/PLOT/result_plot/umap_sample_status.pdf", plot = p, width = 15, height = 5)
ggsave("/path/result/PLOT/result_plot/umap_sample_status.png", plot = p, width = 15, height = 5)

plot_size(22,8)

tsne_Status <- DimPlot(sce, reduction = "umap", group.by = "Status",

    pt.size = 0.4, label.size = 4,label.box = F,cols=brewer.pal(12,"Paired"))



sample <- DimPlot(sce, reduction = "umap", group.by = "orig.ident1",

    pt.size = 0.4, label.size = 4,label.box = F)



p <- sample| tsne_Status
p
ggsave("/path/result/PLOT/result_plot/umap_sample_status.pdf", plot = p, width = 15, height = 5)
ggsave("/path/result/PLOT/result_plot/umap_sample_status.png", plot = p, width = 15, height = 5)

plot_size(20,8)



umap_Type <- DimPlot(sce, reduction = "umap", group.by = "Type",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols= colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

p <-  umap_Type

p


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" )  


umap_subtype <- DimPlot(sce, reduction = "umap", group.by = "subtype",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

pp <-  umap_subtype

pp



ggsave("/path/result/PLOT/result_plot/umap_tsne_subtype.pdf", plot = pp, width = 15, height = 5)
ggsave("/path/result/PLOT/result_plot/umap_tsne_subtype.png", plot = pp, width = 15, height = 5)

plot_size(20,8)

tsne_Type <- DimPlot(sce, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label.size = 4,label.box = F,cols=colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

umap_Type <- DimPlot(sce, reduction = "umap", group.by = "Type",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols= colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

p <- tsne_Type | umap_Type

p


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" )  

tsne_subtype <- DimPlot(sce, reduction = "tsne", group.by = "subtype",

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

umap_subtype <- DimPlot(sce, reduction = "umap", group.by = "subtype",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

pp <- tsne_subtype | umap_subtype

pp
ggsave("/path/result/PLOT/result_plot/umap_tsne_type.pdf", plot = p, width = 15, height = 5)
ggsave("/path/result/PLOT/result_plot/umap_tsne_type.png", plot = p, width = 15, height = 5)


ggsave("/path/result/PLOT/result_plot/umap_tsne_subtype.pdf", plot = pp, width = 15, height = 5)
ggsave("/path/result/PLOT/result_plot/umap_tsne_subtype.png", plot = pp, width = 15, height = 5)

scales::show_col(colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

asymptomatic<- subset(sce,Status =="asymptomatic case" )

asymptomatic_sct <- SCTransform(asymptomatic, vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>% 
                  RunPCA(npcs = 50, verbose = F) %>%
                  RunHarmony(group.by.vars = "orig.ident",assay.use = "SCT",max.iter.harmony = 20)

dim = 25
asymptomatic_sct <- FindNeighbors(asymptomatic_sct,dims=1:dim,assay='SCT',reduction = 'harmony') %>%
    FindClusters(resolution = 0.1) %>%
    RunTSNE(dims = 1:dim,reduction = 'harmony')%>%
    RunUMAP(dims = 1:dim,reduction = 'harmony')

plot_size(20,8)
DotPlot(asymptomatic_sct,c(     

    'PPBP','GP9',"PF4",#Platelet  "CST3",
    "TCF4", "UGCG" ,"SERPINF1", "CCDC50",#pDCs
    "NCAM1","CD160","NKG7","CD247",#NK
    "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",#Monocytes   
    "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC  
    "CD3D","CD3E","CD3G", "CD8A","CD8B","GNLY", #CD8 T
    "IL7R","CD4","FOXP3",#CD4 T    
    "CD19","MS4A1","CD79A","CD79B","BLK" #B cell   


# "FOXP3","IL2RA","TRAC",#"CD4",             #CD4+ T cell
#"GZMA","FYN"            #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",      
   
    
    
                      
         ),group.by='Type',assay='SCT')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

saveRDS(asymptomatic_sct,'/path/result/RDS/asymptomatic_sct_reduction.rds',compress = F)



asymptomatic_sct<-readRDS('/path/result/RDS/asymptomatic_sct_reduction.rds')

asymptomatic_sct@meta.data





plot_size(20,8)

tsne_Type <- DimPlot(asymptomatic_sct, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label.size = 4,label.box = F,cols=colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

umap_Type <- DimPlot(asymptomatic_sct, reduction = "umap", group.by = "Type",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols= colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

p1 <- tsne_Type | umap_Type

p1


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" )  

tsne_subtype <- DimPlot(asymptomatic_sct, reduction = "tsne", group.by = "subtype",

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

umap_subtype <- DimPlot(asymptomatic_sct, reduction = "umap", group.by = "subtype",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

pp1 <- tsne_subtype | umap_subtype

pp1


Healthy<- subset(sce,Status =="Healthy" )

Healthy_sct <- SCTransform(Healthy, vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>% 
                  RunPCA(npcs = 50, verbose = F) %>%
                  RunHarmony(group.by.vars = "orig.ident",assay.use = "SCT",max.iter.harmony = 20)

dim = 25
Healthy_sct <- FindNeighbors(Healthy_sct,dims=1:dim,assay='SCT',reduction = 'harmony') %>%
    FindClusters(resolution = 0.1) %>%
    RunTSNE(dims = 1:dim,reduction = 'harmony')%>%
    RunUMAP(dims = 1:dim,reduction = 'harmony')

plot_size(20,8)
DotPlot(Healthy_sct,c(     

    'PPBP','GP9',"PF4",#Platelet  "CST3",
    "TCF4", "UGCG" ,"SERPINF1", "CCDC50",#pDCs
    "NCAM1","CD160","NKG7","CD247",#NK
    "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",#Monocytes   
    "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC  
    "CD3D","CD3E","CD3G", "CD8A","CD8B","GNLY", #CD8 T
    "IL7R","CD4","FOXP3",#CD4 T    
    "CD19","MS4A1","CD79A","CD79B","BLK" #B cell   


# "FOXP3","IL2RA","TRAC",#"CD4",             #CD4+ T cell
#"GZMA","FYN"            #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",      
   
    
    
                      
         ),group.by='Type',assay='SCT')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

saveRDS(Healthy_sct,'/path/result/RDS/Healthy_sct_reduction.rds',compress = F)

Healthy_sct<-readRDS('/path/result/RDS/Healthy_sct_reduction.rds')





plot_size(20,8)

tsne_Type <- DimPlot(Healthy_sct, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label.size = 4,label.box = F,cols=colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

umap_Type <- DimPlot(Healthy_sct, reduction = "umap", group.by = "Type",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols= colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

p1 <- tsne_Type | umap_Type

p1


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" )  

tsne_subtype <- DimPlot(Healthy_sct, reduction = "tsne", group.by = "subtype",

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

umap_subtype <- DimPlot(Healthy_sct, reduction = "umap", group.by = "subtype",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

pp1 <- tsne_subtype | umap_subtype

pp1



Moderate<- subset(sce,Status =="Moderate" )

Moderate_sct <- SCTransform(Moderate, vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>% 
                  RunPCA(npcs = 50, verbose = F) %>%
                  RunHarmony(group.by.vars = "orig.ident",assay.use = "SCT",max.iter.harmony = 20)

dim = 25
Moderate_sct <- FindNeighbors(Moderate_sct,dims=1:dim,assay='SCT',reduction = 'harmony') %>%
    FindClusters(resolution = 0.1) %>%
    RunTSNE(dims = 1:dim,reduction = 'harmony')%>%
    RunUMAP(dims = 1:dim,reduction = 'harmony')

plot_size(20,8)
DotPlot(Moderate_sct,c(     

    'PPBP','GP9',"PF4",#Platelet  "CST3",
    "TCF4", "UGCG" ,"SERPINF1", "CCDC50",#pDCs
    "NCAM1","CD160","NKG7","CD247",#NK
    "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",#Monocytes   
    "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC  
    "CD3D","CD3E","CD3G", "CD8A","CD8B","GNLY", #CD8 T
    "IL7R","CD4","FOXP3",#CD4 T    
    "CD19","MS4A1","CD79A","CD79B","BLK" #B cell   


# "FOXP3","IL2RA","TRAC",#"CD4",             #CD4+ T cell
#"GZMA","FYN"            #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",      
   
    
    
                      
         ),group.by='Type',assay='SCT')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

saveRDS(Moderate_sct,'/path/result/RDS/Moderate_sct_reduction.rds',compress = F)

Moderate_sct<-readRDS('/path/result/RDS/Moderate_sct_reduction.rds')




plot_size(20,8)

tsne_Type <- DimPlot(Moderate_sct, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label.size = 4,label.box = F,cols=colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

umap_Type <- DimPlot(Moderate_sct, reduction = "umap", group.by = "Type",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols= colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

p1 <- tsne_Type | umap_Type

p1


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" )  

tsne_subtype <- DimPlot(Moderate_sct, reduction = "tsne", group.by = "subtype",

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

umap_subtype <- DimPlot(Moderate_sct, reduction = "umap", group.by = "subtype",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

pp1 <- tsne_subtype | umap_subtype

pp1



Recover<- subset(sce,Status =="Recover" )

Recover_sct <- SCTransform(Recover, vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>% 
                  RunPCA(npcs = 50, verbose = F) %>%
                  RunHarmony(group.by.vars = "orig.ident",assay.use = "SCT",max.iter.harmony = 20)

dim = 25
Recover_sct <- FindNeighbors(Recover_sct,dims=1:dim,assay='SCT',reduction = 'harmony') %>%
    FindClusters(resolution = 0.1) %>%
    RunTSNE(dims = 1:dim,reduction = 'harmony')%>%
    RunUMAP(dims = 1:dim,reduction = 'harmony')

plot_size(20,8)
DotPlot(Recover_sct,c(     

    'PPBP','GP9',"PF4",#Platelet  "CST3",
    "TCF4", "UGCG" ,"SERPINF1", "CCDC50",#pDCs
    "NCAM1","CD160","NKG7","CD247",#NK
    "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",#Monocytes   
    "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC  
    "CD3D","CD3E","CD3G", "CD8A","CD8B","GNLY", #CD8 T
    "IL7R","CD4","FOXP3",#CD4 T    
    "CD19","MS4A1","CD79A","CD79B","BLK" #B cell   


# "FOXP3","IL2RA","TRAC",#"CD4",             #CD4+ T cell
#"GZMA","FYN"            #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",      
   
    
    
                      
         ),group.by='Type',assay='SCT')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

saveRDS(Recover_sct,'/path/result/RDS/Recover_sct_reduction.rds',compress = F)

Recover_sct<-readRDS('/path/result/RDS/Recover_sct_reduction.rds')





plot_size(20,8)

tsne_Type <- DimPlot(Recover_sct, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label.size = 4,label.box = F,cols=colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

umap_Type <- DimPlot(Recover_sct, reduction = "umap", group.by = "Type",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols= colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

p1 <- tsne_Type | umap_Type

p1


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" )  

tsne_subtype <- DimPlot(Recover_sct, reduction = "tsne", group.by = "subtype",

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

umap_subtype <- DimPlot(Recover_sct, reduction = "umap", group.by = "subtype",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

pp1 <- tsne_subtype | umap_subtype

pp1



Severe<- subset(sce,Status =="Severe" )

Severe_sct <- SCTransform(Severe, vars.to.regress = c("percent.MT",'orig.ident'), verbose = FALSE) %>% 
                  RunPCA(npcs = 50, verbose = F) %>%
                  RunHarmony(group.by.vars = "orig.ident",assay.use = "SCT",max.iter.harmony = 20)

dim = 25
Severe_sct <- FindNeighbors(Severe_sct,dims=1:dim,assay='SCT',reduction = 'harmony') %>%
    FindClusters(resolution = 0.1) %>%
    RunTSNE(dims = 1:dim,reduction = 'harmony')%>%
    RunUMAP(dims = 1:dim,reduction = 'harmony')

plot_size(20,8)
DotPlot(Severe_sct,c(     

    'PPBP','GP9',"PF4",#Platelet  "CST3",
    "TCF4", "UGCG" ,"SERPINF1", "CCDC50",#pDCs
    "NCAM1","CD160","NKG7","CD247",#NK
    "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",#Monocytes   
    "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC  
    "CD3D","CD3E","CD3G", "CD8A","CD8B","GNLY", #CD8 T
    "IL7R","CD4","FOXP3",#CD4 T    
    "CD19","MS4A1","CD79A","CD79B","BLK" #B cell   


# "FOXP3","IL2RA","TRAC",#"CD4",             #CD4+ T cell
#"GZMA","FYN"            #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",      
   
    
    
                      
         ),group.by='Type',assay='SCT')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

saveRDS(Severe_sct,'/path/result/RDS/Severe_sct_reduction.rds',compress = F)

"asymptomatic case","Healthy","Moderate","Recover","Severe"

Severe_sct<-readRDS('/path/result/RDS/Severe_sct_reduction.rds')





plot_size(20,8)

tsne_Type <- DimPlot(Severe_sct, reduction = "tsne", group.by = "Type",

    pt.size = 0.4, label.size = 4,label.box = F,cols=colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

umap_Type <- DimPlot(Severe_sct, reduction = "umap", group.by = "Type",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols= colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))

p1 <- tsne_Type | umap_Type

p1


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" )  

tsne_subtype <- DimPlot(Severe_sct, reduction = "tsne", group.by = "subtype",

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

umap_subtype <- DimPlot(Severe_sct, reduction = "umap", group.by = "subtype",label = T,

    pt.size = 0.4, label.size = 4,label.box = F,cols=COLOR)

pp1 <- tsne_subtype | umap_subtype

pp1

sce@meta.data


for(i in unique(sce$Type)){    
    obj <- subset(sce,idents=i)
    #obj@active.assay <- 'RNA'
    tmp <- try(marker <- FindMarkers(
        obj,ident.1 = "asymptomatic case",ident.2="Healthy",ident.2="Healthy",ident.2="Healthy",ident.2="Healthy",
                                     group.by='Status',slot='data',assay = 'SCT',logfc.threshold =0,min.pct = 0.25,test.use='wilcox'))

if('try-error' %in% class(tmp)){
  next
}
    print(i)
    write.csv(marker,paste0('/path/result/analy/new_type/',i,'.csv'))
    
        
}





Idents(sce) <- sce$Status1
unique(sce$Status1)
table(sce@meta.data$Status1)



#for(i in unique(sce$Status1)){    
   
    
    marker <- FindMarkers(
        sce,ident.1 = "COVID-19",ident.2="Healthy",
                                     group.by='Status1',slot='data',assay = 'SCT',logfc.threshold =0,min.pct = 0.25,test.use='wilcox')

 

    #print(i)
    write.csv(marker,paste0('/path/result/analy/new_Volcano/hc_covid/',"all_deg",'.csv'))
    
        
#}

path = '/path/result/analy/new_Volcano/hc_covid'
files <- list.files(path)
inter_val <- c()
for(i in files){

    name <- unlist(strsplit(i,split='.',fixed=T))[1]

    df <- read.csv(paste0(path,'/',i),header=T)

    up_gene_diff<-length(which(df$avg_log2FC>0))
    up_gene5.pct<-ceiling(up_gene_diff*0.2)
  
    down_gene_diff<-length(which(df$avg_log2FC<0))
    down_gene5.pct<-ceiling(down_gene_diff*0.2)
    
    gene_up <- df[which(df$avg_log2FC>0),]
    gene_up_5.pct <- gene_up[order(-gene_up$avg_log2FC),][1:up_gene5.pct,1]

    gene_down <- df[which(df$avg_log2FC<0),]
    gene_down_5.pct <- gene_down[order(gene_down$avg_log2FC),][1:down_gene5.pct,1]


    All_gene_5.pct <- df[which(df$X %in% c(gene_up_5.pct,gene_down_5.pct)),]


    FC_up <- All_gene_5.pct[which(All_gene_5.pct$avg_log2FC>0),]$avg_log2FC
    FC_down <- All_gene_5.pct[which(All_gene_5.pct$avg_log2FC<0),]$avg_log2FC

    cut_off_pvalue <- 0.05 
    cut_off_logFC_up <- 0.25   
    cut_off_logFC_down <- -0.25
    
    df$change <- ifelse(df$p_val < cut_off_pvalue &(df$avg_log2FC > cut_off_logFC_up | df$avg_log2FC < cut_off_logFC_down),
                       ifelse(df$avg_log2FC > cut_off_logFC_up,'Up','Down'),'Stable')

    df$xinterceptleft <- cut_off_logFC_down
    df$xinterceptright <- cut_off_logFC_up
    df$yintercept <- -log(cut_off_pvalue,10)
    df$type <- name
    df$`-log10(P.value) `<- -log(df$p_val,10)
    
   
    colnames(df)[1] <-'gene'
    inter_val <- rbind(inter_val,df)
    print(i)    
    write.csv(df,paste0('/path/result/analy/new_Volcano/hc_covid/','Volcano_',i),row.names=F,quote=F)
}
write.csv(inter_val,'/path/result/analy/new_Volcano/hc_covid/All_Volcano.csv',row.names=F,quote=F)







Idents(sce) <- sce$Type
unique(sce$Type)

table(sce@meta.data$Type)

table(sce@meta.data$Type,sce@meta.data$Status1)



<div class="girk">
obj[[]]</div><i class="fa fa-lightbulb-o "></i>



<div class="burk">

for(i in unique(sce$Type)){    
    obj <- subset(sce,idents=i)
    #obj@active.assay <- 'RNA'
    tmp <- try(marker <- FindMarkers(
        obj,ident.1 = "COVID-19",ident.2="Healthy",
                                     group.by='Status1',slot='data',assay = 'SCT',logfc.threshold =0,min.pct = 0.25,test.use='wilcox'))

 
if('try-error' %in% class(tmp)){
  next
}
    print(i)
    write.csv(marker,paste0('/path/result/analy/new_type/',i,'.csv'))
    
        
}

path = '/path/result/analy/new_type'
files <- list.files(path)
inter_val <- c()
for(i in files){
    name <- unlist(strsplit(i,split='.',fixed=T))[1]
    #读入文件
    df <- read.csv(paste0(path,'/',i),header=T)
        #上调基因降序排列，选取avg_log2FC前20%的基因。  
    up_gene_diff<-length(which(df$avg_log2FC>0))
    up_gene5.pct<-ceiling(up_gene_diff*0.2)
  
    down_gene_diff<-length(which(df$avg_log2FC<0))
    down_gene5.pct<-ceiling(down_gene_diff*0.2)
    
    gene_up <- df[which(df$avg_log2FC>0),]
    gene_up_5.pct <- gene_up[order(-gene_up$avg_log2FC),][1:up_gene5.pct,1]
    #下调基因升序排列，选取avg_log2FC前5%的基因 。 
    gene_down <- df[which(df$avg_log2FC<0),]
    gene_down_5.pct <- gene_down[order(gene_down$avg_log2FC),][1:down_gene5.pct,1]

    #把上下调的基因行提取出来
    All_gene_5.pct <- df[which(df$X %in% c(gene_up_5.pct,gene_down_5.pct)),]

    #上下调基因的avg_log2FC分别提取出来
    FC_up <- All_gene_5.pct[which(All_gene_5.pct$avg_log2FC>0),]$avg_log2FC
    FC_down <- All_gene_5.pct[which(All_gene_5.pct$avg_log2FC<0),]$avg_log2FC

    cut_off_pvalue <- 0.05 #统计显著性
    cut_off_logFC_up <- 0.25   #差异倍数值
    cut_off_logFC_down <- -0.25
    
    #根据阈值参数，上调的基因设为up，下调基因设为down，无差异设为stable，保存在change列中。
    df$change <- ifelse(df$p_val < cut_off_pvalue &(df$avg_log2FC > cut_off_logFC_up | df$avg_log2FC < cut_off_logFC_down),
                       ifelse(df$avg_log2FC > cut_off_logFC_up,'Up','Down'),'Stable')
    #补全xinterceptleft、xinterceptright、yintercept、还有细胞亚型名
    df$xinterceptleft <- cut_off_logFC_down
    df$xinterceptright <- cut_off_logFC_up
    df$yintercept <- -log(cut_off_pvalue,10)
    df$type <- name
    df$`-log10(P.value) `<- -log(df$p_val,10)
    
   
    #修改第一个列名
    colnames(df)[1] <-'gene'
    inter_val <- rbind(inter_val,df)
    print(i)    
    write.csv(df,paste0('/path/result/analy/new_Volcano/type/','Volcano_',i),row.names=F,quote=F)
}
write.csv(inter_val,'/path/result/analy/new_Volcano/Type_Volcano.csv',row.names=F,quote=F)

path = '/path/result/analy/new_Volcano/type_plot_data/'
files <- list.files(path)

for(i in 1:length(files)){
    
   type<-strsplit(files[i],split = ".csv")  
  filepath<-paste0("/path/result/analy/new_Volcano/type/",
                   as.character(type),".csv")
  Volcano_data<-read.csv(filepath,header=T)
    
    
    cut_off_logFC_down= -0.25
    cut_off_logFC_up=0.25
    cut_off_pvalue=0.05
    
    
      Volcano_plot<-ggplot(
    # 数据、映射、颜色
    Volcano_data, aes(x = avg_log2FC, y = -log(p_val,10), colour=change)) +
    geom_point(alpha=0.8, size=1.8) +
    scale_color_manual(values=c("#82B2D2", "grey","#FA7F6F"))+
    # 辅助线
    
    geom_vline(xintercept=c(cut_off_logFC_down,cut_off_logFC_up),lty=4,col="black",lwd=0.5) +
    geom_hline(yintercept = -log(cut_off_pvalue,10),lty=4,col="black",lwd=0.5) +
    # 坐标轴
    labs(x="log2(fold change)",
         y="-log10 (p-value)")+
    theme_bw()+
    # 图例
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank())
  
    
    print(i)
  Volcano_plotpath<-paste0("/path/result/analy/new_Volcano/type_plot/",
                           as.character(type),".png")
  Volcano_plotpath1<-paste0("/path/result/analy/new_Volcano/type_plot/",
                           as.character(type),".pdf") 
  ggsave(Volcano_plotpath, plot = Volcano_plot, width = 8, height = 6)
  ggsave(Volcano_plotpath1, plot = Volcano_plot, width = 8, height = 6)  
    
    }

Idents(sce) <- sce$subtype
unique(sce$subtype)

#根据细胞亚型subset细胞
for(i in unique(sce$subtype)){    
    obj <- subset(sce,idents=i)
    #obj@active.assay <- 'RNA'
    tmp <- try(marker <- FindMarkers(
        obj,ident.1 = "COVID-19",ident.2="Healthy",
                                     group.by='Status1',slot='data',assay = 'SCT',logfc.threshold =0,min.pct = 0.25,test.use='wilcox'))
# 设置logfc.threshold = log(2)参数过滤掉那些在两个不同组之间平均表达的差异倍数低于0.25的基因
 # 设置min.pct = 0.25参数过滤掉那些在50%以下细胞中检测到的基因
 
if('try-error' %in% class(tmp)){
  next
}
    print(i)
    write.csv(marker,paste0('/path/result/analy/new_subtype/',i,'.csv'))
    
        
}


path = '/path/result/analy/new_subtype'
files <- list.files(path)
inter_val <- c()
for(i in files){
    #获取细胞亚型名称
    name <- unlist(strsplit(i,split='.',fixed=T))[1]
    #读入文件
    df <- read.csv(paste0(path,'/',i),header=T)
        #上调基因降序排列，选取avg_log2FC前5%的基因。  
    up_gene_diff<-length(which(df$avg_log2FC>0.25))
    up_gene5.pct<-ceiling(up_gene_diff*0.2)
  
    down_gene_diff<-length(which(df$avg_log2FC<0.25))
    down_gene5.pct<-ceiling(down_gene_diff*0.2)
    
    gene_up <- df[which(df$avg_log2FC>0.25),]
    gene_up_5.pct <- gene_up[order(-gene_up$avg_log2FC),][1:up_gene5.pct,1]
    #下调基因升序排列，选取avg_log2FC前5%的基因 。 
    gene_down <- df[which(df$avg_log2FC<0.25),]
    gene_down_5.pct <- gene_down[order(gene_down$avg_log2FC),][1:down_gene5.pct,1]

    #把上下调的基因提取出来
    All_gene_5.pct <- df[which(df$X %in% c(gene_up_5.pct,gene_down_5.pct)),]

    #上下调基因的avg_log2FC分别提取出来
    FC_up <- All_gene_5.pct[which(All_gene_5.pct$avg_log2FC>0.25),]$avg_log2FC
    FC_down <- All_gene_5.pct[which(All_gene_5.pct$avg_log2FC<0.25),]$avg_log2FC

    cut_off_pvalue <- 0.05 #统计显著性
    cut_off_logFC_up <- 0.25   #差异倍数值
    cut_off_logFC_down <- -0.25
    
    #根据阈值参数，上调的基因设为up，下调基因设为down，无差异设为stable，保存在change列中。
    df$change <- ifelse(df$p_val < cut_off_pvalue &(df$avg_log2FC > cut_off_logFC_up | df$avg_log2FC < cut_off_logFC_down),
                       ifelse(df$avg_log2FC > cut_off_logFC_up,'Up','Down'),'Stable')
    #补全xinterceptleft、xinterceptright、yintercept、还有细胞亚型名
    df$xinterceptleft <- cut_off_logFC_down
    df$xinterceptright <- cut_off_logFC_up
    df$yintercept <- -log(cut_off_pvalue,10)
    df$subtype <- name
    df$`-log10(P.value) `<- -log(df$p_val,10)
    
   
    #修改第一个列名
    colnames(df)[1] <-'gene'
    inter_val <- rbind(inter_val,df)
    print(i)    
    write.csv(df,paste0('/path/result/analy/new_Volcano/subtype/','Volcano_',i),row.names=F,quote=F)
}
write.csv(inter_val,'/path/result/analy/new_Volcano/Subtype_Volcano.csv',row.names=F,quote=F)




path = '/path/result/analy/new_Volcano/subtype_plot_data/'
files <- list.files(path)

for(i in 1:length(files)){
    
   subtype<-strsplit(files[i],split = ".csv")  
  filepath<-paste0("/path/result/analy/new_Volcano/subtype/",
                   as.character(subtype),".csv")
  Volcano_data<-read.csv(filepath,header=T)
    
      cut_off_logFC_down= -0.25
    cut_off_logFC_up=0.25
    cut_off_pvalue=0.05
    
    
    
      Volcano_plot<-ggplot(
    # 数据、映射、颜色
    Volcano_data, aes(x = avg_log2FC, y = -log(p_val,10), colour=change)) +
    geom_point(alpha=0.8, size=1.8) +
    scale_color_manual(values=c("#82B2D2", "grey","#FA7F6F"))+
    # 辅助线
    geom_vline(xintercept=c(cut_off_logFC_down,cut_off_logFC_up),lty=4,col="black",lwd=0.5) +
    geom_hline(yintercept = -log(cut_off_pvalue,10),lty=4,col="black",lwd=0.5) +
    # 坐标轴
    labs(x="log2(fold change)",
         y="-log10 (p-value)")+
    theme_bw()+
    # 图例
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank())
  
    
    print(i)
  Volcano_plotpath<-paste0("/path/result/analy/new_Volcano/subtype_plot/",
                           as.character(subtype),".png")
  Volcano_plotpath1<-paste0("/path/result/analy/new_Volcano/subtype_plot/",
                           as.character(subtype),".pdf") 
  ggsave(Volcano_plotpath, plot = Volcano_plot, width = 8, height = 6)
  ggsave(Volcano_plotpath1, plot = Volcano_plot, width = 8, height = 6)  
    
    }

plot_size(12,8)
p1 <- DotPlot(sce, c('PPBP','GP9',"PF4",#Platelet  "CST3",
    "TCF4", "UGCG" ,"SERPINF1", "CCDC50",#pDCs
    "NCAM1","CD160","NKG7","CD247","GZMA","CCL5",#NK
    "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",#Monocytes   
    "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC  
    "CD3D","CD3E","CD3G", "CD8A","CD8B","GNLY", #CD8 T
    "IL7R","LTB","TRAC","CCR7","CD40LG","CD28","CD4",#CD4 T    
    "CD19","MS4A1","CD79A","CD79B","BLK","IGHM"#B cell 
                    ) ,dot.scale=6,# 缩放点的大小，类似于cex
     group.by='Type',assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
p1

ggsave("/path/result/PLOT/result_plot/DotPlot_type.pdf", plot = p1, width = 12, height = 8)
ggsave("/path/result/PLOT/result_plot/DotPlot_type.png", plot = p1, width = 12, height = 8)

plot_size(20,8)
DotPlot(sce,c(     

    'ITGA2B','GP9',"CD41",#Platelet  "CST3",
    "TCF4", "UGCG" ,"SERPINF1", "CCDC50",#pDCs
    "KLRD1","CD94",#NK  "CD16A",
    "CD16A","CD14","LYZ",#Monocytes   
    "LILRA4","CD83","CDF85g","CD1C","HLA-DQB2",#cDC   CD16A  LYZ
   "CD3E",# T

    "CD20","MS4A1"#B cell   


# "FOXP3","IL2RA","TRAC",#"CD4",             #CD4+ T cell
#"GZMA","FYN"            #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",      
   
    
    
                      
         ),group.by='Type',assay='RNA')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))





plot_size(20,8)
DotPlot(sce,c(     

    'PPBP','GP9',"PF4",#Platelet  "CST3",
    "TCF4", "UGCG" ,"SERPINF1", "CCDC50",#pDCs
    "NCAM1","CD160","NKG7","CD247","GZMA","CCL5",#NK
    "FCGR3A","CD14","FCN1","APOBEC3A","CD68","CSF1R",#Monocytes   
    "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC  
    "CD3D","CD3E","CD3G", "CD8A","CD8B","GNLY", #CD8 T
    "IL7R","LTB","TRAC","CCR7","CD40LG","CD28","CD4",#CD4 T    
    "CD19","MS4A1","CD79A","CD79B","BLK","IGHM","CD74"#B cell   


# "FOXP3","IL2RA","TRAC",#"CD4",             #CD4+ T cell
#"GZMA","FYN"            #CD8+ T cell  "CD3D","CD3E","CD3G","TRAC",      
   
    
    
                      
         ),group.by='Type',assay='RNA')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))

plot_size(12,15)
p1 <- DotPlot(sce, c(
      'PPBP','GP9',"PF4","NRGN",#Platelet  "CST3",
      "JCHAIN","IGKC","IGHA1","IGHG1",#Plasma  "IGKV4−1",
     "PTGDS","ITM2C","PLD4","IRF7","TCF4", "UGCG" ,#pDC
     "NCAM1","CD160","NKG7","CD247",#NK
     "FCGR3A" ,"CSF1R" ,"CDKN1C", "LST1" ,"C1QA", "WARS", "LILRB2", "PSAP" ,"AIF1" , #CD16    
     "CD14", "S100A8","S100A9","LYZ","S100A12","VCAN",#CD14 Monocyte
     "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC 
    "CD3D","CD3E","CD3G", "CD8A","CD8B",
    "CCR7", "LTB","PASK","IL32" ,#CD8+ naive：
     "FASLG","CD44","GNLY","CCL5","GZMH","CST7","KLRB1",#  CD8+ effector ：
     "IL7R" ,"SELL","LRRN3","GPR183",#CD4+ naive T cell   "CCR7", 
     "CD4","GZMK","TNFSF13B","TRAC","CD40LG","CD28","ICOS",#CD4+ effector T cell  "CCR7","SELL","LTB",
    
     'IGHD',"CD34","CD38","TCL1A" ,"BTG1","SIGLEC6",# Naive B    
     "CD19",  "CD79A", "MS4A1", "FGR", "FCRL5"# Memory "CD20",
    
                    ) ,dot.scale=4,# 缩放点的大小，类似于cex
     group.by='subtype',assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
p1

plot_size(20,8)
DotPlot(sce,c(    
    
      'PPBP','GP9',"PF4","NRGN",#Platelet  "CST3",
      "JCHAIN","IGKC","IGHA1","IGHG1",#Plasma  "IGKV4−1",
     "PTGDS","ITM2C","PLD4","IRF7","TCF4", "UGCG" ,#pDC
     "NCAM1","CD160","NKG7","CD247",#NK
     "FCGR3A" ,"CSF1R" ,"CDKN1C", "LST1" ,"C1QA", "WARS", "LILRB2", "PSAP" ,"AIF1" , #CD16    
     "CD14", "S100A8","S100A9","LYZ","S100A12","VCAN",#CD14 Monocyte
     "CD1C","CD1E","RNASE6","FCER1A","HLA-DQB2",#cDC 
    "CD3D","CD3E","CD3G", "CD8A","CD8B",
    "CCR7", "LTB","PASK","IL32" ,#CD8+ naive：
     "FASLG","CD44","GNLY","CCL5","GZMH","CST7","KLRB1",#  CD8+ effector ：
     "IL7R" ,"SELL","LRRN3","GPR183",#CD4+ naive T cell   "CCR7", 
     "CD4","GZMK","TNFSF13B","TRAC","CD40LG","CD28","ICOS",#CD4+ effector T cell  "CCR7","SELL","LTB",
    
     'IGHD',"CD34","CD38","TCL1A" ,"BTG1","SIGLEC6",# Naive B    
     "CD19",  "CD79A", "MS4A1", "FGR", "FCRL5","CD27"# Memory "CD20",
    
                      
         ),group.by='subtype',assay='RNA')+ 
theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))
#+scale_color_gradientn(colours = c("white",'#FFCC33')) #颜色

sample_all<-readRDS('/path/result/RDS/after_reduction.rds')

clusters_markers <- FindAllMarkers(sample_all,only.pos = F,test.use = 'wilcox',assay='SCT',slot='data',
               logfc.threshold = 0.5,min.pct = 0.25)

clusters_markers_top10 <- clusters_markers %>% group_by(cluster) %>% top_n(n=100,wt=avg_log2FC)
clusters_markers_top10

clusters_markers_tail10 <- clusters_markers %>% group_by(cluster) %>% top_n(n=100,wt= -avg_log2FC)
clusters_markers_tail10

write.csv(clusters_markers_top10,
         '/path/result/marker/FindAllMarkers_clusters_markers_top100.csv',
         row.names=F)
write.csv(clusters_markers_tail10,
         '/path/result/marker/FindAllMarkers_clusters_markers_tail100.csv',
         row.names=F)
write.csv(clusters_markers,
         '/path/result/marker/clusters_markers.csv',
         row.names=F)

clusters_markers_top5<-fread('/path/result/marker/FindAllMarkers_clusters_markers_top5.csv',header=T,data.table=F)

clusters_markers_top5

#超过30k细胞集有可能由于数据集过大而导致画不出热图，可以采取随机取样的方式将细胞数量downsample到30k以下
set.seed(42)
sample_all_subset <- subset(sample_all, downsample = 300)

plot_size(20,13)
DoHeatmap(sample_all_subset, features = clusters_markers_top5$gene, disp.min=-2.5, disp.max=2.5)

sample_all_subset

sample_all

plot_size(20,15)
maker_heatmap<-DoHeatmap(sample_all_subset, features = clusters_markers_top5$gene,label = F,assay = "SCT",slot = "scale.data")

maker_heatmap

mat <- GetAssayData(sce, slot = "counts")
mat <- log2(mat + 1)

gene_features <- clusters_markers_top5
cluster_info <- sort(sce$subtype)

cluster_info 

mat <- as.matrix(mat[clusters_markers_top5$gene, names(cluster_info)])

mat

write.table(mat,"/path/result/marker/result1_heatmap.txt",quote=F,sep="\t")

cluster_info

#cluster_rows= FALSE: 不作行聚类
#cluster_columns= FALSE: 不作列聚类
#show_column_names=FALSE: 不展示列名
#show_row_names=FALSE: 不展示行名，基因数目不多时候可以考虑设置为TRUE
plot_size(20,15)
Heatmap(mat,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T)

names(jdb_color_maps)<-c("B memory","B naive","CD4 T memory","CD4 T naive","CD8 T effector","CD8 T naive","cDC",
                         "Monocyte classical","Monocyte nonclassical","NK","pDC","Plasma","Platelet")

col <- jdb_color_maps[1:13]
col

names(col) <- levels(cluster_info)

levels(cluster_info)

names(col) 

#cluster_rows= FALSE: 不作行聚类
#cluster_columns= FALSE: 不作列聚类
#show_column_names=FALSE: 不展示列名
#show_row_names=FALSE: 不展示行名，基因数目不多时候可以考虑设置为TRUE
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info)

#HeatmapAnnotation对细胞的列进行注释，而rowAnnotation函数可以对行进行注释
top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # 设置填充色
                       labels = levels(cluster_info), 
                      labels_gp = gpar(cex = 0.8, col = "white"))) # 设置字体

top_anno 

Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno, # 在热图上边增加注释
        column_title = NULL ) # 不需要列标题

mark_gene <- c("IL7R","CCR7","IL7R","S100A4","CD14","LYZ","MS4A1","CD8A","FCGR3A","MS4A7","GNLY","NKG7","FCER1A", "CST3","PPBP")
gene_pos <- which(rownames(mat) %in% mark_gene)

row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                     labels = mark_gene))


Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL)

Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL,
        heatmap_legend_param = list(
          title = "log2(count 1)",
          title_position = "leftcenter-rot"
        ))

sce@meta.data

maker_gene<-list()
maker_gene$gene<-c( 
    "BANK1","MS4A1",#B cell 
    "IL7R","LDHB",#CD4 T   
    "CD8A","GNLY", #CD8 T    
    "CD1C","FCER1A",#cDC 
    "CD14","CD68",#Monocyte
    "NKG7","CD247",#NK
    "TCF4", "UGCG" ,#pDC
    'PPBP','GP9'#Platelet  "CST3"
 )   
maker_gene$type<-c("B cell",
                   "B cell",
                   "CD4 T",
                   "CD4 T",
                   "CD8 T",
                   "CD8 T",
                   "cDC",
                   "cDC",
                   "Monocyte",
                   "Monocyte",
                   "NK",
                   "NK",
                   "pDC",
                   "pDC",
                   "Platelet",
                   "Platelet")
maker_gene<-as.data.frame(maker_gene)


maker_gene$gene=as.character(maker_gene$gene)

maker_gene$gene

sce$Type=factor(sce$Type,levels = sort(unique(sce$Type)))
Idents(sce)="Type"

 sort(unique(sce$Type))

plot_size(10,35) 
VlnPlot(sce, features =maker_gene$gene, pt.size = 0, ncol = 1)+
#pt.size参数表示点的大小，一个点就是一个细胞，一般可以直接设置为0，即不显示点
  scale_x_discrete("")+
  theme(
    axis.text.x.bottom = element_blank()
  )

vln.df=as.data.frame(sce[["RNA"]]@data[maker_gene$gene,])
vln.df$gene=rownames(vln.df)
#vln.df<-rownames_to_column(vln.df,rownames(vln.df))
vln.df=reshape2::melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
vln.df

sce@meta.data["CB"]<-NA
sce@meta.data["CB"]<-rownames(sce@meta.data)
sce@meta.data

anno=sce@meta.data[,c("CB","Type")]
anno

vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = maker_gene$gene )#为了控制画图的基因顺序
vln.df

write.table(vln.df,"/path/code/result/data.txt",quote = F,row.names = F)

#col_celltypefine<-c("#FF77FF","#9955FF","#00DDDD","#00DD77","#DDFF77","#FFBB00","#EE7700","#FF0000")
#colors<- colorRampPalette(brewer.pal(n = 7, name ="Oranges"))(200)
#ColorCount <- length(unique(maker_gene$gene))
#getPalette <- colorRampPalette(brewer.pal(9,"Set1"))

color<-c("#984EA3","#E41A1C","#4DAF4A","#FF7F00","#377EB8","#FFFF33","#F781BF","#A65628")


vln.df

plot_size(8,7)
p1<-ggplot(vln.df,aes(Type,exp,fill=Type))+
geom_violin(linetype="blank",scale = "width",adjust = 1, trim = TRUE)+  
#linetype="blank" 把小提琴图的外缘轮廓去除  scale = "width"  设置相同宽度
  facet_grid(vln.df$gene~.,scales = "free_y")+  #分面
  #scale_fill_manual(values = getPalette(ColorCount))+
#scale_fill_manual(values =color)+
scale_fill_manual(values = colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))+
  #scale_x_discrete(limits = c("low", "med", "high"), #scale_x_discrete 修改 R 中的 x 轴名称，之前的名字
                   #labels= c("Low", "Medium", "High"))+ #修改后的名字
  scale_y_continuous("")+#改变y轴的刻度范围
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 90,hjust = 1,vjust = 1,color = 'black', family = 'Arial',size=12),
   panel.grid = element_blank(),
      # panel.grid.major = element_blank(),# 删除主要网格线
    #panel.grid.minor = element_blank(),# 删除次要网格线
      panel.background = element_blank(),#删除图像背景
    strip.text.y  = element_text(angle = 0,hjust = 0.5,vjust = 0.5,size=10),
      #strip.text.y = element_text(angle = 0),#strip.text.y调整分面文本的字体属性
      strip.background = element_blank(),#基因后面的背景 
       panel.border = element_blank(),#分割页面的外边框线
      plot.background = element_blank(),
      axis.ticks.x = element_blank(),  #删除轴刻度
       axis.ticks.y = element_blank(),  #删除轴刻度
      axis.text.y = element_blank(), #左面坐标刻度
     axis.title = element_blank(),#删除轴文本    
    legend.position = "left"
  )
p1
#panel.spacing.y调整分面的间距

#strip.background.y调整分面方块的背景
#expand_limits(0,0)快速设置在x和y轴在 (0,0) 处的截距项


# 画最右边的基因注释线段，用geom_segment来画，当然也可以用geom_tile()和geom_bar()来画
a=count(maker_gene$type)
# 用geom_segment，我们需要计算线段的起始点和终止点
# 计算每个基因group的成员数量，然后进行累加，每个累加值就是线段的终点，最后减去每个组的成员数量，就是线段的起始位点
a$yend=cumsum(a$freq)  #cumsum () R语言中的函数用于计算作为参数传递的向量的累积和
a$ystart=a$yend-a$freq
# 计算y轴的text的位置，就是group名称的位置，应该是每个线段的中点
a$label_position=a$ystart+(a$yend-a$ystart)/2
colnames(a)[1]<-"type"
a
p2=ggplot(a,aes(x=1,y=ystart,color=type))+
  geom_segment(aes(xend=1,yend=yend),size=2)+
  scale_y_continuous(position = "right",
                     breaks = a$label_position,
                     labels = a$type,expand = c(0,0))+
 # scale_color_manual(values = getPalette(ColorCount))+
scale_color_manual(values = colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))+
  scale_x_continuous(expand = c(0,0))+
  facet_grid(type~.,scales = 'free')+
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0,hjust = 1,vjust = 1,color = 'black', family = 'Arial',size=12),
        axis.ticks = element_blank(),
        legend.position = "none")
p2

p1+p2+plot_spacer()+plot_layout(ncol = 2, widths  = c(3, 0.03),heights = c(3,0.03))
#widths 和 heights 参数提供了各个矩形作图区域的长和宽的比例
#plot_spacer()     空白的占位图
dev.copy2pdf( file="/path/result/PLOT/result_plot/violin_type_maker.pdf",width = 11, height = 8)

a

plot_size(5,3)
#geom_tile()完成细胞注释
p3<-ggplot(a,
           aes(x=type,y=freq,fill=type))+
           geom_tile()+
 # scale_fill_brewer(palette = "Set1")+
scale_fill_manual(values = colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  #facet_grid(.~type,scales = "free",space = "free",switch = 'x')
  facet_grid(.~type,scales = 'free',space = "free",switch = 'x')+
  theme(panel.background = element_blank(),
        strip.text.x = element_text(angle = 0,hjust = 0.5,vjust =0.5,size=13),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
p3

#pdf("/path/result/PLOT/result_plot/violin_type_maker.pdf", width = 10, height = 8)


plot_size(11,7)
a=p1+p2+p3+plot_spacer()+plot_layout(ncol = 2, widths  = c(3, .03),heights = c(3,.03))
a
#ggsave("/path/result/PLOT/result_plot/violin_type_maker.pdf", plot = a, width = 11, height = 8)
#ggsave("/path/result/PLOT/result_plot/violin_type_maker.png", plot = a, width = 11, height = 8)
#dev.copy2pdf(file = "/path/result/PLOT/result_plot/violin_type_maker.png", paper = "a4r")
#dev.copy2pdf(file = "/path/result/PLOT/result_plot/violin_type_maker.pdf", paper = "a4r")
#dev.off()

dev.copy2pdf( file="/path/result/PLOT/result_plot/violin_type_maker.pdf",width = 11, height = 8)

a
dev.copy2pdf()

getwd()

'CCR7','NELL2','CD55','KLF2')

sce<-readRDS("/path/result/RDS/all_subtype_annotated.rds")

sce@meta.data

sce

options(repr.plot.width = 10, repr.plot.height = 4)
p1<-VlnPlot(sce,features = c("STAT3"),group.by = "Type",raster = FALSE,pt.size = 0,
        cols = c("#984EA3","#E41A1C","#4DAF4A","#FF7F00","#377EB8","#FFFF33","#F781BF","#A65628"))
p1
ggsave("/path/result/PLOT/result_plot/violin_subtype_STAT3.pdf", plot = p1, width = 13, height = 5)

options(repr.plot.width = 10, repr.plot.height = 4)
p1<-VlnPlot(sce,features = c("STAT3"),group.by = "Status",raster = FALSE,pt.size = 0,
        cols = c(  "#B6B6B7","#BDD7E7","#6BAED6","#08519C","#3182BD"))
p1
ggsave("/path/result/PLOT/result_plot/violin_Status_STAT3.pdf", plot = p1, width = 13, height = 5)

options(repr.plot.width = 10, repr.plot.height = 4)

levels(sce$Type)
levels=c("Monocyte","B","CD4","CD8","cDC",
                                              "NK","pDC",
                                              "Platelet")
p2<-VlnPlot(sce,features = c("IL10RB"),group.by = "Type",raster = FALSE,pt.size = 0,
        cols = c('#FF850E','#1F83B4','#24A052','#95AF33','#FFBC46','#CD4771','#AC4BB2','#6F63BB'))
p2
ggsave("/path/result/PLOT/result_plot/violin_subtype_IL10RB.pdf", plot = p2, width = 13, height = 5)







subtype<-as.data.frame(sce@meta.data$subtype)
subtype

sce@meta.data

 'PPBP','GP9',#Platelet  
      "JCHAIN","IGHA1"#Plasma  
     "ITM2C","PLD4",#pDC
     "NKG7","CD247",#NK
     "FCGR3A" ,"CSF1R" ,"CDKN1C" #CD16    
     "CD14", "VCAN",#CD14 Monocyte
     "CD1C","FCER1A",#cDC 
   "CD8B",
    "CCR7", #CD8+ naive：
    "GNLY","GZMH",#  CD8+ effector ：
     "IL7R" ,"SELL",#CD4+ naive T cell   "CCR7", 
     "GZMK","TRAC",#CD4+ effector T cell  "CCR7","SELL","LTB",
    
     'IGHD',"TCL1A" ,# Naive B    
     "CD79A", "MS4A1",# Memory "CD20",

maker_gene<-list()
maker_gene$gene<-c(
       
      "CD79A", "MS4A1",# Memory "CD20",
     'IGHD',"TCL1A" ,# Naive B  
     "JCHAIN","IGHA1",#Plasma 
     "GZMK","TRAC",#CD4+ memory T cell  "CCR7","SELL","LTB",
    "IL7R" ,"SELL",#CD4+ naive T cell   "CCR7", 
    "GNLY","GZMH",#  CD8+ effector ：
    "CD8B","TCF7",#CD8+ naive：
     "CD1C","FCER1A",#cDC 
    "CD14", "VCAN",#Monocyte classical
     "FCGR3A" ,"CDKN1C", #CD16  
    "NKG7", "CD247",#NK
      "ITM2C","PLD4",#pDC
   
       'PPBP','GP9'#Platelet 
 )   
maker_gene$subtype<-c(
        
                  "B memory",
                  "B memory",
                  "B naive" ,
                  "B naive" ,
                  "B_Plasma",
                  "B_Plasma",
                  "CD4 memory",
                  "CD4 memory",
                  "CD4 naive",
                  "CD4 naive",
                  "CD8 effector",
                  "CD8 effector",
                  "CD8 naive",
                  "CD8 naive",
                  "cDC",
                  "cDC",
                  "Monocyte classical",
                  "Monocyte classical",
                  "Monocyte nonclassical",
                  "Monocyte nonclassical",
                  "NK",
                  "NK",
                  "pDC",
                  "pDC",           
                  "Platelet",
                  "Platelet"
                  )
maker_gene<-as.data.frame(maker_gene)



maker_gene<-list()
maker_gene$gene<-c(
       
      "CD79A", "MS4A1",# Memory "CD20",
     'FCER2',"IGHD" ,# Naive B  
     "JCHAIN","IGHA1",#Plasma 
     "CD3G","IL7R",#CD4+ memory T cell  "CCR7","SELL","LTB",
    "TCF7" ,"SELL",#CD4+ naive T cell   "CCR7", 
    "GNLY","GZMH",#  CD8+ effector ：
    "CD8B","CCR7",#CD8+ naive：
     "CD1C","FCER1A",#cDC 
    "CD14", "VCAN",#Monocyte classical
     "FCGR3A" ,"CDKN1C", #CD16  
    "NKG7", "CD247",#NK
      "ITM2C","PLD4",#pDC
   
       'PPBP','GP9'#Platelet 
 )   
maker_gene$subtype<-c(
        
                  "B memory",
                  "B memory",
                  "B naive" ,
                  "B naive" ,
                  "B_Plasma",
                  "B_Plasma",
                  "CD4 memory",
                  "CD4 memory",
                  "CD4 naive",
                  "CD4 naive",
                  "CD8 effector",
                  "CD8 effector",
                  "CD8 naive",
                  "CD8 naive",
                  "cDC",
                  "cDC",
                  "Monocyte classical",
                  "Monocyte classical",
                  "Monocyte nonclassical",
                  "Monocyte nonclassical",
                  "NK",
                  "NK",
                  "pDC",
                  "pDC",           
                  "Platelet",
                  "Platelet"
                  )
maker_gene<-as.data.frame(maker_gene)


maker_gene$gene=as.character(maker_gene$gene)
maker_gene$gene

sce$subtype=factor(sce$subtype,levels = unique(maker_gene$subtype))
Idents(sce)="subtype"

unique(maker_gene$subtype)

plot_size(10,70) 
VlnPlot(sce, features =maker_gene$gene, pt.size = 0, ncol = 1)+
#pt.size参数表示点的大小，一个点就是一个细胞，一般可以直接设置为0，即不显示点
  scale_x_discrete("")+
  theme(
    axis.text.x.bottom = element_blank()
  )

vln.df=as.data.frame(sce[["RNA"]]@data[maker_gene$gene,])
vln.df$gene=rownames(vln.df)
#vln.df<-rownames_to_column(vln.df,rownames(vln.df))
vln.df=reshape2::melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
vln.df

sce@meta.data["CB"]<-NA
sce@meta.data["CB"]<-rownames(sce@meta.data)
sce@meta.data["subtype"]<-subtype[,1]
sce@meta.data

anno=sce@meta.data[,c("CB","subtype")]
anno

vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = maker_gene$gene )#为了控制画图的基因顺序
vln.df$subtype[which(vln.df$subtype=="Plasma")]<-"B_Plasma"
vln.df

#col_celltypefine<-c("#FF77FF","#9955FF","#00DDDD","#00DD77","#DDFF77","#FFBB00","#EE7700","#FF0000")
#colors<- colorRampPalette(brewer.pal(n = 7, name ="Oranges"))(200)
ColorCount <- length(unique(maker_gene$gene))
getPalette <- colorRampPalette(brewer.pal(9,"Set1"))


COLOR<-c("#C30E6F","#F3AD5E","#30983C","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#0871B4" ,"#9CC8E0" )  


#COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
#         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ) 

plot_size(10,13)
p1<-ggplot(vln.df,aes(subtype,exp,fill=subtype))+
geom_violin(linetype="blank",scale = "width",adjust = 1, trim = TRUE)+  
#linetype="blank" 把小提琴图的外缘轮廓去除  scale = "width"  设置相同宽度
  facet_grid(vln.df$gene~.,scales = "free_y")+  #分面
  scale_fill_manual(values =COLOR)+
  #scale_x_discrete(limits = c("low", "med", "high"), #scale_x_discrete 修改 R 中的 x 轴名称，之前的名字
                   #labels= c("Low", "Medium", "High"))+ #修改后的名字
  scale_y_continuous("")+#改变y轴的刻度范围
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 90,hjust = 1,vjust = 1,color = 'black', family = 'Arial',size=12),
   panel.grid = element_blank(),
      # panel.grid.major = element_blank(),# 删除主要网格线
    #panel.grid.minor = element_blank(),# 删除次要网格线
      panel.background = element_blank(),#删除图像背景
    strip.text.y  = element_text(angle = 0,hjust = 0.5,vjust = 0.5,size=10),
      #strip.text.y = element_text(angle = 0),#strip.text.y调整分面文本的字体属性
      strip.background = element_blank(),#基因后面的背景 
       panel.border = element_blank(),#分割页面的外边框线
      plot.background = element_blank(),
      axis.ticks.x = element_blank(),  #删除轴刻度
       axis.ticks.y = element_blank(),  #删除轴刻度
      axis.text.y = element_blank(), #左面坐标刻度
     axis.title = element_blank(),#删除轴文本    
    legend.position = "left"
  )
p1
#panel.spacing.y调整分面的间距

#strip.background.y调整分面方块的背景
#expand_limits(0,0)快速设置在x和y轴在 (0,0) 处的截距项


# 画最右边的基因注释线段，用geom_segment来画，当然也可以用geom_tile()和geom_bar()来画
a=count(maker_gene$subtype)
# 用geom_segment，我们需要计算线段的起始点和终止点
# 计算每个基因group的成员数量，然后进行累加，每个累加值就是线段的终点，最后减去每个组的成员数量，就是线段的起始位点
a$yend=cumsum(a$freq)  #cumsum () R语言中的函数用于计算作为参数传递的向量的累积和
a$ystart=a$yend-a$freq
# 计算y轴的text的位置，就是group名称的位置，应该是每个线段的中点
a$label_position=a$ystart+(a$yend-a$ystart)/2
colnames(a)[1]<-"subtype"
a
p2=ggplot(a,aes(x=1,y=ystart,color=subtype))+
  geom_segment(aes(xend=1,yend=yend),size=2)+
  scale_y_continuous(position = "right",
                     breaks = a$label_position,
                     labels = a$subtype,expand = c(0,0))+
  scale_color_manual(values = COLOR)+
  scale_x_continuous(expand = c(0,0))+
  facet_grid(subtype~.,scales = 'free')+
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0,hjust = 1,vjust = 1,color = 'black', family = 'Arial',size=12),
        axis.ticks = element_blank(),
        legend.position = "none")
p2

p1+p2+plot_spacer()+plot_layout(ncol = 2, widths  = c(3, 0.03),heights = c(3,0.03))
#widths 和 heights 参数提供了各个矩形作图区域的长和宽的比例
#plot_spacer()     空白的占位图

a$type<-c("B cell","B cell","B cell","CD4 T","CD4 T","CD8 T","CD8 T","cDC","Monocyte","Monocyte","NK", "pDC","Platelet")

a



#COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
#         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ,"#9CC8E0" )   

plot_size(10,5)
#geom_tile()完成细胞注释
p3<-ggplot(a,
           aes(x=subtype,y=freq,fill=type))+
           geom_tile()+
 # scale_fill_brewer(palette = "Set1")+
  #scale_fill_manual(values = COLOR)+
scale_fill_manual(values = colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  #facet_grid(.~subtype,scales = "free",space = "free",switch = 'x')+
  facet_grid(.~type,scales = 'free',space = "free",switch = 'x')+
  theme(panel.background = element_blank(),
        strip.text.x = element_text(angle = 0,hjust = 0.5,vjust =0.1,size=13),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
p3



plot_size(13,10)
a=p1+p2+p3+plot_spacer()+plot_layout(ncol = 2, widths  = c(3, .03),heights = c(3,.03))
a
#ggsave("/path/result/PLOT/result_plot/violin_type_maker.pdf", plot = a, width = 10, height = 8)
ggsave("/path/result/PLOT/result_plot/violin_subtype_maker.png", plot = a, width = 13, height = 10)









plot_size(15,5)  #"CD247"
FeaturePlot(sce,feature=c("SON","IL10RB","DNMT1"


                 ),reduction = 'umap',ncol=4)

plot_size(15,5)  #"CD247"
FeaturePlot(sce,feature=c("CD14","FCGR3A","FCGR3B","FCG3","FCGR3","IGFR3"


                 ),reduction = 'umap',ncol=4)

MALAT1
MTRNR2L12
S100A8
S100A9
MALAT1
IGHG2
IGLC2
IGHA1
FTL
MALAT1


plot_size(15,10)  #"CD247"
FeaturePlot(sce,feature=c("MALAT1",
"MTRNR2L12",
"S100A8",
"S100A9",
"MALAT1",
"IGHG2",
"IGLC2",
"IGHA1",
"FTL",
"MALAT1",
                          "SAT1",

"MT2A",
"CCL3",

"FTH1",
"S100A12",
"MALAT1",
"MT2A"


                 ),reduction = 'umap',ncol=4)



plot_size(15,10)  #"CD247"
FeaturePlot(sce,feature=c("PF4" ,"TCF4" ,"CD247" ,"FCN1" ,"CD1C" ,"GNLY", "IL7R" ,"MS4A1"
                 ),reduction = 'umap',ncol=4)

sce@meta.data

table.cluster<-table(sce@active.ident)
##计算每个样本的细胞数
table.sample<-table(sce$orig.ident1)
##计算细胞比例
prop_allcell<-prop.table(table(sce@active.ident))
cell_num<-table(sce@active.ident, sce$orig.ident1)
cell_numtable<-as.data.frame(table(sce@active.ident, sce$orig.ident1))
colnames(cell_numtable)<-c("Celltype","Sample","Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_num.csv',row.names = F)


cell_numtable

 COLOR<-c("#24A052","#FF850E","#95AF33","#1F83B4","#FFBC46","#CD4771","#6F63BB","#AC4BB2") 

#COLOR<-c("#24A052","#FF850E","#95AF33","#1F83B4","#FFBC46","#CD4771","#6F63BB","#AC4BB2") 
plot_size(15,8)
##细胞比例图##
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=Celltype))+
  geom_bar(stat="identity",position="fill",width = .9)+
  ggtitle("")+
  theme_bw()+
  scale_x_discrete("")+
  scale_y_continuous("Total cell proportion",expand = c(0,0),labels = scales::label_percent(),position = "left")+  
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  #scale_fill_brewer(palette = "Set1")+
#scale_fill_manual(values = colorRampPalette(paletteer_d("ggthemes::Classic_Cyclic"))(8))+
scale_fill_manual(values = COLOR)+
   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size=11,color ="black"),
         axis.text.y=element_text(size=11,color ="black"),
         panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")) +      #让横轴上的标签倾斜45度
  xlab("Sample")    
# +coord_flip() #让条形图横过来,本来向上的条形图向右

cell_plot 

ggsave("/path/result/PLOT/result_plot/bar_type_sample.pdf", plot = cell_plot , width = 15, height = 8)
ggsave("/path/result/PLOT/result_plot/bar_type_sample.png", plot = cell_plot , width = 15, height = 8)

ggplot(bar.df,aes(x=name))+
  geom_bar(aes(fill=Celltype),position = "stack",width = .7)+
  scale_x_discrete("",position = "top")+  #x轴在上
  scale_fill_manual(values = color_cluster)+
  scale_y_reverse("Total cell number",expand = c(0,0),position = "right",limits=c(5000,0))+  #y轴在右
  theme(
    panel.grid = element_blank(),
    legend.position = "none", # 不加图例
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+
  coord_flip()

sce<-readRDS("/path/result/RDS/all_subtype_annotated.rds")

sce@meta.data

table.cluster<-table(sce@meta.data$subtype)
##计算每个样本的细胞数
table.sample<-table(sce$orig.ident1)
##计算细胞比例
prop_allcell<-prop.table(table(sce@meta.data$subtype))
cell_num<-table(sce@meta.data$subtype, sce$orig.ident1)
cell_numtable<-as.data.frame(table(sce@meta.data$subtype, sce$orig.ident1))
colnames(cell_numtable)<-c("subtype","Sample","Number")
##输出cell_numtable
#write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/subcell_num.csv',row.names = F)
cell_numtable


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ,"#9CC8E0" ) 
plot_size(15,8)
##细胞比例图##
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=subtype))+
  geom_bar(stat="identity",position="fill",width = .9)+
  ggtitle("")+
  theme_bw()+
  scale_x_discrete("")+
  scale_y_continuous("Total cell proportion",expand = c(0,0),labels = scales::label_percent(),position = "left")+  
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = COLOR)+
   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size=11,color ="black"),
         panel.grid = element_blank(),
         axis.text.y=element_text(size=11,color ="black"),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")) +      #让横轴上的标签倾斜45度
  xlab("Sample")    
# +coord_flip() #让条形图横过来,本来向上的条形图向右

cell_plot 

ggsave("/path/result/PLOT/result_plot/bar_subtype_sample.pdf", plot = cell_plot , width = 15, height = 8)
ggsave("/path/result/PLOT/result_plot/bar_subtype_sample.png", plot = cell_plot , width = 15, height = 8)





table.cluster<-table(sce@meta.data$Status)
##计算每个样本的细胞数
table.sample<-table(sce$seurat_clusters)
##计算细胞比例
prop_allcell<-prop.table(table(sce@meta.data$Status))
cell_num<-table(sce@meta.data$Status, sce$seurat_clusters)
cell_numtable<-as.data.frame(table(sce@meta.data$Status, sce$seurat_clusters))
colnames(cell_numtable)<-c("Status","clusters","Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_seurat_clusters_num.csv',row.names = F)


color<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99")

plot_size(10,5)
cell_plot<-ggplot(cell_numtable,aes(clusters,Number,fill=Status))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
 # scale_fill_brewer(palette = "Set1")+
#scale_fill_brewer(brewer.pal(5,"Paired"))+
  scale_fill_manual(values = color)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +      #让横轴上的标签倾斜45度
  xlab("clusters")    
cell_plot



ggsave("/path/result/PLOT/result_plot/bar_type_clusters.pdf", plot = cell_plot , width = 10, height = 5)
ggsave("/path/result/PLOT/result_plot/bar_type_clusters.png", plot = cell_plot , width = 10, height = 5)

color<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99")

plot_size(10,5)
cell_plot<-ggplot(cell_numtable,aes(clusters,Number,fill=Status))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
 # scale_fill_brewer(palette = "Set1")+
scale_fill_brewer(brewer.pal(5,"Paired"))+
  #scale_fill_manual(values = color)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +      #让横轴上的标签倾斜45度
  xlab("clusters")    
cell_plot



ggsave("/path/result/PLOT/result_plot/blue_bar_type_clusters.pdf", plot = cell_plot , width = 10, height = 5)
ggsave("/path/result/PLOT/result_plot/blue_bar_type_clusters.png", plot = cell_plot , width = 10, height = 5)

ggsave("/path/result/analy/CELL/cell_numtable/clusters.pdf", plot = cell_plot, width = 10, height = 6)
ggsave("/path/result/analy/CELL/cell_numtable/clusters.png", plot = cell_plot, width = 10, height = 6)









table.cluster<-table(sce@meta.data$orig.ident1)
##计算每个样本的细胞数
table.sample<-table(sce$seurat_clusters)
##计算细胞比例
prop_allcell<-prop.table(table(sce@meta.data$orig.ident1))
cell_num<-table(sce@meta.data$orig.ident1, sce$seurat_clusters)
cell_numtable<-as.data.frame(table(sce@meta.data$orig.ident1, sce$seurat_clusters))
colnames(cell_numtable)<-c("orig.ident1","clusters","Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_seurat_clusters_orig.ident1_num.csv',row.names = F)


plot_size(10,5)
cell_plot<-ggplot(cell_numtable,aes(clusters,Number,fill=orig.ident1))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
 # scale_fill_brewer(palette = "Set1")+
   theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +      #让横轴上的标签倾斜45度
  xlab("clusters")    
cell_plot

ggsave("/path/result/PLOT/result_plot/bar_clusters_sample.pdf", plot = cell_plot , width = 10, height = 5)
ggsave("/path/result/PLOT/result_plot/bar_clusters_sample.png", plot = cell_plot , width = 10, height = 5)

table.cluster<-table(sce@active.ident)
##计算每个样本的细胞数
table.sample<-table(sce$Status)
##计算细胞比例
prop_allcell<-prop.table(table(sce@active.ident))
cell_num<-table(sce@active.ident, sce$Status)
cell_numtable<-as.data.frame(table(sce@active.ident, sce$Status))
colnames(cell_numtable)<-c("Celltype","Sample","Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_Status_num.csv',row.names = F)
cell_numtable


 COLOR<-c("#24A052","#FF850E","#95AF33","#1F83B4","#FFBC46","#CD4771","#6F63BB","#AC4BB2") 
plot_size(5,8)
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  #scale_fill_brewer(palette = "Set1")+
scale_fill_manual(values = COLOR)+
   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size=12,color ="black"),
         axis.text.y=element_text(size=11,color ="black"),) +      #让横轴上的标签倾斜45度
  xlab("Sample")    
cell_plot

ggsave("/path/result/PLOT/result_plot/bar_type_status.pdf", plot = cell_plot , width = 5, height = 8)
ggsave("/path/result/PLOT/result_plot/bar_type_status.png", plot = cell_plot , width = 5, height = 8)

plot_size(8,5)
##细胞比例图##
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=Celltype))+
  geom_bar(stat="identity",position="fill",width = .7)+
  ggtitle("")+
  theme_bw()+
  scale_x_discrete("")+
  scale_y_continuous("Total cell proportion",expand = c(0,0),labels = scales::label_percent(),position = "right",)+  
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_brewer(palette = "Set1")+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size=15,color = 'black', family = 'Times New Roman'),
         axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size=15,color = 'black', family = 'Times New Roman'),
         panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black"),
        legend.position = "right",
        legend.justification = c(0.5,0.5),
       legend.text = element_text(family = 'Times New Roman')) +      #让横轴上的标签倾斜45度
  xlab("Sample") +   
 coord_flip() #让条形图横过来,本来向上的条形图向右

cell_plot 

metadata <-sce@meta.data 
metadata <- mutate(metadata,name=factor(metadata$Status, levels=c("asymptomatic case","Healthy","Moderate","Recover","Severe")))
text.df <- as.data.frame(table(metadata$name))
color_cluster=c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999')
names(color_cluster)=c("B","CD4 T","CD8 T","cDC","Monocyte","NK","pDC","Platelet")
metadata

#col_celltypefine<-c("#FF77FF","#9955FF","#00DDDD","#00DD77","#DDFF77","#FFBB00","#EE7700","#FF0000")
colors<- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(40)
#color_cluster=c("#99a9cc","#f89e81","#acd485","#dd9bc5","#f6d573","#84c7b3")
color_cluster=c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999')
color_cluster

plot_size(20,5)
cell_plot1<-ggplot(metadata,aes(x=name,fill=Type))+
  geom_bar(position="stack",width = .7)+
  ggtitle("")+
  theme_bw()+
scale_x_discrete(position = "top")+  #x轴在上
  scale_fill_manual(values = color_cluster)+
scale_y_reverse(expand = c(0,0),position = "right",limits=c(62000,0))+  #y轴在右
  theme( axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size=15,color = 'black', family = 'Times New Roman'),
        #让横轴上的标签倾斜0度
         axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size=15,color = 'black', family = 'Times New Roman'),
      axis.ticks.length=unit(0.1,'cm'),
         panel.grid = element_blank(),
         legend.position = "none")+
      # legend.justification = c(0.5,0.5),
      # legend.text = element_text(family = 'Times New Roman'))+
  guides(fill=guide_legend(title=NULL))+
# scale_fill_brewer(palette = "Set1")+
 #scale_fill_discrete(=c("#E41A1C",'#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'),
    #name="Experimental\nCondition",  #修改图例名称
    #labels = c("1","2","3","4","5","6","7","8"),  #修改图例标签
    #breaks = c("B","CD4 T","CD8 T","cDC","Monocyte","NK","pDC","Platelet"))+#图里之前的标签
  xlab("")   + ylab("")+
  coord_flip() 
cell_plot1+cell_plot

table.cluster<-table(sce@active.ident)
##计算每个样本的细胞数
table.sample<-table(sce$Status)
##计算细胞比例
prop_allcell<-prop.table(table(sce@active.ident))
cell_num<-table(sce@active.ident, sce$Status)
cell_numtable<-as.data.frame(table(sce@active.ident, sce$Status))
cell_numtable$table_Status<-rep(table(sce@meta.data$Status), each = 8,len = 40)
cell_numtable$table.Type<-table(sce@meta.data$Type)
colnames(cell_numtable)<-c("Celltype","Sample","Number","Status_Total_Number","Type_Total_Number")
##输出cell_numtable
#write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_Status_num.csv',row.names = F)
cell_numtable$Type_Total_Number=factor(cell_numtable$Number)
cell_numtable

table.cluster<-table(sce@meta.data$subtype)
##计算每个样本的细胞数
table.sample<-table(sce$Status)
##计算细胞比例
prop_allcell<-prop.table(table(sce@meta.data$subtype))
cell_num<-table(sce@meta.data$subtype, sce$Status)
cell_numtable<-as.data.frame(table(sce@meta.data$subtype, sce$Status))
colnames(cell_numtable)<-c("subtype","Sample","Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/subcell_num.csv',row.names = F)
cell_numtable


COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ,"#9CC8E0" ) 
plot_size(5,8)
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=subtype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = COLOR)+
   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size=12,color ="black"),
          axis.text.y=element_text(size=11,color ="black"))+      #让横轴上的标签倾斜45度
  xlab("Sample")    
cell_plot

ggsave("/path/result/PLOT/result_plot/bar_subtype_5status.pdf", plot = cell_plot , width = 5, height = 8)
ggsave("/path/result/PLOT/result_plot/bar_subtype_5status.png", plot = cell_plot , width = 5, height = 8)



COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ,"#9CC8E0" ) 
plot_size(8,5)
##细胞比例图##
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=subtype))+
  geom_bar(stat="identity",position="fill",width = .7)+
  ggtitle("")+
  theme_bw()+
  scale_x_discrete("")+
  scale_y_continuous("Total cell proportion",expand = c(0,0),labels = scales::label_percent(),position = "right",)+  
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
   scale_fill_manual(values = COLOR)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size=15,color = 'black', family = 'Times New Roman'),
         axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size=15,color = 'black', family = 'Times New Roman'),
         panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black"),
        legend.position = "right",
        legend.justification = c(0.5,0.5),
       legend.text = element_text(family = 'Times New Roman')) +      #让横轴上的标签倾斜45度
  xlab("Sample") +   
 coord_flip() #让条形图横过来,本来向上的条形图向右

cell_plot 

COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ,"#9CC8E0" ) 

metadata <-sce@meta.data 
metadata <- mutate(metadata,name=factor(metadata$Status, levels=c("asymptomatic case","Healthy","Moderate","Recover","Severe")))
text.df <- as.data.frame(table(metadata$name))
#color_cluster=c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999')
COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ,"#9CC8E0" ) 
names(COLOR)=c(

    'B memory','B naive','CD4 T memory','CD4 T naive','CD8 T effector','CD8 T naive','cDC','Monocyte classical','Monocyte nonclassical','NK','pDC','Plasma','Platelet'

)
metadata

#col_celltypefine<-c("#FF77FF","#9955FF","#00DDDD","#00DD77","#DDFF77","#FFBB00","#EE7700","#FF0000")
colors<- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(40)
#color_cluster=c("#99a9cc","#f89e81","#acd485","#dd9bc5","#f6d573","#84c7b3")
color_cluster=c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999')
color_cluster

plot_size(20,5)
cell_plot1<-ggplot(metadata,aes(x=name,fill=subtype))+
  geom_bar(position="stack",width = .7)+
  ggtitle("")+
  theme_bw()+
scale_x_discrete(position = "top")+  #x轴在上
    scale_fill_manual(values = COLOR)+
scale_y_reverse(expand = c(0,0),position = "right",limits=c(62000,0))+  #y轴在右
  theme( axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size=15,color = 'black', family = 'Times New Roman'),
        #让横轴上的标签倾斜0度
         axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size=15,color = 'black', family = 'Times New Roman'),
      axis.ticks.length=unit(0.1,'cm'),
         panel.grid = element_blank(),
         legend.position = "none")+
      # legend.justification = c(0.5,0.5),
      # legend.text = element_text(family = 'Times New Roman'))+
  guides(fill=guide_legend(title=NULL))+
# scale_fill_brewer(palette = "Set1")+
 #scale_fill_discrete(=c("#E41A1C",'#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF'),
    #name="Experimental\nCondition",  #修改图例名称
    #labels = c("1","2","3","4","5","6","7","8"),  #修改图例标签
    #breaks = c("B","CD4 T","CD8 T","cDC","Monocyte","NK","pDC","Platelet"))+#图里之前的标签
  xlab("")   + ylab("")+
  coord_flip() 
cell_plot1+cell_plot

cell_numtable<-as.data.frame(table(sce@active.ident, sce$Status))
cell_numtable$table_Status<-rep(table(sce@meta.data$Status), each = 8,len = 40)
cell_numtable$table.Type<-table(sce@meta.data$Type)
colnames(cell_numtable)<-c("Celltype","Sample","Number","Status_Total_Number","Type_Total_Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_Status_num.csv',row.names = F)
cell_numtable$Type_Total_Number=factor(cell_numtable$Number)
cell_numtable

table.cluster<-table(sce@meta.data$subtype)
##计算每个样本的细胞数
table.sample<-table(sce$Status)
##计算细胞比例
prop_allcell<-prop.table(table(sce@meta.data$subtype))
cell_num<-table(sce@meta.data$subtype, sce$Status)
cell_numtable<-as.data.frame(table(sce@meta.data$subtype, sce$Status))
cell_numtable$table_Status<-rep(table(sce@meta.data$Status), each = 13,len = 65)
cell_numtable$table.subtype<-table(sce@meta.data$subtype)
colnames(cell_numtable)<-c("subtype","Sample","Number","Status_Total_Number","Type_Total_Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_Status_num.csv',row.names = F)
cell_numtable$Type_Total_Number=factor(cell_numtable$Number)
cell_numtable











table.cluster<-table(sce@active.ident)
##计算每个样本的细胞数
table.sample<-table(sce$Status1)
##计算细胞比例
prop_allcell<-prop.table(table(sce@active.ident))
cell_num<-table(sce@active.ident, sce$Status1)
cell_numtable<-as.data.frame(table(sce@active.ident, sce$Status1))
colnames(cell_numtable)<-c("Celltype","Sample","Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_Status1_num.csv',row.names = F)



 COLOR<-c("#24A052","#FF850E","#95AF33","#1F83B4","#FFBC46","#CD4771","#6F63BB","#AC4BB2") 

plot_size(3.5,8)
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  #scale_fill_brewer(palette = "Set1")+
scale_fill_manual(values = COLOR)+
   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size=12,color ="black"),
      axis.text.y=element_text(size=11,color ="black")) +      #让横轴上的标签倾斜45度
  xlab("Sample")  +
theme(legend.position = "right")
#scale_fill_manual(values= c("#E41A1C",'#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999'),
  # name="Experimental\nCondition",  #修改图例名称
   # labels = c("1","2","3","4","5","6","7","8"),  #修改图例标签
   # breaks = c("B","CD4 T","CD8 T","cDC","Monocyte","NK","pDC","Platelet")) #图里之前的标签
 
cell_plot

ggsave("/path/result/PLOT/result_plot/bar_type_covid_hc.pdf", plot = cell_plot , width = 3.5, height = 8)
ggsave("/path/result/PLOT/result_plot/bar_type_covid_hc.png", plot = cell_plot , width = 3.5, height = 8)



#saveRDS(merged,'/path/EQTL/result/all_subtype_annotated.rds')

#merged@meta.data$subtype[which(merged@meta.data$Type =='pDC')] <- 'pDC'

table.cluster<-table(sce@meta.data$subtype)
##计算每个样本的细胞数
table.sample<-table(sce$Status1)
##计算细胞比例
prop_allcell<-prop.table(table(sce@meta.data$subtype))
cell_num<-table(sce@meta.data$subtype, sce$Status1)
cell_numtable<-as.data.frame(table(sce@meta.data$subtype, sce$Status1))
colnames(cell_numtable)<-c("subtype","Sample","Number")
##输出cell_numtable
write.csv(cell_numtable,'/path/result/analy/CELL/cell_numtable/cell_Status1_num.csv',row.names = F)
cell_numtable

COLOR<-c("#C30E6F","#F3AD5E","#78AFD5","#EB736C","#249089","#82BF83","#EEEC96",
         "#623C89","#EC7A1F","#C1B0CA","#E0091B","#30983C","#0871B4" ,"#9CC8E0" ) 


plot_size(4,8)
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=subtype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
   scale_fill_manual(values = COLOR)+
   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size=12,color ="black"),
      axis.text.y=element_text(size=11,color ="black")) +      #让横轴上的标签倾斜45度
  xlab("Sample")  +
theme(legend.position = "right")
#scale_fill_manual(values= c("#E41A1C",'#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999'),
  # name="Experimental\nCondition",  #修改图例名称
   # labels = c("1","2","3","4","5","6","7","8"),  #修改图例标签
   # breaks = c("B","CD4 T","CD8 T","cDC","Monocyte","NK","pDC","Platelet")) #图里之前的标签
 
cell_plot

ggsave("/path/result/PLOT/result_plot/bar_subtype_status.pdf", plot = cell_plot , width = 3.5, height = 8)
ggsave("/path/result/PLOT/result_plot/bar_subtype_status.png", plot = cell_plot , width = 3.5, height = 8)

plot_size(5,8)
cell_plot<-ggplot(cell_numtable,aes(Sample,Number,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_brewer(palette = "Set1")+
   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size=12,color ="black"),
         axis.text.y=element_text(size=11,color ="black")) +      #让横轴上的标签倾斜45度
  xlab("Sample")    
cell_plot

ggsave("/path/result/PLOT/result_plot/bar_type_status.pdf", plot = cell_plot , width = 5, height = 8)
ggsave("/path/result/PLOT/result_plot/bar_type_status.png", plot = cell_plot , width = 5, height = 8)









cell_numtable<-as.data.frame(table(sce$Type, sce$orig.ident1))
colnames(cell_numtable)<-c("Celltype","Sample","Number")

cell_numtable

cell_numtable[["Status"]]<-c("NA")
cell_numtable$Status[1:40]<-("asymptomatic case")
cell_numtable$Status[41:128]<-("Healthy")
cell_numtable$Status[129:232]<-("Moderate")
cell_numtable$Status[233:328]<-("Recover")
cell_numtable$Status[329:424]<-("Severe")


celltypeall<-c()
for(i in 1:length(unique(cell_numtable$Status))){
  celltype1<-subset(cell_numtable,Status==unique(cell_numtable$Status)[i],select = c("Celltype","Sample","Number","Status"))
    for(j in 1:length(unique(celltype1$Sample))){
        celltype2<-subset(celltype1,Sample==unique(celltype1$Sample)[j],select = c("Celltype","Sample","Number","Status"))
        prop<-prop.table(celltype2$Number)
        celltype2<-cbind(celltype2,prop)
       celltypeall<-rbind(celltypeall,celltype2)
}
     #celltypel<-rbind(celltype2,celltype1)
     
    }
write.csv(celltypeall,'/path/result/analy/DEG_box/Status_cell_num_prop.csv',row.names = F)
celltypeall

col_celltypefine<-c("#FF77FF","#9955FF","#00DDDD","#00DD77","#DDFF77","#FFBB00","#EE7700","#FF0000")
#colors<- colorRampPalette(brewer.pal(n = 7, name ="Oranges"))(200)

result_path="/path/result/analy/DEG_box/plot/"          


plot_size(8,5)
celltype<-c("B" ,"CD4 T" ,"CD8 T" ,"cDC" ,"Monocyte","NK","pDC","Platelet")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
  p4<-ggplot(celltype3,aes(Status,prop,fill=Status))+  #,color=Status
    
   geom_boxplot(cex=1)+
    scale_fill_manual(values = col_celltypefine)+
    
    ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(method = "anova")+
    ggtitle(celltype[i])+
    theme(plot.title = element_text(hjust = 0.5))
 # ggsave(paste0("/path/result/analy/DEG_box/plot/",celltype[i],".pdf"),plot = p1,width = 5, height = 4)
 
}
p4

plot_size(8,5)
celltype<-c("Monocyte")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
  p4<-ggplot(celltype3,aes(Status,prop,fill=Status))+  #,color=Status
    
   geom_boxplot(cex=1)+
    scale_fill_manual(values = col_celltypefine)+
    
    ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(method = "test")+
    ggtitle(celltype[i])+
    theme(plot.title = element_text(hjust = 0.5))
 # ggsave(paste0("/path/result/analy/DEG_box/plot/",celltype[i],".pdf"),plot = p1,width = 5, height = 4)
 
}
p4

plot_size(8,5)
celltype<-c("B" ,"CD4 T" ,"CD8 T" ,"cDC" ,"NK","pDC","Platelet","Monocyte")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
  p5<-ggplot(celltype3,aes(Status,prop))+  #,color=Status
    
    geom_violin(aes(fill=Status),cex=1)+
    scale_fill_manual(values = col_celltypefine)+
    geom_boxplot(width=0.1,cex=1)+
    ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(method = "anova")+
    ggtitle(celltype[i])+
    theme(plot.title = element_text(hjust = 0.5))
  #ggsave(paste0("/path/result/analy/DEG_box/plot/",celltype[i],".pdf"),plot = p1,width = 5, height = 4)
 
}
p5

plot_size(8,5)
celltype<-c("B" ,"CD4 T" ,"CD8 T" ,"cDC" ,"NK","pDC","Platelet","Monocyte")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
  p5<-ggplot(celltype3,aes(Status,prop))+  #,color=Status
    
    geom_violin(aes(fill=Status),cex=1)+
    scale_fill_manual(values = col_celltypefine)+
    geom_boxplot(width=0.1,cex=1)+
    ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(method = "test")+
    ggtitle(celltype[i])+
    theme(plot.title = element_text(hjust = 0.5))
  #ggsave(paste0("/path/result/analy/DEG_box/plot/",celltype[i],".pdf"),plot = p1,width = 5, height = 4)
 
}
p5

plot_size(10,7)
celltype<-c("B" ,"CD4 T" ,"CD8 T" ,"cDC" ,"Monocyte","NK","pDC","Platelet")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
p<- ggplot(data=celltype3)+ 
  geom_boxplot(mapping=aes(Status,prop,fill=Status ), #箱线图
               alpha = 0.5,
               size=1,      #调整箱子外边框的宽度
               width = 0.5)+ #调整箱子的宽度
    #geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
    
  geom_jitter(mapping=aes(Status,prop,colour=Status), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("asymptomatic case","Healthy","Moderate","Recover","Severe"), 
                     values=c("#85B22E","#5F80B4","#E29827","#922927","#FF77FF"))+ #颜色
  geom_signif(mapping=aes(x=Status,y=prop), # 不同组别的显著性
              comparisons = list(  c("asymptomatic case", "Healthy"), # 哪些组进行比较
                                #  c("asymptomatic case", "Moderate"),
                                #  c("asymptomatic case", "Recover"),
                                #  c("asymptomatic case", "Severe"),
                                   c("Healthy", "Moderate"),
                                   c("Healthy", "Recover"),
                                   c("Healthy", "Severe")),
                               #   c("Moderate", "Recover"),
                               #   c("Moderate", "Severe"),
                               #   c("Recover", "Severe")),
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075), # 设置显著性线的位置高度
              size=1, # 修改差异线的粗细
              textsize = 3, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型      t.test() 比较两组(参数)
                                                 #wilcox.test() 比较两组(非参数)
                                                 #aov()或anova() 比较多组(参数)
                                                 #kruskal.test() 比较多组(非参数)
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  labs(title=celltype[i],x="Status",y="Prop")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_text(size = 15, 
                                    #family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5, 
                                    angle = 0),
          
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )+
theme(axis.text.x = element_text(angle = 45, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))
    
 print(p)  
    }


#celltypeal<-c()
#for(i in 1:length(unique(cell_numtable$Sample))){
  celltype1<-subset(cell_numtable,Sample==unique(cell_numtable$Sample)[i],select = c("Celltype","Sample","Number"))
  prop<-prop.table(celltype1$Number)
  celltype1<-cbind(celltype1,prop)
  celltypeall<-rbind(celltypeall,celltype1)
}
write.csv(celltypeall,'/path/result/analy/DEG_box/Status_cell_num_prop.csv',row.names = F)
celltypeall

options(repr.plot.width = 18, repr.plot.height = 8)
ggplot(celltypeall,aes(x=Status,y=prop))+
geom_boxplot(aes(fill=Status))+
facet_grid(.~ Celltype)+
 ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(method = "anova")+
    ggtitle(celltype)+
    theme(plot.title = element_text(hjust = 0.5))+
theme(axis.text.x = element_text(angle = 45, 
                                 vjust = 0.5, hjust=0.5,size=10),axis.text.y = element_text(size=12))


    

   
   
 

options(repr.plot.width = 18, repr.plot.height = 8)
ggplot(celltypeall,aes(x=Celltype,y=prop))+
geom_boxplot(aes(fill=Status))+
facet_grid(.~ Status)+
 ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(method = "anova") 
ggtitle(celltype)+
    theme(plot.title = element_text(hjust = 0.5))+
theme(axis.text.x = element_text(angle = 45, 
                                 vjust = 0.5, hjust=0.5,size=10),axis.text.y = element_text(size=12))

celltypeall
write.table(celltypeall,"celltypeall.txt",quote=F,sep="\t")
getwd()

#运行成对T-test
stat.test <- celltypeall%>%
  group_by(Status) %>%
  pairwise_t_test(
    prop~ Celltype, paired = TRUE, 
    p.adjust.method = "bonferroni"
    ) %>%
  select(-df, -statistic, -p) # Remove details
stat.test



options(repr.plot.width = 18, repr.plot.height = 8)
ggplot(celltypeall,aes(x=Celltype,y=prop))+
geom_boxplot(aes(fill=Status))+
facet_grid(.~ Status)+
 ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(method = "anova") 
ggtitle(celltype)+
    theme(plot.title = element_text(hjust = 0.5))+
theme(axis.text.x = element_text(angle = 45, 
                                 vjust = 0.5, hjust=0.5,size=10),axis.text.y = element_text(size=12))

stat.test <- stat.test %>% add_xy_position(x = "Status")
stat.test 

bxp <- ggboxplot(
  celltypeall, x = "Status", y = "prop",
  color = "Celltype", palette = "jco"
  )
bxp

## 添加显著性标记
bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.05
  )

#隐藏不显著项
bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE, tip.length = 0
  )



plot_size(10,7)
celltype<-c("pDC")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
p<- ggplot(data=celltype3)+ 
  geom_boxplot(mapping=aes(Status,prop,fill=Status ), #箱线图
               alpha = 0.5,
               size=1,      #调整箱子外边框的宽度
               width = 0.5)+ #调整箱子的欢度
    #geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
    
  geom_jitter(mapping=aes(Status,prop,colour=Status), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("asymptomatic case","Healthy","Moderate","Recover","Severe"), 
                     values=c("#85B22E","#5F80B4","#E29827","#922927","#FF77FF"))+ #颜色
  geom_signif(mapping=aes(x=Status,y=prop), # 不同组别的显著性
              comparisons = list( # c("asymptomatic case", "Healthy"), # 哪些组进行比较
                                 c("asymptomatic case", "Moderate"),
                                  c("asymptomatic case", "Recover"),
                                  c("asymptomatic case", "Severe"),
                                   c("Healthy", "Moderate"),
                                   c("Healthy", "Recover"),
                                   c("Healthy", "Severe")),
                                  #c("Moderate", "Recover"),
                                 # c("Moderate", "Severe"),
                                 # c("Recover", "Severe")),
              map_signif_level=T, #T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.1), # 设置显著性线的位置高度
              size=1, # 修改差异线的粗细
              textsize = 3, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型      t.test() 比较两组(参数)
                                                 #wilcox.test() 比较两组(非参数)
                                                 #aov()或anova() 比较多组(参数)
                                                 #kruskal.test() 比较多组(非参数)
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  labs(title=celltype[i],x="Status",y="Prop")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_text(size = 15, 
                                    #family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5, 
                                    angle = 0),
          
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )+
theme(axis.text.x = element_text(angle = 45, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))
    
    
    }
p


plot_size(10,7)
celltype<-c("Platelet")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
p<- ggplot(data=celltype3)+ 
  geom_boxplot(mapping=aes(Status,prop,fill=Status ), #箱线图
               alpha = 0.5,
               size=1,      #调整箱子外边框的宽度
               width = 0.5)+ #调整箱子的欢度
    #geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
    
  geom_jitter(mapping=aes(Status,prop,colour=Status), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("asymptomatic case","Healthy","Moderate","Recover","Severe"), 
                     values=c("#85B22E","#5F80B4","#E29827","#922927","#FF77FF"))+ #颜色
  geom_signif(mapping=aes(x=Status,y=prop), # 不同组别的显著性
              comparisons = list(  c("asymptomatic case", "Healthy"), # 哪些组进行比较
                                # c("asymptomatic case", "Moderate"),
                                  #c("asymptomatic case", "Recover"),
                                 # c("asymptomatic case", "Severe"),
                                   #c("Healthy", "Moderate"),
                                   #c("Healthy", "Recover"),
                                   c("Healthy", "Severe"),
                                 # c("Moderate", "Recover"),
                                  c("Moderate", "Severe"),
                                  c("Recover", "Severe")),
              map_signif_level=T, #T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.1,0.15), # 设置显著性线的位置高度
              size=1, # 修改差异线的粗细
              textsize = 3, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型      t.test() 比较两组(参数)
                                                 #wilcox.test() 比较两组(非参数)
                                                 #aov()或anova() 比较多组(参数)
                                                 #kruskal.test() 比较多组(非参数)
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  labs(title=celltype[i],x="Status",y="Prop")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_text(size = 15, 
                                    #family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5, 
                                    angle = 0),
          
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )+
theme(axis.text.x = element_text(angle = 45, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))
    
    
    }
p


plot_size(10,7)
celltype<-c("CD8 T")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
p<- ggplot(data=celltype3)+ 
  geom_boxplot(mapping=aes(Status,prop,fill=Status ), #箱线图
               alpha = 0.5,
               size=1,      #调整箱子外边框的宽度
               width = 0.5)+ #调整箱子的欢度
    #geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
    
  geom_jitter(mapping=aes(Status,prop,colour=Status), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("asymptomatic case","Healthy","Moderate","Recover","Severe"), 
                     values=c("#85B22E","#5F80B4","#E29827","#922927","#FF77FF"))+ #颜色
  geom_signif(mapping=aes(x=Status,y=prop), # 不同组别的显著性
              comparisons = list(  c("asymptomatic case", "Healthy"), # 哪些组进行比较
                                 #c("asymptomatic case", "Moderate"),
                                  #c("asymptomatic case", "Recover"),
                                  #c("asymptomatic case", "Severe"),
                                  # c("Healthy", "Moderate"),
                                  # c("Healthy", "Recover"),
                                   c("Healthy", "Severe")),
                                  #c("Moderate", "Recover"),
                                 # c("Moderate", "Severe"),
                                  #c("Recover", "Severe")),
              map_signif_level=T, #T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(0.4,0.43,0.4,0.45,0.5,0.55,0.65,0.7,0.75,0.8), # 设置显著性线的位置高度
              size=1, # 修改差异线的粗细
              textsize = 3, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型      t.test() 比较两组(参数)
                                                 #wilcox.test() 比较两组(非参数)
                                                 #aov()或anova() 比较多组(参数)
                                                 #kruskal.test() 比较多组(非参数)
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  labs(title=celltype[i],x="Status",y="Prop")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_text(size = 15, 
                                    #family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5, 
                                    angle = 0),
          
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )+
theme(axis.text.x = element_text(angle = 45, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))
    
    
    }
p




plot_size(10,7)
celltype<-c("CD4 T")
for (i in 1:length(celltype)) {
  celltype3<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
p<- ggplot(data=celltype3)+ 
  geom_boxplot(mapping=aes(Status,prop,fill=Status ), #箱线图
               alpha = 0.5,
               size=1,      #调整箱子外边框的宽度
               width = 0.5)+ #调整箱子的欢度
    #geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
    
  geom_jitter(mapping=aes(Status,prop,colour=Status), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("asymptomatic case","Healthy","Moderate","Recover","Severe"), 
                     values=c("#85B22E","#5F80B4","#E29827","#922927","#FF77FF"))+ #颜色
  geom_signif(mapping=aes(x=Status,y=prop), # 不同组别的显著性
              comparisons = list(  # c("asymptomatic case", "Healthy"), # 哪些组进行比较
                                #  c("asymptomatic case", "Moderate"),
                                  # c("asymptomatic case", "Recover"),
                                  # c("asymptomatic case", "Severe"),
                                   # c("Healthy", "Moderate"),
                                  # c("Healthy", "Recover"),
                                   c("Healthy", "Severe"),
                                  #c("Moderate", "Recover"),
                                  c("Moderate", "Severe")),
                                 # c("Recover", "Severe")),
              map_signif_level=T, #T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(0.5,0.55,0.4,0.45,0.5,0.55,0.65,0.7,0.75,0.8), # 设置显著性线的位置高度
              size=1, # 修改差异线的粗细
              textsize = 3, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型      t.test() 比较两组(参数)
                                                 #wilcox.test() 比较两组(非参数)
                                                 #aov()或anova() 比较多组(参数)
                                                 #kruskal.test() 比较多组(非参数)
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  labs(title=celltype[i],x="Status",y="Prop")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_text(size = 15, 
                                    #family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5, 
                                    angle = 0),
          
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )+
theme(axis.text.x = element_text(angle = 45, 
                                 vjust = 0.5, hjust=0.5,size=15),axis.text.y = element_text(size=15))
    
    
    }
p


save.image("/path/RData/5080_plot_result1.RData")


