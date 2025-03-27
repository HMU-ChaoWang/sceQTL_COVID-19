library(magrittr)
library(data.table)

dim(gwas)

gwas<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.txt")

eqtl<- read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_cis/B_cis_eqtls.txt")
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
data.frame(snpexp)

B_exp_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_exp/B.txt")
B_exp_anno<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")


head(gwas)
head(eqtl)
head(snp)
head(genotype_case)
head(B_exp_case)
head(B_exp_anno)

####整理snp表型为TCGA#### 
for(i in 1:nrow(snpexp)){
  for(j in 1:ncol(snpexp)){
    if(is.na(snpexp[i,j])){
      snpexp[i,j]<-paste(0,0,sep=" ")
    }else if(snpexp[i,j]==0){
      snpexp[i,j]<-paste(snp[i,"ref"],snp[i,"ref"],sep = " ")
    }else if(snpexp[i,j]==1){
      snpexp[i,j]<-paste(snp[i,"ref"],snp[i,"alt"],sep = " ")
    }else if(snpexp[i,j]==2){
      snpexp[i,j]<-paste(snp[i,"alt"],snp[i,"alt"],sep = " ")
    }
      
    }
  }

ch<-subset(snp,!snp$alt %in% c("T","C","G","A"))
for(i in rownames(ch)){
  ref<-snp[i,"ref"]
  alt<-substr(snp[i,"alt"],1,1)
  for(j in 1:ncol(snpexp)){
    if(is.na(genotype_case[i,j])){
      snpexp[i,j]<-paste(0,0,sep=" ")
    }else if(genotype_case[i,j]==0){
      snpexp[i,j]<-paste(ref,ref,sep = " ")
    }else if(genotype_case[i,j]==1){
      snpexp[i,j]<-paste(ref,alt,sep = " ")
    }else if(genotype_case[i,j]==2){
      snpexp[i,j]<-paste(alt,alt,sep = " ")
    }
  }
  
}

gene<-as.vector(unique(eqtl$gene)) #获取所有基因名
length(gene)

gene<-as.vector(unique(eqtl$gene))
gene[1:2000]

#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[1:2000])
{
   mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/1_2000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[1:2000])
{
     mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
     exp1<-B_exp_case[which( B_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/1_2000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(eqtl$gene))[1:2000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/genename1_2000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(eqtl$gene))
gene[2001:4000]

#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[2001:4000])
{
   mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/2001_4000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[2001:4000])
{
     mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
     exp1<-B_exp_case[which( B_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/2001_4000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(eqtl$gene))[2001:4000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/genename2001_4000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(eqtl$gene))
gene[4001:6000]

#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()
gene<-as.vector(unique(eqtl$gene))
for(i in gene[4001:6000])
{
   mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/4001_6000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[4001:6000])
{
     mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
     exp1<-B_exp_case[which( B_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/4001_6000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(eqtl$gene))[4001:6000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/genename4001_6000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(eqtl$gene))
gene[6001:8000]

#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[6001:8000])
{
   mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/6001_8000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[6001:8000])
{
     mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
     exp1<-B_exp_case[which( B_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/6001_8000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(eqtl$gene))[6001:8000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/genename6001_8000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(eqtl$gene))
gene[8001:10000]

#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[8001:10000])
{
   mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/8001_10000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[8001:10000])
{
     mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
     exp1<-B_exp_case[which( B_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/8001_10000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(eqtl$gene))[8001:10000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/genename8001_10000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(eqtl$gene))
gene[10001:11008]
length(gene)

#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[10001:11008])
{
   mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/10001_11008/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[10001:11008])
{
     mapsnp<-subset(eqtl,eqtl$gene==i)[,1]
     exp1<-B_exp_case[which( B_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/10001_11008/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}


gene<-as.vector(unique(eqtl$gene))[10001:11008]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/genename10001_11008.txt",quote=F,row.names=F,col.names=F,sep="\t")


#while循环每个基因进行plink
#!!注意，如果是用windows传过来的文件，一定要去除回车符\r，否则没法在linux上while

sed -i 's/\r//g' /data_alluser/CXY/CXY/zsn_coloc/genename10001_11663.txt
cat /data_alluser/CXY/CXY/zsn_coloc/genename10001_11663.txt | while read line
do
/data_alluser/CXY/CXY/app/Plink/plink \
--file /data_alluser/CXY/CXY/zsn_coloc/map/10001_11663/${line} \
--make-bed \
--out /data_alluser/CXY/CXY/zsn_coloc/bbf/10001_11663/${line} \
--noweb --no-parents
done

ln -s ./ output 
cat /path/singlecell_eqtl/result/TWAS_DATA/celltype/B/map/genename8001_10000.txt | while read line
do
Rscript /data_alluser/CXY/CXY/TWAS/fusion_twas-master/FUSION.compute_weights.R \
--bfile /path/singlecell_eqtl/result/TWAS_DATA/celltype/B/bbf/8001_10000/${line} \
--tmp ./output/ \
--out /path/singlecell_eqtl/result/TWAS_DATA/celltype/B/out/1/${line} \
--models top1,lasso,enet \
--crossval 5 \
--PATH_plink /data_alluser/CXY/CXY/app/plink_linux_x86_64_20210606/plink \
--PATH_gcta /data_alluser/CXY/CXY/TWAS/fusion_twas-master/gcta_nr_robust \
--PATH_gemma /data_alluser/CXY/CXY/TWAS/QZ/gemma-0.98.5-linux-static-AMD64/gemma-0.98.5 \
--hsq_p 0.05 \
--noclean
done

file<-list.files("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/out/weight",pattern="*.RDat")
file
length(file)

a<-unique(file)
length(a)
a

write.table(a,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/out/B_filename_weight_3.txt",col.names=F,row.names=F,quote=F)

filename_weight_gene<-fread("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/out/pro_pos_B_filename_weight_3.txt")
filename_weight_gene

gene_anno<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")
gene_anno

merge<-merge(filename_weight_gene,gene_anno,by.x="ID",by.y="id",all=F)
merge

merge<-merge[,c(2,3,1,4,5,6)]
colnames(merge)<-c("PANEL","WGT","ID","CHR","P0","P1")
merge

write.table(merge,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/B/out/B_merge_weight_3.txt",row.names=F,quote=F,sep="\t")


sed -i 's/\r//g' /path/singlecell_eqtl/result/TWAS_DATA/celltype/B/out/B_merge_weight_3.pos   

sed -i 's/\r//g' /path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.sumstats


 #! /bin/bash
for((i=15;i<=20;i++))
 do
Rscript /data_alluser/CXY/CXY/TWAS/fusion_twas-master/FUSION.assoc_test.R \
--sumstats /path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.sumstats \
--weights /path/singlecell_eqtl/result/TWAS_DATA/celltype/B/out/B_merge_weight_3.pos \
--weights_dir /path/singlecell_eqtl/result/TWAS_DATA/celltype/B/out/weight/ \
--ref_ld_chr /data_alluser/CXY/CXY/LD_SNP/1000G.EUR. \
--chr ${i} \
--out /path/singlecell_eqtl/result/TWAS_DATA/celltype/B/result/B.${i}.dat
 done

setwd("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B1/out/result")
TWAS<-list.files("/path/singlecell_eqtl/result/TWAS_DATA/celltype/B1/out/result",pattern="*.dat")
TWAS



merge<-list()
mer<-list()
for(i in 1:length(TWAS)){

  twas<-read.csv(file=TWAS[i],header=TRUE,sep="\t")
  me<-rbind(merge,twas)
  mer<-rbind(mer,me)
  print(c(i,dim(twas)[1],dim(mer)[1]))
 }


TWAS_Bonferroni <- p.adjust(as.numeric(mer$TWAS.P),method = "bonferroni")
jiaozheng<-cbind(mer,TWAS_Bonferroni)

TWAS_ner<-jiaozheng[which(as.numeric(jiaozheng$TWAS_Bonferroni)<0.05),]

dim(TWAS_ner)

write.table(TWAS_ner,"B_eqtl_twas.txt",sep="\t",quote=F,row.names=F,col.names=T)

write.table(jiaozheng,"B_TAWS_jiaozheng_P.txt",quote=F,row.names=F,col.names=T,sep="\t")

gwas<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.txt")
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
exp_anno<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
data.frame(snpexp)

cd4_eqtl<- read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_cis/CD4_T_cis_eqtls.txt")
cd4_exp_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_exp/CD4_T.txt")

head(gwas)
head(snp)
head(genotype_case)
head(exp_anno)

head(cd4_eqtl)
head(cd4_exp_case)

#上面跑过就不跑了
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
####整理snp表型为TCGA#### 
for(i in 1:nrow(snpexp)){
  for(j in 1:ncol(snpexp)){
    if(is.na(snpexp[i,j])){
      snpexp[i,j]<-paste(0,0,sep=" ")
    }else if(snpexp[i,j]==0){
      snpexp[i,j]<-paste(snp[i,"ref"],snp[i,"ref"],sep = " ")
    }else if(snpexp[i,j]==1){
      snpexp[i,j]<-paste(snp[i,"ref"],snp[i,"alt"],sep = " ")
    }else if(snpexp[i,j]==2){
      snpexp[i,j]<-paste(snp[i,"alt"],snp[i,"alt"],sep = " ")
    }
      
    }
  }

ch<-subset(snp,!snp$alt %in% c("T","C","G","A"))
for(i in rownames(ch)){
  ref<-snp[i,"ref"]
  alt<-substr(snp[i,"alt"],1,1)
  for(j in 1:ncol(snpexp)){
    if(is.na(genotype_case[i,j])){
      snpexp[i,j]<-paste(0,0,sep=" ")
    }else if(genotype_case[i,j]==0){
      snpexp[i,j]<-paste(ref,ref,sep = " ")
    }else if(genotype_case[i,j]==1){
      snpexp[i,j]<-paste(ref,alt,sep = " ")
    }else if(genotype_case[i,j]==2){
      snpexp[i,j]<-paste(alt,alt,sep = " ")
    }
  }
  
}

gene<-as.vector(unique(cd4_eqtl$gene)) #获取所有基因名
length(gene)

gene<-as.vector(unique(cd4_eqtl$gene))
gene[1:2000]

#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[1:2000])
{
   mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/1_2000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[1:2000])
{
     mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
     exp1<-cd4_exp_case[which( cd4_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/1_2000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(cd4_eqtl$gene))[1:2000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/genename1_2000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd4_eqtl$gene))
gene[2001:4000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[2001:4000])
{
   mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/2001_4000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表  4
pedlist<-list()
for(i in gene[2001:4000])
{
     mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
     exp1<-cd4_exp_case[which( cd4_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/2001_4000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(cd4_eqtl$gene))[2001:4000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/genename2001_4000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd4_eqtl$gene))
gene[4001:6000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[4001:6000])
{
   mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/4001_6000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表  4
pedlist<-list()
for(i in gene[4001:6000])
{
     mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
     exp1<-cd4_exp_case[which( cd4_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/4001_6000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(cd4_eqtl$gene))[4001:6000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/genename4001_6000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd4_eqtl$gene))
gene[6001:8000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[6001:8000])
{
   mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/6001_8000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表  4
pedlist<-list()
for(i in gene[6001:8000])
{
     mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
     exp1<-cd4_exp_case[which( cd4_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/6001_8000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(cd4_eqtl$gene))[6001:8000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/genename6001_8000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd4_eqtl$gene))
gene[8001:10000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[8001:10000])
{
   mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/8001_10000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表  4
pedlist<-list()
for(i in gene[8001:10000])
{
     mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
     exp1<-cd4_exp_case[which( cd4_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/8001_10000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(cd4_eqtl$gene))[8001:10000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/genename8001_10000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd4_eqtl$gene))
gene[10001:11304]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[10001:11304])
{
   mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/10001_11304/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表  4
pedlist<-list()
for(i in gene[10001:11304])
{
     mapsnp<-subset(cd4_eqtl,cd4_eqtl$gene==i)[,1]
     exp1<-cd4_exp_case[which( cd4_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/10001_11304/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(cd4_eqtl$gene))[10001:11304]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD4/map/genename10001_11304.txt",quote=F,row.names=F,col.names=F,sep="\t")




gwas<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.txt")
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
exp_anno<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
data.frame(snpexp)

cd8_eqtl<- read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_cis/CD8_T_cis_eqtls.txt")
cd8_exp_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_exp/CD8_T.txt")

head(gwas)
head(snp)
head(genotype_case)
head(exp_anno)

head(cd8_eqtl)
head(cd8_exp_case)

gene<-as.vector(unique(cd8_eqtl$gene)) #获取所有基因名
length(gene)

gene<-as.vector(unique(cd8_eqtl$gene))
gene[1:2000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[1:2000])
{
   mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/1_2000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[1:2000])
{
     mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
     exp1<-cd8_exp_case[which( cd8_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/1_2000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(cd8_eqtl$gene))[1:2000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/genename1_2000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd8_eqtl$gene))
gene[2001:4000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[2001:4000])
{
   mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/2001_4000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[2001:4000])
{
     mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
     exp1<-cd8_exp_case[which( cd8_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/2001_4000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}


gene<-as.vector(unique(cd8_eqtl$gene))[2001:4000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/genename2001_4000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd8_eqtl$gene))
gene[4001:6000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[4001:6000])
{
   mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/4001_6000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[4001:6000])
{
     mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
     exp1<-cd8_exp_case[which( cd8_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/4001_6000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}


gene<-as.vector(unique(cd8_eqtl$gene))[4001:6000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/genename4001_6000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd8_eqtl$gene))
gene[6001:8000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[6001:8000])
{
   mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/6001_8000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[6001:8000])
{
     mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
     exp1<-cd8_exp_case[which( cd8_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/6001_8000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}


gene<-as.vector(unique(cd8_eqtl$gene))[6001:8000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/genename6001_8000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd8_eqtl$gene))
gene[8001:10000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[8001:10000])
{
   mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/8001_10000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[8001:10000])
{
     mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
     exp1<-cd8_exp_case[which( cd8_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/8001_10000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}


gene<-as.vector(unique(cd8_eqtl$gene))[8001:10000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/genename8001_10000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cd8_eqtl$gene))
gene[10001:11228]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[10001:11228])
{
   mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/10001_11228/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[10001:11228])
{
     mapsnp<-subset(cd8_eqtl,cd8_eqtl$gene==i)[,1]
     exp1<-cd8_exp_case[which( cd8_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/10001_11228/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}


gene<-as.vector(unique(cd8_eqtl$gene))[10001:11228]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/CD8/map/genename10001_11228.txt",quote=F,row.names=F,col.names=F,sep="\t")








gwas<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.txt")
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
exp_anno<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
data.frame(snpexp)

cDC_eqtl<- read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_cis/cDC_cis_eqtls.txt")
cDC_exp_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_exp/cDC.txt")

head(gwas)
head(snp)
head(genotype_case)
head(exp_anno)

head(cDC_eqtl)
head(cDC_exp_case)

gene<-as.vector(unique(cDC_eqtl$gene)) #获取所有基因名
length(gene)

gene<-as.vector(unique(cDC_eqtl$gene))
gene[1:2000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[1:2000])
{
   mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/1_2000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[1:2000])
{
     mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
     exp1<-cDC_exp_case[which( cDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/1_2000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(cDC_eqtl$gene))[1:2000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/genename1_2000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cDC_eqtl$gene))
gene[2001:4000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[2001:4000])
{
   mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/2001_4000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[2001:4000])
{
     mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
     exp1<-cDC_exp_case[which( cDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/2001_4000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}
gene<-as.vector(unique(cDC_eqtl$gene))[2001:4000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/genename2001_4000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cDC_eqtl$gene))
gene[4001:6000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[4001:6000])
{
   mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/4001_6000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[4001:6000])
{
     mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
     exp1<-cDC_exp_case[which( cDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/4001_6000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}
gene<-as.vector(unique(cDC_eqtl$gene))[4001:6000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/genename4001_6000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cDC_eqtl$gene))
gene[6001:8000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[6001:8000])
{
   mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/6001_8000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[6001:8000])
{
     mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
     exp1<-cDC_exp_case[which( cDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/6001_8000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}
gene<-as.vector(unique(cDC_eqtl$gene))[6001:8000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/genename6001_8000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cDC_eqtl$gene))
gene[8001:10000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[8001:10000])
{
   mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/8001_10000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[8001:10000])
{
     mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
     exp1<-cDC_exp_case[which( cDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/8001_10000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}
gene<-as.vector(unique(cDC_eqtl$gene))[8001:10000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/genename8001_10000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(cDC_eqtl$gene))
gene[10001:10235]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[10001:10235])
{
   mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/10001_10235/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[10001:10235])
{
     mapsnp<-subset(cDC_eqtl,cDC_eqtl$gene==i)[,1]
     exp1<-cDC_exp_case[which( cDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/10001_10235/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}
gene<-as.vector(unique(cDC_eqtl$gene))[10001:10235]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/cDC/map/genename10001_10235.txt",quote=F,row.names=F,col.names=F,sep="\t")




gwas<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.txt")
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
exp_anno<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
data.frame(snpexp)

Monocyte_eqtl<- read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_cis/Monocyte_cis_eqtls.txt")
Monocyte_exp_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_exp/Monocyte.txt")

head(gwas)
head(snp)
head(genotype_case)
head(exp_anno)

head(Monocyte_eqtl)
head(Monocyte_exp_case)

gene<-as.vector(unique(Monocyte_eqtl$gene)) #获取所有基因名
length(gene)

gene<-as.vector(unique(Monocyte_eqtl$gene))
gene[1:2000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[1:2000])
{
   mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/1_2000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[1:2000])
{
     mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
     exp1<-Monocyte_exp_case[which( Monocyte_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/1_2000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Monocyte_eqtl$gene))[1:2000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/genename1_2000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Monocyte_eqtl$gene))
gene[2001:4000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[2001:4000])
{
   mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/2001_4000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[2001:4000])
{
     mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
     exp1<-Monocyte_exp_case[which( Monocyte_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/2001_4000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Monocyte_eqtl$gene))[2001:4000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/genename2001_4000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Monocyte_eqtl$gene))
gene[4001:6000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[4001:6000])
{
   mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/4001_6000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[4001:6000])
{
     mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
     exp1<-Monocyte_exp_case[which( Monocyte_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/4001_6000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Monocyte_eqtl$gene))[4001:6000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/genename4001_6000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Monocyte_eqtl$gene))
gene[6001:8000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[6001:8000])
{
   mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/6001_8000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[6001:8000])
{
     mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
     exp1<-Monocyte_exp_case[which( Monocyte_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/6001_8000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Monocyte_eqtl$gene))[6001:8000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/genename6001_8000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Monocyte_eqtl$gene))
gene[8001:10000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[8001:10000])
{
   mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/8001_10000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[8001:10000])
{
     mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
     exp1<-Monocyte_exp_case[which( Monocyte_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/8001_10000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Monocyte_eqtl$gene))[8001:10000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/genename8001_10000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Monocyte_eqtl$gene))
gene[10001:11542]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[10001:11542])
{
   mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/10001_11542/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[10001:11542])
{
     mapsnp<-subset(Monocyte_eqtl,Monocyte_eqtl$gene==i)[,1]
     exp1<-Monocyte_exp_case[which( Monocyte_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/10001_11542/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Monocyte_eqtl$gene))[10001:11542]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Monocyte/map/genename10001_11542.txt",quote=F,row.names=F,col.names=F,sep="\t")




gwas<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.txt")
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
exp_anno<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
data.frame(snpexp)

NK_eqtl<- read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_cis/NK_cis_eqtls.txt")
NK_exp_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_exp/NK.txt")

head(gwas)
head(snp)
head(genotype_case)
head(exp_anno)

head(NK_eqtl)
head(NK_exp_case)

gene<-as.vector(unique(NK_eqtl$gene)) #获取所有基因名
length(gene)

gene<-as.vector(unique(NK_eqtl$gene))
gene[1:2000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[1:2000])
{
   mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/1_2000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[1:2000])
{
     mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
     exp1<-NK_exp_case[which( NK_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/1_2000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(NK_eqtl$gene))[1:2000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/genename1_2000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(NK_eqtl$gene))
gene[2001:4000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[2001:4000])
{
   mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/2001_4000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[2001:4000])
{
     mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
     exp1<-NK_exp_case[which( NK_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/2001_4000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(NK_eqtl$gene))[2001:4000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/genename2001_4000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(NK_eqtl$gene))
gene[4001:6000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[4001:6000])
{
   mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/4001_6000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[4001:6000])
{
     mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
     exp1<-NK_exp_case[which( NK_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/4001_6000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(NK_eqtl$gene))[4001:6000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/genename4001_6000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(NK_eqtl$gene))
gene[6001:8000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[6001:8000])
{
   mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/6001_8000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[6001:8000])
{
     mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
     exp1<-NK_exp_case[which( NK_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/6001_8000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(NK_eqtl$gene))[6001:8000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/genename6001_8000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(NK_eqtl$gene))
gene[8001:9991]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[8001:9991])
{
   mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/8001_9991/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[8001:9991])
{
     mapsnp<-subset(NK_eqtl,NK_eqtl$gene==i)[,1]
     exp1<-NK_exp_case[which( NK_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/8001_9991/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(NK_eqtl$gene))[8001:9991]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/NK/map/genename8001_9991.txt",quote=F,row.names=F,col.names=F,sep="\t")




gwas<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.txt")
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
exp_anno<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
data.frame(snpexp)

pDC_eqtl<- read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_cis/pDC_cis_eqtls.txt")
pDC_exp_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_exp/pDC.txt")

head(gwas)
head(snp)
head(genotype_case)
head(exp_anno)

head(pDC_eqtl)
head(pDC_exp_case)

gene<-as.vector(unique(pDC_eqtl$gene)) #获取所有基因名
length(gene)

gene<-as.vector(unique(pDC_eqtl$gene))
gene[1:2000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[1:2000])
{
   mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/1_2000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[1:2000])
{
     mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
     exp1<-pDC_exp_case[which( pDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/1_2000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(pDC_eqtl$gene))[1:2000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/genename1_2000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(pDC_eqtl$gene))
gene[2001:4000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[2001:4000])
{
   mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/2001_4000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[2001:4000])
{
     mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
     exp1<-pDC_exp_case[which( pDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/2001_4000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}


gene<-as.vector(unique(pDC_eqtl$gene))[2001:4000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/genename2001_4000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(pDC_eqtl$gene))
gene[4001:6000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[4001:6000])
{
   mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/4001_6000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[4001:6000])
{
     mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
     exp1<-pDC_exp_case[which( pDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/4001_6000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(pDC_eqtl$gene))[4001:6000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/genename4001_6000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(pDC_eqtl$gene))
gene[6001:8000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[6001:8000])
{
   mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/6001_8000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[6001:8000])
{
     mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
     exp1<-pDC_exp_case[which( pDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/6001_8000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(pDC_eqtl$gene))[6001:8000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/genename6001_8000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(pDC_eqtl$gene))
gene[8001:10000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[8001:10000])
{
   mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/8001_10000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[8001:10000])
{
     mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
     exp1<-pDC_exp_case[which( pDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/8001_10000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(pDC_eqtl$gene))[8001:10000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/genename8001_10000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(pDC_eqtl$gene))
gene[10001:10155]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[10001:10155])
{
   mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/10001_10155/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[10001:10155])
{
     mapsnp<-subset(pDC_eqtl,pDC_eqtl$gene==i)[,1]
     exp1<-pDC_exp_case[which( pDC_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/10001_10155/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(pDC_eqtl$gene))[10001:10155]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/pDC/map/genename10001_10155.txt",quote=F,row.names=F,col.names=F,sep="\t")




gwas<-fread("/path/singlecell_eqtl/result/TWAS_DATA/need/ALL_gwas_TWAS.txt")
snp<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/snp_anno.txt")
exp_anno<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/gene_anno_unique_merge.txt")
genotype_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/TWAS_genotype_transpose.txt")
snpexp<-genotype_case
data.frame(snpexp)

Platelet_eqtl<- read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_cis/Platelet_cis_eqtls.txt")
Platelet_exp_case<-read.delim("/path/singlecell_eqtl/result/TWAS_DATA/need/type_exp/Platelet.txt")

head(gwas)
head(snp)
head(genotype_case)
head(exp_anno)

head(Platelet_eqtl)
head(Platelet_exp_case)

gene<-as.vector(unique(Platelet_eqtl$gene)) #获取所有基因名
length(gene)

gene<-as.vector(unique(Platelet_eqtl$gene))
gene[1:2000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[1:2000])
{
   mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/1_2000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[1:2000])
{
     mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
     exp1<-Platelet_exp_case[which( Platelet_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/1_2000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Platelet_eqtl$gene))[1:2000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/genename1_2000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Platelet_eqtl$gene))
gene[2001:4000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[2001:4000])
{
   mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/2001_4000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[2001:4000])
{
     mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
     exp1<-Platelet_exp_case[which( Platelet_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/2001_4000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}


gene<-as.vector(unique(Platelet_eqtl$gene))[2001:4000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/genename2001_4000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Platelet_eqtl$gene))
gene[4001:6000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[4001:6000])
{
   mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/4001_6000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[4001:6000])
{
     mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
     exp1<-Platelet_exp_case[which( Platelet_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/4001_6000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Platelet_eqtl$gene))[4001:6000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/genename4001_6000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Platelet_eqtl$gene))
gene[6001:8000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[6001:8000])
{
   mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/6001_8000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[6001:8000])
{
     mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
     exp1<-Platelet_exp_case[which( Platelet_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/6001_8000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Platelet_eqtl$gene))[6001:8000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/genename6001_8000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Platelet_eqtl$gene))
gene[8001:10000]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[8001:10000])
{
   mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/8001_10000/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[8001:10000])
{
     mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
     exp1<-Platelet_exp_case[which( Platelet_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/8001_10000/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Platelet_eqtl$gene))[8001:10000]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/genename8001_10000.txt",quote=F,row.names=F,col.names=F,sep="\t")


gene<-as.vector(unique(Platelet_eqtl$gene))
gene[10001:10258]
#一个文件夹不要放太多文件)
#制作map列表
maplist<-list()
maplist1<-list()
maplist2<-list()
maplist3<-list()

for(i in gene[10001:10258])
{
   mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
 
     for( j in 1:length(mapsnp ))
  {
       snp1<-snp[which( snp[,1]==mapsnp[j] ),]
       mapdf<-data.frame(chr=snp1[,2],snp=snp1[,1],x=0,pos=snp1[,3])
       maplist<-rbind(maplist,i)
       maplist1<-rbind(maplist1,mapdf)
       maplist2<-cbind(maplist1,maplist)
   }
    
   maplist3[[i]]<-maplist2[which(maplist2[,5]==i),-5]
  write.table(maplist3[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/10001_10258/",i,".map",sep=""),quote=F,row.names=F,col.names=F,sep=" ")

}

#制作ped列表
pedlist<-list()
for(i in gene[10001:10258])
{
     mapsnp<-subset(Platelet_eqtl,Platelet_eqtl$gene==i)[,1]
     exp1<-Platelet_exp_case[which( Platelet_exp_case[,1]==i ),]
     peddf<-data.frame(FamilyID=colnames(exp1[,2:dim(exp1)[2]]) ,PersonalID=colnames(exp1[,2:dim(exp1)[2]]),
                       Sex=0,genexp=as.vector(t(exp1[,2:dim(exp1)[2]])))
        for( j in 1:length(mapsnp ))
           {
                 snpexp1<-snpexp[which( snpexp[,1]==mapsnp[j] ),]
                 peddf<-cbind(peddf,t(snpexp1[,2:dim(snpexp1)[2]]))
                 

            }
      pedlist[[i]]<-peddf
     write.table(pedlist[[i]],file=paste("/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/10001_10258/",i,".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
}

gene<-as.vector(unique(Platelet_eqtl$gene))[10001:10258]
write.table(gene,"/path/singlecell_eqtl/result/TWAS_DATA/celltype/Platelet/map/genename10001_10258.txt",quote=F,row.names=F,col.names=F,sep="\t")





