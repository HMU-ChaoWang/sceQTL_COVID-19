library(data.table)
library(MatrixEQTL)

file_path<-"/path1/"
expfiles_dir<-"/path2/"
result_path<-"/path3/"

allfiles<-list.files(expfiles_dir)


# 导入数据
SNP_file_name = (paste0(file_path,"genotype.txt")) # 获取SNP文件位置
snps_location_file_name = (paste0(file_path,"eQTL_snpsloc.txt"))
gene_location_file_name = (paste0(file_path,"genesloc.txt"))
#covariates_file_name =(paste0(file_path,"Covariates.txt"))


# 加载基因型数据
snps = SlicedData$new() # 创建SNP文件为S4对象（S4对象是该包的默认输入方式）
snps$fileDelimiter = "\t"      # 指定SNP文件的分隔符
snps$fileOmitCharacters = "NA" # 定义缺失值
snps$fileSkipRows = 1          # 跳过第一行（适用于第一行是列名的情况）
snps$fileSkipColumns = 1       # 跳过第一列（适用于第一列是SNP ID的情况
snps$fileSliceSize = 2000      # 每次读取2000条数据
snps$LoadFile( SNP_file_name ) # 载入SNP文件


#加载协变量
#cvrt = SlicedData$new()
#cvrt$fileDelimiter = "\t"
#cvrt$fileOmitCharacters = "NA"
#cvrt$fileSkipRows = 1
#cvrt$fileSkipColumns = 1
#cvrt$fileSliceSize = 2000
#cvrt$LoadFile( covariates_file_name )


# 加载位置文件
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, sep= "",fill = TRUE,na.strings = "NA")
genepos<-genepos[complete.cases(genepos),]


##2:5;1
for(i in 1:1){
  
  cellsubtype_files<- list.files(paste0(expfiles_dir,allfiles[i]))
  
  for(j in 2:length(cellsubtype_files)){
    
    celltype_exp_dir<-paste0(expfiles_dir,allfiles[i],"/",cellsubtype_files[j])
    
    # 设置路径和模型信息
    useModel = modelLINEAR # 三种模型可选(modelANOVA or modelLINEAR or modelLINEAR_CROSS)
    
    expression_file_name<-celltype_exp_dir
    
    # 将输出文件设置为临时文件
    output_file_name_cis = tempfile()
    output_file_name_tra = tempfile() 
    
    # 定义gene-SNP associations的显著性P值
    pvOutputThreshold_cis = 1e-2
    pvOutputThreshold_tra = 1e-20
    
    # 判断gene-SNP属于cis的距离阈值
    cisDist = 1e5
    errorCovariance = numeric() # 定义误差项的协方差矩阵 (用的很少)
    
    # 加载表达数据
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000
    gene$LoadFile( expression_file_name )
    
    # eQTL分析
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
      cisDist = cisDist,
      snpspos = snpspos,
      genepos = genepos,
      pvalue.hist = "qqplot",
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)
    
    unlink(output_file_name_tra)
    unlink(output_file_name_cis)
    
    #导出文件
    resultfilename<-strsplit(cellsubtype_files[j],split = ".txt")[[1]][1]
    
    fwrite(me$cis$eqtls,paste0(result_path,allfiles[i],"/",resultfilename,"_ciseQTL.txt"),sep = "\t",quote = F,row.names = F)
    fwrite(me$trans$eqtls,paste0(result_path,allfiles[i],"/",resultfilename,"_transeQTL.txt"),sep = "\t",quote = F,row.names = F)
  }
}



