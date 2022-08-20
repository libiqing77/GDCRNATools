##' @title Merge RNA/miRNAs raw counts data
##' @description Merge raw counts data that is downloaded from GDC to a 
##'     single expression matrix
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @param path path to downloaded files for merging
##' @param data.type one of \code{'RNAseq'} and \code{'miRNAs'}
##' @param organized logical, whether the raw counts data have already
##'     been organized into a single folder (eg., data downloaded by the
##'     'GenomicDataCommons' method are already organized). 
##'     Default is \code{FALSE}.
##' @return A dataframe or numeric matrix of raw counts data with rows 
##'     are genes or miRNAs and columns are samples
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Merge RNA expression data #######
##' metaMatrix <- gdcParseMetadata(project.id='TARGET-RT', 
##'     data.type='RNAseq')
##' \dontrun{rnaExpr <- gdcRNAMerge(metadata=metaMatrix, path='RNAseq/', 
##'     data.type='RNAseq')}
gdcRNAMerge <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=F, RNA_type=F){

  #通过合并path,还有sample sheet前两列得到每一个文件的完整路径
  filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                         fsep = .Platform$file.sep)
  #判断需要合并的是什么数类型，如果是RNAseq执行下面的代码
  if (data.type=='RNAseq') {
    message ('############### Merging RNAseq data ################\n',
             '### This step may take a few minutes ###\n')
    #根据需要的RNAseq表达谱数据类型，提取相应列的数据
    if(mRNA_expr_type=="STAR"){
      column=4
    }else if(mRNA_expr_type=="TPM"){
      column=7
    }else if(mRNA_expr_type=="FPKM"){
      column=8
    }else if(mRNA_expr_type=="FPKM_UQ"){
      column=9
    }
    #通过lapply循环去读每一个样本的表达，然后通过cbind合并成矩阵
    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) 
      read.table(fl,skip=6,sep="\t")[,column]))
    #获取第一个文件的第一列作为矩阵的行名
    ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
    gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
    type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
    #排除掉_PAR_Y为后缀的转录本
    index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
    rnaMatrix=rnaMatrix[index,]
    #去掉Ensembl ID后面的.和数字，eg.ENSG00000000003.13
    rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
    gene_symbol=gene_symbol[index]
    type=type[index]
    #将sample sheet的sample id这一列作为表达矩阵的列名
    colnames(rnaMatrix) <- metadata$sample

    
    #统计样本数和基因数
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    if(RNA_type){
      rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
    }
    
    if(symbol){
      rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
    }
    #输出样本数和基因数
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    #返回最后的基因表达矩阵
    return (rnaMatrix)
    
  }else if (data.type=='miRNAs') { #如果需要合并的是miRNA的数据，执行下面代码
           message ('############### Merging miRNAs data ###############\n')
        
        mirMatrix <- lapply(filenames, function(fl) cleanMirFun(fl))
        #mirs <- sort(unique(names(unlist(mirMatrix))))
        mirs <- rownames(mirbase)
        mirMatrix <- do.call('cbind', lapply(mirMatrix, 
            function(expr) expr[mirs]))
        
        rownames(mirMatrix) <- mirbase$v21[match(mirs,rownames(mirbase))]
        colnames(mirMatrix) <- metadata$sample
        
        mirMatrix[is.na(mirMatrix)] <- 0
        
        nSamples = ncol(mirMatrix)
        nGenes = nrow(mirMatrix)
        
        message (paste('Number of samples: ', nSamples, '\n', sep=''),
            paste('Number of miRNAs: ', nGenes, '\n', sep=''))
        
        return (mirMatrix)
    }else{  #如果data.type不是上面提到的两种，就报错，停止执行
    stop('data type error!')
  }
}

                                         
                                             
cleanMirFun <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    
    mirs <- unlist(lapply(strsplit(expr$Group.1, ',', fixed=TRUE),
        function(mir) mir[2]))
    
    expr <- expr[,-1]
    names(expr) <- mirs
    #rownames(expr) <- mirs
    return(expr)
}
