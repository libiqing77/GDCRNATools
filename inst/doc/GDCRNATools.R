## ----stable version, eval=FALSE, message=FALSE, warning=FALSE------------
#  ## try http:// if https:// URLs are not supported
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("GDCRNATools")

## ----development version, eval=FALSE, message=FALSE, warning=FALSE-------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("GDCRNATools", version = "devel")

## ----load, eval=TRUE, message=FALSE, warning=FALSE-----------------------
library(GDCRNATools)

## ----load data q, message=FALSE, warning=FALSE, eval=TRUE----------------
library(DT)

### load RNA counts data
data(rnaCounts)

### load miRNAs counts data
data(mirCounts)

## ----normalization q, message=FALSE, warning=FALSE, eval=TRUE------------
####### Normalization of RNAseq data #######
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

####### Normalization of miRNAs data #######
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)

## ----parse meta2 q, message=FALSE, warning=FALSE, eval=TRUE--------------
####### Parse and filter RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
datatable(as.data.frame(metaMatrix.RNA[1:5,]), extensions = 'Scroller',
        options = list(scrollX = TRUE, deferRender = TRUE, scroller = TRUE))


## ----deg q, message=FALSE, warning=FALSE, eval=TRUE----------------------
DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
datatable(as.data.frame(DEGAll), 
        options = list(scrollX = TRUE, pageLength = 5))

### All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')

### DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')

### DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')

## ----ce q, message=TRUE, warning=FALSE, eval=TRUE------------------------
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)

datatable(as.data.frame(ceOutput), 
          options = list(scrollX = TRUE, pageLength = 5))

## ----sig q, message=FALSE, warning=FALSE, eval=TRUE----------------------
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 
    & ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]

## ----edges q, message=FALSE, warning=FALSE, eval=TRUE--------------------
### Export edges
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
datatable(as.data.frame(edges), 
        options = list(scrollX = TRUE, pageLength = 5))

## ----nodes q, message=FALSE, warning=FALSE, eval=TRUE--------------------
### Export nodes
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
datatable(as.data.frame(nodes), 
        options = list(scrollX = TRUE, pageLength = 5))

## ----auto rna, eval=FALSE, message=FALSE, warning=FALSE------------------
#  project <- 'TCGA-CHOL'
#  rnadir <- paste(project, 'RNAseq', sep='/')
#  mirdir <- paste(project, 'miRNAs', sep='/')
#  
#  ####### Download RNAseq data #######
#  gdcRNADownload(project.id     = 'TCGA-CHOL',
#                 data.type      = 'RNAseq',
#                 write.manifest = FALSE,
#                 method         = 'gdc-client',
#                 directory      = rnadir)
#  
#  ####### Download mature miRNA data #######
#  gdcRNADownload(project.id     = 'TCGA-CHOL',
#                 data.type      = 'miRNAs',
#                 write.manifest = FALSE,
#                 method         = 'gdc-client',
#                 directory      = mirdir)
#  

## ----auto clinical, eval=FALSE, message=FALSE, warning=FALSE-------------
#  ####### Download clinical data #######
#  clinicaldir <- paste(project, 'Clinical', sep='/')
#  gdcClinicalDownload(project.id     = 'TCGA-CHOL',
#                      write.manifest = FALSE,
#                      method         = 'gdc-client',
#                      directory      = clinicaldir)
#  

## ----parse meta2, message=FALSE, warning=FALSE, eval=TRUE----------------
####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)

## ----parse meta3, message=FALSE, warning=FALSE, eval=TRUE----------------
####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)

## ----merge RNAseq, message=FALSE, warning=FALSE, eval=FALSE--------------
#  ####### Merge RNAseq data #######
#  rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA,
#                           path      = rnadir, # the folder in which the data stored
#                           organized = FALSE, # if the data are in separate folders
#                           data.type = 'RNAseq')
#  
#  ####### Merge miRNAs data #######
#  mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
#                           path      = mirdir, # the folder in which the data stored
#                           organized = FALSE, # if the data are in separate folders
#                           data.type = 'miRNAs')

## ----merge clinical, message=FALSE, warning=FALSE, eval=FALSE------------
#  ####### Merge clinical data #######
#  clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
#  clinicalDa[1:6,5:10]

## ----normalization, message=FALSE, warning=FALSE, eval=FALSE-------------
#  ####### Normalization of RNAseq data #######
#  rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
#  
#  ####### Normalization of miRNAs data #######
#  mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)

## ----deg, message=FALSE, warning=FALSE, eval=FALSE-----------------------
#  DEGAll <- gdcDEAnalysis(counts     = rnaCounts,
#                          group      = metaMatrix.RNA$sample_type,
#                          comparison = 'PrimaryTumor-SolidTissueNormal',
#                          method     = 'limma')

## ----data, message=FALSE, warning=FALSE, eval=TRUE-----------------------
data(DEGAll)

## ----extract, message=FALSE, warning=FALSE, eval=TRUE--------------------
### All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')

### DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')

### DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')

## ----volcano, fig.align='center', fig.width=5, message=FALSE, warning=FALSE, eval=TRUE----
gdcVolcanoPlot(DEGAll)

## ----barplot, fig.align='center', fig.height=6, message=FALSE, warning=FALSE, eval=TRUE----
gdcBarPlot(deg = deALL, angle = 45, data.type = 'RNAseq')

## ----heatmap, message=FALSE, warning=FALSE, eval=FALSE-------------------
#  degName = rownames(deALL)
#  gdcHeatmap(deg.id = degName, metadata = metaMatrix.RNA, rna.expr = rnaExpr)

## ----ce, message=FALSE, warning=FALSE, eval=FALSE------------------------
#  ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC),
#                            pc          = rownames(dePC),
#                            lnc.targets = 'starBase',
#                            pc.targets  = 'starBase',
#                            rna.expr    = rnaExpr,
#                            mir.expr    = mirExpr)

## ----ce 2, message=FALSE, warning=FALSE, eval=TRUE-----------------------
### load miRNA-lncRNA interactions
data(lncTarget)

### load miRNA-mRNA interactions
data(pcTarget)
pcTarget[1:3]

## ----ce 22, message=FALSE, warning=FALSE, eval=FALSE---------------------
#  ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC),
#                            pc          = rownames(dePC),
#                            lnc.targets = lncTarget,
#                            pc.targets  = pcTarget,
#                            rna.expr    = rnaExpr,
#                            mir.expr    = mirExpr)

## ----message=FALSE, warning=FALSE, eval=FALSE----------------------------
#  ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 &
#      ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]
#  
#  edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
#  nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
#  
#  write.table(edges, file='edges.txt', sep='\t', quote=F)
#  write.table(nodes, file='nodes.txt', sep='\t', quote=F)

## ----cor plot, eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE-------
#  gdcCorPlot(gene1    = 'ENSG00000251165',
#             gene2    = 'ENSG00000091831',
#             rna.expr = rnaExpr,
#             metadata = metaMatrix.RNA)

## ----shiny cor plot, eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
#  shinyCorPlot(gene1    = rownames(deLNC),
#               gene2    = rownames(dePC),
#               rna.expr = rnaExpr,
#               metadata = metaMatrix.RNA)

## ----survival, message=FALSE, warning=FALSE, eval=FALSE------------------
#  ####### CoxPH analysis #######
#  survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL),
#                                    method   = 'coxph',
#                                    rna.expr = rnaExpr,
#                                    metadata = metaMatrix.RNA)

## ----survival2, message=FALSE, warning=FALSE, eval=FALSE-----------------
#  ####### KM analysis #######
#  survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL),
#                                    method   = 'KM',
#                                    rna.expr = rnaExpr,
#                                    metadata = metaMatrix.RNA,
#                                    sep      = 'median')

## ----km plot, eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE----
#  gdcKMPlot(gene     = 'ENSG00000136193',
#            rna.expr = rnaExpr,
#            metadata = metaMatrix.RNA,
#            sep      = 'median')

## ----shiny km plot, eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE----
#  shinyKMPlot(gene = rownames(deALL), rna.expr = rnaExpr,
#              metadata = metaMatrix.RNA)

## ----enrichment, message=FALSE, warning=FALSE, eval=FALSE----------------
#  enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL), simplify = TRUE)

## ----enrichment data, message=FALSE, warning=FALSE, eval=TRUE------------
data(enrichOutput)

## ----go bar, fig.height=8, fig.width=15.5, message=FALSE, warning=FALSE, eval=TRUE----
gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)

## ----go bubble, echo=TRUE, fig.height=8, fig.width=12.5, message=FALSE, warning=FALSE, eval=TRUE----
gdcEnrichPlot(enrichOutput, type='bubble', category='GO', num.terms = 10)

## ----shiny pathview, message=FALSE, warning=FALSE, eval=FALSE------------
#  library(pathview)
#  
#  deg <- deALL$logFC
#  names(deg) <- rownames(deALL)

## ----pathway, message=FALSE, warning=FALSE, eval=TRUE--------------------
pathways <- as.character(enrichOutput$Terms[enrichOutput$Category=='KEGG'])
pathways

## ----shiny pathview2, eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
#  shinyPathview(deg, pathways = pathways, directory = 'pathview')

## ----sessionInfo---------------------------------------------------------
sessionInfo()

