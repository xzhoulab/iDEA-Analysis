##-------------------------------------------------------------
## Real Data Analysis: zingeR-DESeq2 analysis on cell type pairs on Usoskin et al
##-------------------------------------------------------------
rm(list=ls())

args     <- as.numeric(commandArgs(TRUE))
ict1 <- args[1]
ict2 <- args[2]
pmt <- TRUE

## load packages
library(zinbwave)
library(BiocParallel)
library(doParallel)
library(Biobase)
library(edgeR)
library(scales)
library(DESeq2)
library(iCOBRA) # roc
library(limma)
library(genefilter) #filtered pvalues
library(MAST)
library(RColorBrewer)
library(knitr)

## load count data
load("~/Usoskin/esetUsoskin.RData")
cell_type <- droplevels(pData(eset)[,"Level 3"])
batch <- pData(eset)[,"Picking sessions"]
raw_counts <- exprs(eset)
rownames(raw_counts) <- fData(eset)$X1
keep <- rowSums(raw_counts>0)>9
raw_counts <- raw_counts[keep,]
cell_type <- as.character(cell_type)
print( table(cell_type) )

cellType <- unique(cell_type)
sel_index <- which(cell_type==cellType[ict1] | cell_type==cellType[ict2])
count_data <- raw_counts[,sel_index]
cell_type <- cell_type[sel_index]
batch_effect <- batch[sel_index]

cell_type <- as.factor(cell_type)
batch_effect <- as.factor(batch_effect)
# filtering out lowly expressed genes
FilterGenesIndex <- function(raw_count, min_reads=5){
    keep_indx <- apply(raw_count, 1, function(x){length(x[x>min_reads])>=2}  )
    return(keep_indx)
}# end function

index <- FilterGenesIndex(count_data)
count_data <- count_data[index,]
num_cell_type <- length(unique(cell_type))
iCounts <- apply(count_data,  2,  function(x){
    storage.mode(x) <- 'integer'
    return(x) }
)
rm(count_data)
names(batch_effect) <- colnames(iCounts)
num_batch_effect <- length(unique(as.numeric(as.factor(batch_effect))))

##-------------------------------------------------------------
## Fit zingeR-DESeq2 for each cell type comparisons
##-------------------------------------------------------------
res_path <- "./"
colData <- data.frame(cell_type = cell_type, batch=batch_effect)
design <- model.matrix(~ cell_type+batch_effect)
dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = colData, design = ~ cell_type+batch)
weights <- zingeR::zeroWeightsLS(counts = iCounts, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~cell_type+batch, verbose = TRUE)
assays(dse)[["weights"]] <- weights
dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
dse <- estimateDispersions(dse)
dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-(num_cell_type+num_batch_effect-1) )
save(dse,file=paste0(res_path,"/res.uso.",cellType[ict1],"vs",cellType[ict2],".zingeR.DESeq2.dse.RData") )
res <- results(dse,cooksCutoff=FALSE)
save(res, file=paste0(res_path,"/res.uso.",cellType[ict1],"vs",cellType[ict2],".zingeR.DESeq2.RData"))


