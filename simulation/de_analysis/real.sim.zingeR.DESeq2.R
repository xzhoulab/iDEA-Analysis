##-------------------------------------------------------------------------------------
## Simulations: zingeR analysis on simulated data
##-------------------------------------------------------------------------------------

##-------------------------------------------------------------------------------------
## Set parameters ## Here we just use one simulation parameter setting as an example
##-------------------------------------------------------------------------------------
coverage_rate <- 0.1 
alpha1 <- -2
alpha2 <- 1.0
sigmabeta2 <- 10
se <- 0.25

##### load library here
library(zinbwave)

##-------------------------------------------------------------
## fit zingeR-DESeq2
##-------------------------------------------------------------
## load simulated data
data_path = "/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/data_ver9"
load(paste0(data_path,"/simdata.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".simBeta.RData"))
iCounts <- apply(SIMUDATA$counts,  2,  function(x){
    #x <- as.numeric(x)
    storage.mode(x) <- 'integer'
    return(x) }
)
rm(simData)
colData <- data.frame(pheno = pheno)
design <- model.matrix(~pheno)
dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = colData, design = ~pheno)
weights <- zingeR::zeroWeightsLS(counts = iCounts, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~pheno, verbose = TRUE)
assays(dse)[["weights"]] <- weights
dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
dse <- estimateDispersions(dse)
dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-2)
res.ZingeR.DESeq2 = results(dse,cooksCutoff=FALSE)
SIMUDATA$summary.zingeR.DESeq2 <- res.ZingeR.DESeq2
data_path = "/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/data_ver9"

save(SIMUDATA, file=paste0(data_path,"/simdata.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".simBeta.RData"))
