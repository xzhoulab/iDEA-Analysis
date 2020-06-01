##-------------------------------------------------------------------------------------
## Simulations: Simulate Count data and get summary statistics
##-------------------------------------------------------------------------------------

##-------------------------------------------------------------------------------------
## Set parameters ## Here we just use one simulation parameter setting as an example
##-------------------------------------------------------------------------------------
coverage_rate <- 0.1 
alpha1 <- -2
alpha2 <- 5.0
sigmabeta2 <- 10
se <- 0.25

##### load library here
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
library(gamlss)
library(gamlss.tr)
library(zingeR)

##-------------------------------------------------------------
## Simulation
##-------------------------------------------------------------
##### load real data set: estimate paramters from real dataset Chu et al.
counts <-read.csv("/net/mulan/shiquans/dataset/SingleCellRNAseq/Chuetal/GSE75748_sc_cell_type_ec.csv.gz",header=T,row.names=1)
str_split <- strsplit(colnames(counts),"_")
cell_type <- c()
for(i in 1:length(str_split)){
    cell_type <- c(cell_type, str_split[[i]][1])
}
sel_index <- c(which(cell_type=="EC"),which(cell_type=="TB"))
Counts <- counts[, sel_index]
cellType <- as.factor(cell_type[sel_index])
##### Filtering the counts
Counts <- ceiling( Counts[rowSums(Counts>5)>2,] )
##### get paramters from real count
paramsChu = getDatasetZTNB(counts=Counts,design=model.matrix(~cellType))
##### create simulation data list
sim.ChuData <- list()
sim.ChuData$cellType <- as.factor(cellType)
sim.ChuData$counts <- Counts
sim.ChuData$params <- paramsChu
sim.ChuData$libSizes <- colSums(Counts)
rm(Counts)
num_gene <- 10000 ##### number of genes in this simulated dataset
data_path = "/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/data_ver9"
SIMUDATA <- list()
num_cell <- length(sim.ChuData$libSizes)
num_gene <- 10000

#### function to simulate beta effect size
simulate_Beta<-function(num_gene,CR,alpha1,alpha2,se_betahat,sigmabeta2){
    annot1 <- round(num_gene*CR)
    annot0 <- num_gene - annot1
    A <- matrix(1, num_gene, 2)
    annot <- c(rep(1, annot1), rep(0, annot0))
    A[,2] <- (annot-mean(annot) )
    prob <- exp(alpha1+A[,2]*alpha2)/(1+exp(alpha1+A[,2]*alpha2))
    causal <- rbinom(num_gene, 1, prob)
    causal[is.na(causal)] <- 0
    DEind <- which(causal==1)
    beta_true    <- rep(0,num_gene)
    for (i in DEind){
        beta_true[i] = rnorm(1,mean=0,sd =sqrt(sigmabeta2)*se_betahat[i])
    }
    beta_est = beta_true + rnorm(num_gene,mean = 0, sd = se_betahat)
    res = list(betaEst = beta_est, SE = se_betahat,Annot = A,TrueDE = DEind,beta_true = beta_true)
    return(res)
}
se_betahat = rep(se,times = num_gene)
Beta <- simulate_Beta(num_gene,coverage_rate,alpha1,alpha2,se_betahat,sigmabeta2)
A = Beta$Annot
DEind = Beta$TrueDE
fcSim_all = exp(Beta$beta_true)

##### the fold change of DE genes
fcSim = fcSim_all[DEind]

#### simulate the counts using zero-truncated negative binomial
simData <- NBsimSingleCell(foldDiff = fcSim, ind = DEind,
dataset = sim.ChuData$counts, nTags = num_gene,
group = sim.ChuData$cellType,
verbose = TRUE, params = sim.ChuData$params,
lib.size = sim.ChuData$libSizes)

##### to avoid NA
if(sum(is.na(simData$counts))>0){
    simData$counts[is.na(simData$counts)] <- 0
}

rownames(simData$counts) <- paste0("gene", 1:num_gene)
colnames(simData$counts) <- paste0("cell", 1:num_cell)

SIMUDATA <- simData
SIMUDATA$annotation <- A
SIMUDATA$DEind <- DEind
SIMUDATA$libSizes <- sim.ChuData$libSizes
SIMUDATA$cellType <- sim.ChuData$cellType
SIMUDATA$foldchange <- fcSim
pheno <- sim.ChuData$cellType
SIMUDATA$pheno <- pheno
data_path = "/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/data_ver9"
save(SIMUDATA, file=paste0(data_path,"/simdata.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".simBeta.RData"))

