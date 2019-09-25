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

##-------------------------------------------------------------------------------------
## load DE results from zingeR
##-------------------------------------------------------------------------------------
library(limma)
data_path = "/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/data_ver9"
load(paste0(data_path,"/simdata.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".simBeta.RData"))
counts = as.matrix(SIMUDATA$counts)
library(edgeR)
lcpm <- cpm(counts, log=TRUE)
design <- SIMUDATA$design
## Estimate the gene gene correlation
cor = interGeneCorrelation(lcpm,design)
## create gene sets
summary_stats <- SIMUDATA$summary.zingeR.DESeq2
summary_stats$GeneID <- rownames(summary_stats)
summary_stats$GeneID <- rownames(summary_stats)
Annotation <- Annotation  <- SIMUDATA$annotation[,2]
colnames(Annotation) = gsub("[.]","",colnames(Annotation))
colnames(Annotation)[1] = "A0"
gene_sets_index <- list()
gene_sets_index[[1]] <- which(Annotation>0)
stats <- summary_stats$stat
names(stats) <- summary_stats$GeneID
##-------------------------------------------------------------
## fit CAMERA
##-------------------------------------------------------------
res = cameraPR(summary_stats$stat, gene_sets_index,inter.gene.cor=cor$correlation)
res_path <- "/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/results_ver11/power/GSEAMethods"
save(res,file=paste0(res_path,"/CAMERA/simres.power1000.fixannot.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".CAMERA.ModiCor.RData"))
