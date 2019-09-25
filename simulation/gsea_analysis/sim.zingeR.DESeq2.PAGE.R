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

data_path = "/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/data_ver9"
res_path <- "/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/results_ver11/power/GSEAMethods"
load(paste0(data_path,"/simdata.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,"simBeta..RData"))
summary_stats <- SIMUDATA$summary.zingeR.DESeq2
summary_stats$GeneID <- rownames(summary_stats)
summary_stats$GeneID <- rownames(summary_stats)
gene_sets_index = list()
gene_sets_index[[1]] <- which(Annotation>0)

##-------------------------------------------------------------
## fit PAGE (using PGSEA package)
##-------------------------------------------------------------
library(PGSEA)
Coef = data.frame(Coefficient = summary_stats$log2FoldChange)
res = PGSEA(Coef, cl=gene_sets_index, range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)
save(res,file=paste0(res_path,"/PGSEA/simres.power1000.fixannot.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".PGSEA.RData"))

