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
if(!file.exists(paste0(res_path,"/fGSEA/simres.power1000.fixannot.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".fGSEA.RData"))){
    load(paste0(data_path,"/simdata.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".simBeta.RData"))
        print(irpt)
        library(fgsea)
        data(examplePathways)
        gs <- examplePathways[1]
        rm(examplePathways)
        library(DESeq2)
        summary_stats <- SIMUDATA$summary.zingeR.DESeq2
        summary_stats$GeneID <- rownames(summary_stats)
        gs[[1]] <- summary_stats$GeneID[which(SIMUDATA$annotation[,2]>0)]
        stats <- summary_stats$stat
        names(stats) <- summary_stats$GeneID
        rank_stats <- stats[order(stats,decreasing=F)]
##-------------------------------------------------------------
## fit fGSEA
##-------------------------------------------------------------
        res <- fgsea(pathways = gs, stats = rank_stats, nperm=10000)
        save(res, file=paste0(res_path,"/fGSEA/simres.power1000.fixannot.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".fGSEA.RData"))
    rm(res)
}
