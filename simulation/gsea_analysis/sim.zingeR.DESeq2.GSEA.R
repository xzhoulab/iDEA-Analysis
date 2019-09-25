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
## Create data files for GSEA
load(paste0(data_path,"/simdata.tau",alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".simBeta.RData"))
summary_stats <- SIMUDATA$summary.zingeR.DESeq2
summary_stats$GeneID <- rownames(summary_stats)
summary_stats$GeneID <- rownames(summary_stats)
stats <- summary_stats$stat
names(stats) = summary_stats$GeneID
stats <- stats[order(stats, decreasing=T)]
rank_stats = data.frame(GeneID = names(stats), metric = stats)
## Write .rnk files
write.table(rank_stats,file=paste0(res_path,"/GSEAPrerank/simres.power1000.fixannot.tau",-alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".rnk"),quote=F,sep="\t",row.names=F)
## Create .gmt files
A = Annotation
num_gene <- dim(A)[1]
genename = summary_stats$GeneID
index=which(A>0)
gene=genename[index]
data=c("na",gene)
GoGene=as.matrix(data)
GoGene=t(GoGene)
rownames(GoGene)="A1"
## Write .gmt files
write.table(GoGene,file=paste0(res_path,"/GSEAPrerank/Anno/simres.power1000.fixannot.tau",-alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".A1.rpt",irpt,".gmt"),quote=F,row.names=T,col.names=F,sep="\t")

##-------------------------------------------------------------
## Run GSEA using java script
##-------------------------------------------------------------

script = paste0("java -cp ~/gsea-3.0.jar -Xmx5000m xtools.gsea.GseaPreranked -rnk /net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/results_ver11/power/GSEAMethods/GSEAPrerank/RNK/simres.power1000.fixannot.tau",-alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt,".rnk -gmx /net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/results_ver11/power/GSEAMethods/GSEAPrerank/Anno/simres.power1000.fixannot.tau",-alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".A1.rpt",irpt,".gmt -nperm 10000 -scoring_scheme weighted -rpt_label simres.power1000.fixannot.tau",-alpha1,".tau",alpha2,".cr",coverage_rate,".fc",fold_change,".rpt",irpt," -make_sets true -rnd_seed timestamp -set_max 15000 -set_min 0 -zip_report false -out /net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_DEA/iDEA/results_ver11/power/GSEAMethods/GSEAPrerank/result/ -gui false")
system(script)

