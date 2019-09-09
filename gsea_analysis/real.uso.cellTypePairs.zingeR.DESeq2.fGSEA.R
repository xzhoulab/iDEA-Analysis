##-------------------------------------------------------------
## Simulations: fGSEA analysis on Real Dataset, Usoskin et al
##-------------------------------------------------------------
rm(list=ls())
for(icelltype in 1:10){
    for(min in c(20,50)){
        library(doParallel)
        cl <- makeCluster(100)
        registerDoParallel(cl)
        res_path ="~/res.scRNAseq.Usoskin/cellType_pair"
        CellTypePair = c("NP1vsNF1","NP1vsNF2","NP1vsNF3","NP1vsNF4","NP1vsNF5","NP1vsNP2","NP1vsNP3","NP1vsPEP1","NP1vsPEP2","NP1vsTH")
        dirs = list.dirs(recursive=FALSE,res_path)
        dirs = dirs[grepl('.new', dirs)]
        index1=sub("~/res.scRNAseq.Usoskin/cellType_pair/(.*?).new", "\\1", dirs)
        dirs = dirs[which(index1%in%CellTypePair)]
        index1 = index1[which(index1%in%CellTypePair)]
        sub_path = paste0(dirs[icelltype],"/null.",min)
        files <- list.files(path=sub_path, pattern="iDEA11.RData")
        cat(paste("number of output files in this folder: ",length(files), "\n") )
        resA <- foreach(i = 1:length(files) )%dopar%{
            if(file.size(paste0(sub_path,"/",files[i]))!=0){
                load(paste0(sub_path,"/",files[i]))
                reseach <- data.frame(A=res$A[,2])}
        }
        stopCluster(cl)
        dir = dirs[icelltype]
        celltype = index1[icelltype]
        load(paste0(dir,"/res.uso.",celltype,".zingeR.DESseq2.summary.ComA.RData"))
        combined_data <- combined_data[!is.na(combined_data$log2FoldChange),]
        zingeR =combined_data$padj
        names(zingeR) = combined_data$GeneID
        combined_data <- combined_data[!is.na(combined_data$log2FoldChange),]
        # resA  here is the annotation matrix
        annotation <- as.data.frame(resA)
        #annotation<-do.call(cbind.data.frame, resA)
        summary_stats <- as.data.frame(combined_data$stat)
        summary_stats$GeneID <- combined_data$GeneID
        colnames(summary_stats) <- c("stats","GeneID")
        gene_sets <- list()
        for(i in 1:ncol(annotation)){
            gene_sets[[i]] <- summary_stats$GeneID[which(annotation[,i]>0) ]
        }
        names(gene_sets) <- paste0("A",1:ncol(annotation))
        stats <- summary_stats$stat
        names(stats) <- summary_stats$GeneID
        rank_stats <- stats[order(stats, decreasing=F)]
        library(fgsea)
        ##-------------------------------------------------------------
        ## run fGSEA under the null
        ##-------------------------------------------------------------
        res_fgsea <- fgsea(pathways = gene_sets, stats = rank_stats, nperm=10000)
        save_path <- "~/res.scRNAseq.Usoskin/fGSEA/cellType_pair"
        save(res_fgsea,file=paste0(save_path,"/annot.",min,"/res.Uso.",celltype,".null.Annot",min,".fGSEA.RData"))
        #### alternative:
        Annotation <- combined_data[,8:ncol(combined_data)]
        print(dim(Annotation))
        sumAnnot <- apply(Annotation, 2, sum)
        # filtering out at 10 genes are annotated,
        Annotation <- Annotation[, which(sumAnnot>min)]
        Annotation <- as.data.frame(Annotation)
        summary_stats <- as.data.frame(combined_data$stat)
        summary_stats$GeneID <- combined_data$GeneID
        colnames(summary_stats) <- c("stats","GeneID")
        gene_sets <- list()
        for(i in 1:ncol(Annotation)){
            gene_sets[[i]] <- summary_stats$GeneID[which(Annotation[,i]>0) ]
        }
        names(gene_sets) <- paste0("A",1:ncol(Annotation))
        stats <- summary_stats$stat
        names(stats) <- summary_stats$GeneID
        rank_stats <- stats[order(stats, decreasing=F)]
        library(fgsea)
        ##-------------------------------------------------------------
        ## run fGSEA under the alternative
        ##-------------------------------------------------------------
        res_fgsea <- fgsea(pathways = gene_sets, stats = rank_stats, nperm=10000)
        save_path <- "~/res.scRNAseq.Usoskin/fGSEA/cellType_pair"
        save(res_fgsea,file=paste0(save_path,"/annot.",min,"/res.Uso.",celltype,".alternative.Annot",min,".fGSEA.RData"))
        
    }
}
