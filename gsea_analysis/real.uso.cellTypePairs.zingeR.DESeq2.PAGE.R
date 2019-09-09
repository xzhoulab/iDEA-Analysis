##-------------------------------------------------------------
## Simulations: PAGE analysis on Real Dataset, Usoskin et al
##-------------------------------------------------------------
CellTypePair = c("NP1vsNF1","NP1vsNF2","NP1vsNF3","NP1vsNF4","NP1vsNF5","NP1vsNP2","NP1vsNP3","NP1vsPEP1","NP1vsPEP2")
for(icelltype in 1:length(CellTypePair)){
    res_path ="~/res.scRNAseq.Usoskin/cellType_pair"
    dirs = list.dirs(recursive=FALSE,res_path)
    dirs = dirs[grepl('.new', dirs)]
    index1=sub("~/res.scRNAseq.Usoskin/cellType_pair/(.*?).new", "\\1", dirs)
    dirs = dirs[which(index1%in%CellTypePair)]
    index1 = index1[which(index1%in%CellTypePair)]
    dir = dirs[icelltype]
    celltype = index1[icelltype]
    load(paste0(dir,"/res.uso.",celltype,".zingeR.DESseq2.summary.ComA.RData"))
    combined_data <- combined_data[!is.na(combined_data$log2FoldChange),]
    cellsplit <- strsplit(celltype, "vs")
    cell1 = cellsplit[[1]][1]
    cell2 = cellsplit[[1]][2]
    for(iperm in c(1:10)){
        library(limma)
        library(PGSEA)
        ##### Count data
        ipmt = iperm
        ##### ZingeR summary statistics
        res_path = "~/Usoskin/NP1vsOthers/PGSEA/cellType"
        Annotation <- combined_data[,8:ncol(combined_data)]
        sumAnnot <- apply(Annotation, 2, sum)
        Annotation <- Annotation[, which(sumAnnot>50)] ##
        print(dim(Annotation))
        num_gene = dim(Annotation)[1]
        gene_sets_index <- list()
        for(i in 1:ncol(Annotation)){
            GOAnnot = Annotation[,i]
            set.seed(iperm)
            GOAnnot = GOAnnot[sample(1:num_gene)]
            gene_sets_index[[i]] <- which(GOAnnot>0)
            names(gene_sets_index)[i] = colnames(Annotation)[i]
        }
        ###check the names of gene sets
        Coef = data.frame(Coefficient = combined_data[,3])
        ##-------------------------------------------------------------
        ## run PAGE on permuted null
        ##-------------------------------------------------------------
        res = PGSEA(Coef, cl=gene_sets_index, range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)
        save(res,file = paste0(res_path,"/res.Uso.",celltype,".zd.pmt",iperm,".PGSEA.RData"))
    }
    print(dim(Annotation))
    gene_sets_index2 <- list()
    for(i in 1:ncol(Annotation)){
        gene_sets_index2[[i]] <- which(Annotation[,i]>0)
        names(gene_sets_index2)[i] = colnames(Annotation)[i]
    }
    ###check the names of gene sets
    Coef = data.frame(Coefficient = combined_data[,3])
    ##-------------------------------------------------------------
    ## run PAGE on alternative
    ##-------------------------------------------------------------
    res = PGSEA(Coef, cl=gene_sets_index, range = c(0, 15000), ref = NULL, center = FALSE, p.value = TRUE, weighted = TRUE)
    res_path = "~/Usoskin/NP1vsOthers/cellType"
    save(res,file = paste0(res_path,"/res.Uso.",celltype,".zd.PGSEA.RData"))
}

