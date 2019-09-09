##-------------------------------------------------------------
## Real Data Analysis: CAMERA analysis on cell type pairs on Usoskin et al
##-------------------------------------------------------------

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

CellTypePair = c("NP1vsNF1","NP1vsNF2","NP1vsNF3","NP1vsNF4","NP1vsNF5","NP1vsNP2","NP1vsNP3","NP1vsPEP1","NP1vsPEP2","NP1vsTH")
# directories corresponding to the cell type pair
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
    for(iperm in 1:10){
        library(limma)
        library(limma)
        library(SingleCellExperiment)
        load("~/Usoskin/esetUsoskin.RData")
        cellType= droplevels(pData(eset)[,"Level 3"])
        batch = pData(eset)[,"Picking sessions"]
        genename = fData(eset)[,1]
        raw_counts = exprs(eset)
        keep = rowSums(raw_counts>0)>9
        raw_counts=raw_counts[keep,]
        genename = genename[keep]
        counts <- apply(raw_counts,  2,  function(x){
            #x <- as.numeric(x)
            storage.mode(x) <- 'integer'
            return(x) }
        )
        colnames(counts) <- paste0(colnames(raw_counts),"_",cellType)
        names(cellType) <- paste0(colnames(raw_counts),"_",cellType)
        names(batch) <- paste0(colnames(raw_counts),"_",cellType)
        rownames(counts) = genename
        counts = as.matrix(counts)
        colData <- data.frame(cell_type = cellType, batch=batch)
        design <- model.matrix(~ cellType+batch)
        res_path = "~/Usoskin/NP1vsOthers/CAMERA/cellType"
        counts = as.data.frame(counts)
        numCell = dim(counts)[2]+1
        counts$GeneID = rownames(counts)
        merge1 = merge(counts,combined_data,by = "GeneID")
        merge2 = merge1[,2:numCell]
        library(edgeR)
        lcpm <- cpm(merge2, log=TRUE)
        cor = interGeneCorrelation(lcpm,design)
        counts$GeneID <- NULL
        Annotation <- combined_data[,8:ncol(combined_data)]
        sumAnnot <- apply(Annotation, 2, sum)
        Annotation <- Annotation[, which(sumAnnot>50)] 
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
        ##-------------------------------------------------------------
        ## run CAMERA on permuted null
        ##-------------------------------------------------------------
        res = cameraPR(combined_data$stat, gene_sets_index,inter.gene.cor=cor$correlation)
        save(res,file = paste0(res_path,"/res.Uso.",celltype,".zd.pmt",iperm,".CAMERA.ModiCor.RData"))
        rm(res)
    }
    print(dim(Annotation))
    gene_sets_index2 <- list()
    for(i in 1:ncol(Annotation)){
        gene_sets_index2[[i]] <- which(Annotation[,i]>0)
        names(gene_sets_index2)[i] = colnames(Annotation)[i]
    }
    ##-------------------------------------------------------------
    ## run CAMERA on the alternative
    ##-------------------------------------------------------------
    res = cameraPR(combined_data$stat, gene_sets_index2,inter.gene.cor=cor$correlation)
    res_path = "~/Usoskin/NP1vsOthers/cellType"
    save(res,file = paste0(res_path,"/res.Uso.",celltype,".zd.CAMERA.ModiCor.RData"))
}
        
        
    

