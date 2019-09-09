##-------------------------------------------------------------
## Simulations: GSEA analysis on Real Dataset, Usoskin et al
##-------------------------------------------------------------

##-------------------------------------------------------------
## Create data files for GSEA
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
        ##### Count data
        ipmt = iperm
        ##### ZingeR summary statistics
        res_path = "~/Usoskin/NP1vsOthers/GSEAPrerank/cellType"
        Annotation <- combined_data[,8:ncol(combined_data)]
        sumAnnot <- apply(Annotation, 2, sum)
        Annotation <- Annotation[, which(sumAnnot>50)] #### 11930 genes
        print(dim(Annotation))
        stats <- combined_data$stat
        names(stats) = combined_data$GeneID
        stats <- stats[order(stats, decreasing=T)]
        rank_stats = data.frame(GeneID = names(stats), metric = stats)
        ##-------------------------------------------------------------
        ## Write .rnk files
        ##-------------------------------------------------------------
        write.table(rank_stats,file=paste0(res_path,"/Uso.",celltype,".pmt",iperm,".rnk"),quote=F,sep="\t",row.names=F)
        print(dim(Annotation))
        A = Annotation
        num_gene <- dim(A)[1]
        genename = combined_data$GeneID
        for (iannot in 1:dim(A)[2]){
            GOAnnot <- A[,iannot]
            set.seed(iperm)
            GOAnnot_perm<- GOAnnot[sample(1:num_gene)]
            index=which(GOAnnot_perm>0)
            gene=genename[index]
            data=c("na",gene)
            GoGene=as.matrix(data)
            GoGene=t(GoGene)
            rownames(GoGene)=colnames(A)[iannot]
            ##-------------------------------------------------------------
            ## Write .gmt files
            ##-------------------------------------------------------------
            write.table(GoGene,file=paste0(res_path,"/Annopmt/Uso.",celltype,".Anno",iannot,".pmt",iperm,".gmt"),quote=F,row.names=T,col.names=F,sep="\t")
            
        }
    }
    ##-------------------------------------------------------------
    ## Create data files for GSEA under the alternative
    ##-------------------------------------------------------------

res_path="~/Usoskin/NP1vsOthers/GSEAPrerank/cellType"
    stats <- combined_data$stat
    names(stats) = combined_data$GeneID
    stats <- stats[order(stats, decreasing=T)]
    rank_stats = data.frame(GeneID = names(stats), metric = stats)
    write.table(rank_stats,file=paste0(res_path,"/Uso.",celltype,".rnk"),quote=F,sep="\t",row.names=F)
    print(dim(Annotation))
    A = Annotation
    num_gene <- dim(A)[1]
    genename = combined_data$GeneID
    for (iannot in 1:dim(A)[2]){
        GOAnnot <- A[,iannot]
        GOAnnot_perm<- GOAnnot
        index=which(GOAnnot_perm>0)
        gene=genename[index]
        data=c("na",gene)
        GoGene=as.matrix(data)
        GoGene=t(GoGene)
        rownames(GoGene)=colnames(A)[iannot]
        write.table(GoGene,file=paste0(res_path,"/Anno/Uso.",celltype,".Anno",iannot,".txt"),quote=F,row.names=T,col.names=F,sep="\t")
        
    }
}

##-------------------------------------------------------------
## Run GSEA under the alternative
##-------------------------------------------------------------
args    <- as.numeric(commandArgs(TRUE))
icelltype <- args[1]
iperm <- args[2]
script1 = paste0("java -cp ~/gsea-3.0.jar -Xmx8000m xtools.gsea.GseaPreranked -rnk ~/Usoskin/NP1vsOthers/GSEAPrerank/cellType/Uso.",celltype,".rnk -gmx ~/Usoskin/NP1vsOthers/GSEAPrerank/cellType/Uso.",celltype,".gmt -nperm 1000 -scoring_scheme weighted -rpt_label Uso.",celltype," -make_sets true -rnd_seed timestamp -set_max 15000 -set_min 0 -zip_report false -out ~/Usoskin/NP1vsOthers/GSEAPrerank/cellType/ -gui false")
system(script1)

##-------------------------------------------------------------
## Run GSEA under the null
##-------------------------------------------------------------
script2 = paste0("java -cp ~/gsea-3.0.jar -Xmx5000m xtools.gsea.GseaPreranked -rnk ~/Usoskin/NP1vsOthers/GSEAPrerank/cellType/Uso.",celltype,".pmt",iperm,".rnk -gmx ~/Usoskin/NP1vsOthers/GSEAPrerank/cellType/Uso.",celltype,".pmt",iperm,".gmt -nperm 1000 -scoring_scheme weighted -rpt_label Uso.",celltype,".pmt.",iperm," -make_sets true -rnd_seed timestamp -set_max 15000 -set_min 0 -zip_report false -out ~/Usoskin/NP1vsOthers/GSEAPrerank/cellType/ -gui false")
system(script2)


    
