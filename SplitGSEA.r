#### GSEA split
####
for(irpt in 1:20){
    #load(paste0(res_path,"/Data.H9vsDEC.rpt.",irpt,".split.",isplit,".res_",useMethod,".summary.ComHaoA.RData"))
    #### generate GSEA files
    for(isplit in 1:2){
        load(paste0("/net/mulan/home/yingma/iDEA/chu/dataset/H9vsDEC/split/Chu.H9vsDEC.split.test.Counts.rpt",irpt,".RData"))
        data = test_data[[isplit]]
        phenotype = data$cell_type
        count = data$counts
        load("/net/mulan/shiquans/dataset/Anotation/HaoA.GeneSets.c2.c5.c6.c7.larger20.RData")
        common = intersect(rownames(combined_pathway),rownames(count))
        Annotation = combined_pathway[which(rownames(combined_pathway) %in% common),]
        index = which(rownames(count) %in% common)
        count = count[index,]
        Annotation = Annotation[order(match(rownames(Annotation), rownames(count))),]
        sumAnnot <- apply(Annotation, 2, sum)
        # filtering out at 10 genes are annotated,
        Annotation <- Annotation[, which(sumAnnot>20)]
        genename = rownames(count)
        
        #### generate pheno
        pData = data.frame(phenotype = phenotype)
        rownames(pData) = c(paste0(1:length(phenotype),".CEL"))
        library("ArrayTools")
        pData$celltype = phenotype
        output.cls(pData,"phenotype", filename = paste0("/net/mulan/home/yingma/iDEA/chu/dataset/H9vsDEC/GSEA/simPhenotype.split",isplit))
        #### generate count
        colnames(count) = c(paste0(1:length(phenotype),".CEL"))
        counts1 = cbind(genename,count)
        eset = createExpressionSet(pData,counts1)
        output.gct(eset,filename=paste0("/net/mulan/home/yingma/iDEA/chu/dataset/H9vsDEC/GSEA/counts.rpt.",irpt,".split.",isplit,".RData"))
        #### generate Anno files
        for(iannot in 1:dim(Annotation)[2]){
            genename = rownames(Annotation)
            annoGene = which(Annotation[,iannot] == 1)
            gene=genename[annoGene]
            data=c("na",gene)
            GoGene=as.matrix(data)
            GoGene=t(GoGene)
            rownames(GoGene)=colnames(Annotation)[iannot]
            write.table(GoGene,file=paste0("/net/mulan/home/yingma/iDEA/chu/dataset/H9vsDEC/GSEA/Anno.",iannot,"rpt.",irpt,".split.",isplit,".txt"),quote=F,row.names=T,col.names=F,sep="\t")
        }
    }
}

#### Run GSEA split
args    <- as.numeric(commandArgs(TRUE))
irpt <- args[1]
isplit <- args[2]
script = paste0("java -cp ~/gsea-3.0.jar -Xmx5000m xtools.gsea.Gsea -res /net/mulan/home/yingma/iDEA/chu/dataset/H9vsDEC/GSEA2/counts.rpt.",irpt,".split.",isplit,".RData.gct -cls /net/mulan/home/yingma/iDEA/chu/dataset/H9vsDEC/GSEA2/simPhenotype.split",isplit,".cls -gmx /net/mulan/home/yingma/iDEA/chu/dataset/H9vsDEC/GSEA2/AnnoForH9vsDEC.rpt.",irpt,".split.",isplit,".gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label H9vsDEC.rpt.",irpt,".split.",isplit," -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 15000 -set_min 0 -zip_report false -out /net/mulan/home/yingma/iDEA/chu/result/H9vsDEC/GSEA/ -gui false")
system(script)

