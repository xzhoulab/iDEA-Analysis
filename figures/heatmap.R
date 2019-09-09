##-------------------------------------------------------------
## Simulation Plot: Heatmap
##-------------------------------------------------------------

#### counts
library(SingleCellExperiment)
load("~/Usoskin/esetUsoskin.RData")
cellType= droplevels(pData(eset)[,"Level 3"])
batch = pData(eset)[,"Picking sessions"]
raw_counts = exprs(eset)
raw_counts = as.data.frame(raw_counts)
genename = fData(eset)[,1]
rownames(raw_counts) = genename

NP1 = raw_counts[,which(cellType == "NP1")]
NP2 = raw_counts[,which(cellType == "NP2")]
NP3 = raw_counts[,which(cellType == "NP3")]
NF1 = raw_counts[,which(cellType == "NF1")]
NF2 = raw_counts[,which(cellType == "NF2")]
NF3 = raw_counts[,which(cellType == "NF3")]
NF4 = raw_counts[,which(cellType == "NF4")]
NF5 = raw_counts[,which(cellType == "NF5")]
PEP1 = raw_counts[,which(cellType == "PEP1")]
PEP2 = raw_counts[,which(cellType == "PEP2")]
TH = raw_counts[,which(cellType == "TH")]
#### new heatmap 2019. 07.5
select = c("Cav1","Basp1","Camk2a","Carhsp1","Arap1","A3galt2","Cpne4"
,"Chchd10","Clec2l","Clstn3","Daf2","Celf6","Ctxn1","F2rl2",
"Fam70a","Efna1","Dgki","Dgkg","Dkk3","Gm765","Gna14","Fxyd7","Gas6","Kcnq2","Kitl", "Htr1d","Iqsec1","Kcnk1","Kcnh2","Hs3st2","Lgi3","Htra1","Necab3","Mtap7d1", "Lphn1","Mcam","Ndrg2","Lynx1","Nfasc","Ly86","Agtr1a","Lpar3","Lpar5",   "Mrgprd","Nnat","P2rx3","Pkig","Plcb3","Plxnc1","Prkar2b")

counts = cbind(NP1,NP2,NP3,NF1,NF2,NF3,NF4,NF5,PEP1,PEP2,TH)
counts = as.data.frame(counts)
heatdata = counts[which(rownames(counts)%in%select),]
heatdata = heatdata[order(match(rownames(heatdata),select)),]
library(reshape2)
library(ggplot2)
##### using pheatmap
heatdata1 = log(heatdata+1,base=10)
phenotype = as.character(cellType)
phenotype = factor(phenotype,levels = c("NP1","NP2","NP3","NF1","NF2","NF3","NF4","NF5","PEP1","PEP2","TH"))

mat_col <- data.frame(group = phenotype)
colnames(mat_col) = "Cell Type"
rownames(mat_col) <- colnames(heatdata1)
mat_colors <- list(group = c("pink","springgreen4","blue","thistle1","green","coral1","burlywood1","darkorchid1","goldenrod1","lightsteelblue1","rosybrown1"))
names(mat_colors$group) <- levels(phenotype)
##-------------------------------------------------------------
## Write out the .png file
##-------------------------------------------------------------
png("~/manuscript/ver13/heatmap_uso_50.png",width = 11500, height = 9000, res = 600)
p=pheatmap(
mat               = heatdata1,
show_colnames     = FALSE,
show_rownames     = TRUE,
cluster_rows      = T,
cluster_cols      = F,
annotation_col    = mat_col,
annotation_colors = mat_colors,
drop_levels       = F,
fontsize          = 26,
)
p
dev.off()
