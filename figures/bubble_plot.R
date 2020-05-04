##-------------------------------------------------------------
## Simulation Plot: Bubble plot
##-------------------------------------------------------------

##-------------------------------------------------------------
## Create data frame for bubble plot
##-------------------------------------------------------------
load("./plotdata_Figure4e.RData")

#### bubble plot
#### order the data frame by category of the gene sets
CateoryLevels = c("IMMUNOLOGIC SIGNATURES", "CHEMICAL AND GENETIC PERTURBATIONS", "GO BIOLOGICAL PROCESS",
"GO MOLECULAR FUNCTION","GO CELLULAR COMPONENT","ONCOGENIC SIGNATURES","REACTOME","KEGG","PID","BIOCARTA")
plotdata = plotdata[order(match(plotdata$Category, CateoryLevels)),]

#### select the biologically significant gene set in the bubbleplot to show
SigSelectTem = c("GSE29618_PDC_VS_MDC_DN","GO_VASCULATURE_DEVELOPMENT","GO_BLOOD_VESSEL_MORPHOGENESIS","LIM_MAMMARY_STEM_CELL_UP",
                "GO_ENDODERM_DEVELOPMENT","GO_ANCHORING_JUNCTION","ONDER_CDH1_TARGETS_2_DN","GO_EXTRACELLULAR_MATRIX","GO_ANGIOGENESIS",
                "GO_CELL_ADHESION_MOLECULE_BINDING")
Sig = plotdata[which(plotdata$Term %in% SigSelectTem),]
#### you can modify your text label here .i.e. not capitalize the gene set term 
Sig$Term = tolower(Sig$Term)

##-------------------------------------------------------------
## ggplot function for bubble plot
##-------------------------------------------------------------
library(ggplot2)
library("ggrepel")
p4 = ggplot(plotdata, aes(x = IDNum, y = Log,color = Category))+
geom_point(shape = 19,aes(fill = Category, size = Count,shape = 19),alpha=0.8)+
scale_radius()+
labs(x = "", y = expression(paste(bold(-log[10]),bold("("),bolditalic(p),bold("-value)"))))+
theme(plot.margin = margin(1, 1, 1, 1, "cm"),
axis.text.x=element_blank(),
plot.title = element_text(lineheight=.8, face="bold"),
axis.text = element_text(size = 30),
axis.line = element_line(colour = 'black'),
axis.ticks = element_line(colour = 'grey80'),
axis.title = element_text(size = 40, face = 'bold'),
axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
legend.title=element_text(size=24,face = 'bold'),
legend.text=element_text(size=24),
panel.background = element_blank(),
panel.grid.minor = element_blank(),
panel.grid.major = element_line(colour = 'white'),
plot.background = element_blank(),
legend.key = element_rect(color = "transparent", fill = "transparent"))+
geom_hline(yintercept = 1.82, col = 'black',linetype = 2,size=2)+
guides(color = guide_legend(order = 1,override.aes = list(alpha = 1,size=7)),
size = guide_legend(order = 2,override.aes = list(alpha = 1,shape=21)),
fill = FALSE)+ labs(size = "Gene set size")+
scale_color_manual(values=c("salmon","gold2","#42d4f4","#3cb44b","chocolate2","#4363d8","#bfef45","#911eb4","#f032e6","#a9a9a9"))+
scale_fill_manual(values = c("salmon","gold2","#42d4f4","#3cb44b","chocolate2","#4363d8","#bfef45","#911eb4","#f032e6","#a9a9a9"))+
theme(legend.direction = "vertical")+
theme(legend.position = c(0.30, 0.85))+
theme(legend.box = "horizontal")+
theme(legend.title.align = 0)+
geom_text_repel(
data = Sig,
aes(label = Term),
force=1.0, point.padding=unit(1.1,'lines'),
box.padding = unit(0.1, "lines"),
vjust=-1.4,
hjust = 0.2,
size = 8,
direction='y',
nudge_x=0.2,
nudge_y = 0.2,
segment.size=0.2,
col = "black")
##-------------------------------------------------------------
## Write out the .png file
##-------------------------------------------------------------

png("~/Figure4e.png",width = 12500, height = 7500, res = 600)
p4
dev.off()


