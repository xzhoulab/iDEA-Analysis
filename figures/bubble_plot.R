##-------------------------------------------------------------
## Simulation Plot: Bubble plot
##-------------------------------------------------------------

##-------------------------------------------------------------
## Create data frame for bubble plot
##-------------------------------------------------------------

load("~/Usoskin/NP1vsOthers_Bubble_plotdata")
#### bubble plot
Sig=plotdata[which(plotdata$Log>30),]
BP = Sig[which(Sig$Category=="BIOLOGICAL PROCESS"),]
OTH = Sig[which(Sig$Category=="OTHERS"),]
MF =Sig[which(Sig$Category=="MOLECULAR FUNCTION"),]
CC = Sig[which(Sig$Category=="CELLULAR COMPONENT"),]
BP = BP[order(BP$Adj_Pval,decreasing = F),]
OTH = OTH[order(OTH$Adj_Pval,decreasing = F),]
MF = MF[order(MF$Adj_Pval,decreasing = F),]
CC = CC[order(CC$Adj_Pval,decreasing = F),]

names1 = c(unique(BP$Term)[c(2,4,7,8)],
unique(CC$Term)[c(1:2)],
unique(OTH$Term)[c(1,2,4)],
unique(MF$Term)[2])
Sig = Sig[which(Sig$Term %in% names1),]
Sig = Sig[!duplicated(Sig$Term),]
unique(Sig$Term)
##-------------------------------------------------------------
## ggplot function for bubble plot
##-------------------------------------------------------------
library(ggplot2)
library("ggrepel")
p5 = ggplot(plotdata, aes(x = IDNum, y = Log,color = Category))+
geom_point(shape = 19,aes(fill = Category, size = Count,shape = 19),alpha=0.2)+
scale_radius()+
labs(x = "", y = expression(paste(bold(-log[10])," ",bolditalic(p),bold("-value"))))+
scale_size(name   = "Gene set size",
breaks = c(3000,6000,9000),
labels = c("3000","6000","9000"),
range = c(2,10))+
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
scale_color_manual(values=c("salmon","gold2","#42d4f4","#3cb44b"))+
scale_fill_manual(values = c("salmon","gold2","#42d4f4","#3cb44b"))+
theme(legend.direction = "vertical")+
theme(legend.position = c(0.20, 0.92))+
theme(legend.box = "horizontal")+
theme(legend.title.align = 0)+ylim(c(0,60))+
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

png("~/NP1vsOThers_Bubble_plot1.png",width = 12500, height = 7500, res = 600)
p5
dev.off()


