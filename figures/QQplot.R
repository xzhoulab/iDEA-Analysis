##-------------------------------------------------------------
## Simulation Plot: QQ plot
##-------------------------------------------------------------
##-------------------------------------------------------------
## load the data
##-------------------------------------------------------------
##### iDEA and fgsea
res_path = "~/res.Uso"
load(paste0(res_path,"/res.uso.NP1vsothers.null.Louis.iDEA11.RData"))
pvalue_iDEA = results$pvalue_louis
rm(results)
load(paste0(res_path,"/res.uso.NP1vsothers.null.fGSEA.RData"))
pvalue_fGSEA = res_fgsea$pval
##### CAMERA
load("~/Usoskin/NP1vsOthers/CAMERA/res.Uso.NP1vsOthers.zd.null.CAMERA.ModiCor.pvalue.RData")
pvalue_CAMERA = Pvalue
##### PGSEA
load("~/Usoskin/NP1vsOthers/PGSEA/res.Uso.NP1vsOthers.zd.null.PGSEA.pvalue.RData")
pvalue_PGSEA = Pvalue
###### GSEA
load("~/Usoskin/NP1vsOthers/GSEAPrerank/res.Uso.NP1vsOthers.zd.null.GSEA.pvalue.RData")
pvalue_GSEA = combined_GSEA[-which(combined_GSEA==0)]

n1 = length(pvalue_iDEA)
n2 = length(pvalue_fGSEA)
n3 = length(pvalue_CAMERA)
n5 = length(pvalue_PGSEA)
n7 = length(pvalue_GSEA)
lambdaGC = function(pvalue){
    chisq <- qchisq(1-pvalue,1)
    lambda = median(chisq)/qchisq(0.5,1)
    return(lambda)
}
lambda_iDEA = lambdaGC(pvalue_iDEA);lambda_iDEA = round(lambda_iDEA,2)
lambda_fGSEA = lambdaGC(pvalue_fGSEA);lambda_fGSEA = round(lambda_fGSEA,2)
lambda_CAMERA = lambdaGC(pvalue_CAMERA);lambda_CAMERA = round(lambda_CAMERA,2)
lambda_PGSEA = lambdaGC(pvalue_PGSEA);lambda_PGSEA = round(lambda_PGSEA,2)
lambda_GSEA = lambdaGC(pvalue_GSEA);lambda_GSEA = round(lambda_GSEA,2)
##-------------------------------------------------------------
## Create data.frame for qqplot
##-------------------------------------------------------------
createData = function(method,pval,n,lambdagc){
    df = data.frame(
    observed = -log10(sort(pval)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
    para = rep(paste0(method," lambda[gc]"),n),
    Method = rep(method,n),
    lambda = rep(lambdagc,n)
    )
}
df1 = createData("iDEA",pvalue_iDEA,n1,lambda_iDEA) ### #F4A582
df2 = createData("fGSEA",pvalue_fGSEA,n2,lambda_fGSEA)### #B8E186
df3 = createData("CAMERA",pvalue_CAMERA,n3,lambda_CAMERA)### #80B1D3
df5 = createData("PGSEA",pvalue_PGSEA,n5,lambda_PGSEA) ###  #8dd3c7
df7 = createData("GSEA",pvalue_GSEA,n7,lambda_GSEA) ### pink2
df = rbind(df1,df2,df3,df5,df7)
labels = as.expression(c(
bquote(' iDEA'),
bquote(' fGSEA'),
bquote(' CAMERA'),
bquote(' PAGE'),
bquote(' GSEA')
))
log10PE <- expression(paste(bold("Expected ",-log[10]," "),bolditalic(p),bold("-value")))
log10PO <- expression(paste(bold("Observed ",-log[10]," "),bolditalic(p),bold("-value")))
##-------------------------------------------------------------
## ggplot function for qqplot
##-------------------------------------------------------------
library(ggplot2)
p2=ggplot(df)+geom_ribbon(aes(x=expected, ymax=cupper, ymin=clower), fill="grey80", alpha=.6)+geom_abline(intercept = 0, slope = 1, size = 1.5,col="white")+
geom_point(aes(expected, observed,color = Method), size = 1,alpha = 1.0)+
xlab(log10PE)+
ylab(log10PO)+
theme(plot.margin = margin(1, 1, 1, 1, "cm"),
panel.background = element_blank(),
plot.background = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1),
axis.text = element_text(size = 20,color = "black"),
axis.title = element_text(size = 22,face="bold"),
#legend.title=element_text(size = 16,face="bold"),
legend.title=element_blank(),
legend.text=element_text(size = 14.5),
legend.key = element_rect(colour = "transparent", fill = "white"),
legend.key.size = unit(1.5, 'lines'))+
theme(axis.title.x = element_text(vjust=-0.5))+theme(axis.title.y = element_text(vjust=3.0))+scale_y_continuous(breaks = c(0,1,2,3,4,5))+scale_colour_manual(values=c("#F4A582","#B8E186","#80B1D3","#8dd3c7","wheat"),breaks=c("iDEA","fGSEA","CAMERA","PGSEA","GSEA"),labels=labels)+theme(legend.position=c(0.03,0.61),
legend.justification=c(0, 0),
legend.direction="vertical")+scale_shape_discrete(name="para",
breaks= "lambda")+guides(colour = guide_legend(override.aes = list(shape = 15,size = 7.5)))+theme(legend.text.align = 0)
##-------------------------------------------------------------
## Write out the .png file
##-------------------------------------------------------------

png("~/manuscript/ver8/Compare/Uso/GSEAMethods_Uso_NP1vsOthers.qqplot.bold.png",width = 3500, height = 3500, res = 600)
p2
dev.off()

