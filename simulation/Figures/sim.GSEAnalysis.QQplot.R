##--------------------------------------------------------------------------------------
## Simulation Plot: Type I error Control
##--------------------------------------------------------------------------------------
lambdaGC = function(pvalue){
    chisq <- qchisq(1-pvalue,1)
    lambda = median(chisq)/qchisq(0.5,1)
    return(lambda)
}
##---------------------------------------------------------------------------------------
## Parameters for null ## Here we just use one simulation parameter setting as an example
##---------------------------------------------------------------------------------------
coverage_rate <- 0.1 
alpha1 <- -2
alpha2 <- 0


QQplot <- function(pvalue_iDEA, pvalue_fGSEA,pvalue_CAMERA,pvalue_PAGE,pvalue_GSEA){
### create a data frame for qqplot 
n1 = length(pvalue_iDEA); n2 = length(pvalue_fGSEA);n3 = length(pvalue_CAMERA);n4 = length(pvalue_PAGE);n5 = length(pvalue_GSEA)
lambda_iDEA = lambdaGC(pvalue_iDEA);lambda_iDEA = round(lambda_iDEA,2)
lambda_fGSEA = lambdaGC(pvalue_fGSEA);lambda_fGSEA = round(lambda_fGSEA,2)
lambda_CAMERA = lambdaGC(pvalue_CAMERA);lambda_CAMERA = round(lambda_CAMERA,2)
lambda_PAGE = lambdaGC(pvalue_PAGE );lambda_PAGE = round(lambda_PAGE,2)
lambda_GSEA = lambdaGC(pvalue_GSEA);lambda_GSEA = round(lambda_GSEA,2)
createData = function(method,pval,n,lambdagc){
            df = data.frame(
            observed = -log10(sort(pval)),
            expected = -log10(ppoints(n)),
            clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
            cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
            para = rep(paste0(method," lambda[gc]"),n),
            Setting = rep(method,n),
            lambda = rep(lambdagc,n)
            )
        }
df1 = createData("iDEA",pvalue_iDEA,n1,lambda_iDEA) ### #F4A582
df2 = createData("fGSEA",pvalue_fGSEA,n2,lambda_fGSEA)### #B8E186
df3 = createData("CAMERA",pvalue_CAMERA,n3,lambda_CAMERA)### #80B1D3
df4 = createData("PAGE",pvalue_PAGE,n5,lambda_PAGE) ###  #8dd3c7
df5 = createData("GSEA",pvalue_GSEA,n7,lambda_GSEA) ### pink2
df = rbind(df1,df2,df3,df4,df5)
labels = as.expression(c(
        bquote(''~lambda[gc]~'='~.(lambda_iDEA)),
        bquote(''~lambda[gc]~'='~.(lambda_fGSEA)),
        bquote(''~lambda[gc]~'='~.(lambda_CAMERA)),
        bquote(''~lambda[gc]~'='~.(lambda_PAGE)),
        bquote(''~lambda[gc]~'='~.(lambda_GSEA))
        ))
log10PE <- expression(paste(bold("Expected ",-log[10]," "),bolditalic(p),bold("-value")))
log10PO <- expression(paste(bold("Observed ",-log[10]," "),bolditalic(p),bold("-value")))
df$group = paste("Tau0 = ",alpha1," CR = ",coverage_rate)
        string = paste(round(coverage_rate*100,digits= 0),"%",sep="")
df$group = factor(df$group,levels = paste("Tau0 = ",alpha1," CR = ",coverage_rate), labels = as.expression(bquote(tau[0]==.(alpha1)~'CR ='~.(string))))

p = ggplot(df)+geom_ribbon(aes(x=expected, ymax=cupper, ymin=clower), fill="grey80", alpha=.6)+geom_abline(intercept = 0, slope = 1, size = 1.5,col="white")+
        geom_point(aes(expected, observed,color = Setting), size = 1,alpha = 1.0)+
        theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 26,color = "black"),
        axis.text.y = element_text(size = 26,color = "black"),
        axis.title = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size = 15),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(1.5, 'lines'))+scale_y_continuous(breaks = c(0,1,2,3,4,5))+scale_x_continuous(breaks = c(0,1,2,3,4))+scale_colour_manual(values=c("#F4A582","#B8E186","#80B1D3","#8dd3c7","wheat"),breaks=c("iDEA","fGSEA","CAMERA","PAGE","GSEA"),labels=labels)+theme(legend.position=c(0.04,0.5),
        legend.justification=c(0, 0),
        legend.direction="vertical")+scale_shape_discrete(name="para",
        breaks= "lambda")+guides(colour = guide_legend(override.aes = list(shape = 15,size = 7.5)))+theme(legend.text.align = 0)+facet_wrap(~group, labeller = label_parsed,ncol = 1)+theme(strip.text = element_text(size = 22,face = "bold"))
        }
return(p)
}
        
    

