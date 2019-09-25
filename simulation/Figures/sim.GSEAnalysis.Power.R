## functions to ge the detected DE proportio under fdr cutoff
getProp = function(pvalue,lenTrue,cutoff){
    TrueAnno  = paste0("A",c(1:lenTrue))
    order_alt = pvalue[order(pvalue,decreasing = F)]
    fdr <-c()
    for(ifdr in 1:length(order_alt)){
        fdr<-c(fdr, names(order_alt)[ifdr] %in% TrueAnno)
    }
    fdr = cumsum(1-fdr)/c(1:length(order_alt))
    names(fdr) = names(order_alt)
    Prop = sum(names(order_alt)[order_alt<cutoff] %in% TrueAnno) / lenTrue
    return(Prop)
}

##---------------------------------------------------------------------------------------
## Parameters for null ## Here we just use one simulation parameter setting as an example
##---------------------------------------------------------------------------------------
coverage_rate <- 0.1 
alpha1 <- -2
alpha2 <- 0

powerGSEAplot = function(iDEA_pval, fGSEA_pval, CAMERA_pval, PAGE_pval, GSEA_pval){
Prop_iDEA = getProp(iDEA_pval,length(iDEA_pval$alt_iDEA),0.05)
Prop_fGSEA = getProp(fGSEA_pval,length(fGSEA_pval$alt_fGSEA),0.05)
Prop_CAMERA = getProp(CAMERA_pval,length(CAMERA_pval$alt_CAMERA),0.05)
Prop_PAGE = getProp(PAGE_pval,length(PAGE_pval$alt_PAGE),0.05)
Prop_GSEA = getProp(GSEA_pval,length(GSEA_pval$alt_GSEA),0.05)

Setting = as.expression(bquote(''~tau[2]~'='~.(alpha2)~'CR = '~.(coverage_rate)))
############################################ Create data for plot
data=data.frame(Method = c("iDEA","fGSEA","CAMERA","PAGE","GSEA"),Proportion = c(Prop_iDEA,Prop_fGSEA,Prop_CAMERA,Prop_PGSEA,Prop_GSEA))
data$Method = factor(data$Method,levels = c("iDEA","fGSEA","CAMERA","PAGE","GSEA"))
data$group = paste0("Tau = ",alpha2, " CR = ",coverage_rate)

library(ggplot2)
p=ggplot(data, aes(fill=Method, y=Proportion,color = Method,x = Method)) +
geom_bar(position="dodge", stat="identity",colour = "white", size=0.5)+
geom_text(aes(label=paste(round(Proportion*100,digits= 0),"%",sep="")),position=position_dodge(width=0.9), vjust="inward",size=5,colour = "black")+
#scale_y_continuous(breaks = round(seq(0, max(plot1_1$value),by = max(plot1_1$value)/5),0))+
theme(plot.margin = margin(1, 1, 1, 1, "cm"),
#legend.position=c(0.14,0.76),
panel.background = element_blank(),
plot.background = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1),
axis.text = element_text(size = 20,color = "black"),
axis.title.y = element_text(size = 22,face="bold"),
axis.title.x = element_blank(),
axis.text.x=element_blank(),
axis.ticks.x=element_blank(),
legend.title=element_text(size = 18,face="bold"),
legend.text=element_text(size = 16),
legend.key = element_rect(colour = "transparent", fill = "white"),
legend.key.size = unit(1.5, 'lines'))+
#xlab("GSEA Methods")+
ylab("% Truth Detected")+
theme(axis.title.y = element_text(vjust=2.5))+ scale_color_manual(values=c("#F4A582","#B8E186","#80B1D3","#8dd3c7","wheat"),labels = c("iDEA","fGSEA","CAMERA","PAGE","GSEA"))+scale_fill_manual(values=c("#F4A582","#B8E186","#80B1D3","#8dd3c7","wheat"),labels = c("iDEA","fGSEA","CAMERA","PAGE","GSEA"))+theme(legend.position="bottom")+scale_y_continuous(labels = scales::percent_format(accuracy = 2))+facet_wrap(~group,labeller=label_parsed,ncol = 4)+theme(strip.text = element_text(size = 16,face = "bold")) #### change the color in scale_color_manual

return(p)
}


