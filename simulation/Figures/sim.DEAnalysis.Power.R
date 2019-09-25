getFDR = function(pvalue,TrueGene,method){
    if(method == "iDEA"){
        order_alt = pvalue[order(pvalue,decreasing = T)]
    }else{
        order_alt = pvalue[order(pvalue,decreasing = F)]
    }
    fdr <-c()
    for(ifdr in 1:length(order_alt)){
        fdr<-c(fdr, names(order_alt)[ifdr] %in% TrueGene)
    }
    fdr = cumsum(1-fdr)/c(1:length(order_alt))
    names(fdr) = names(order_alt)
    return(fdr)
}
getProp = function(fdr,TrueGene,method,cutoff){
    Prop = sum(names(fdr)[which(fdr<cutoff)] %in% TrueGene) / length(TrueGene)
    return(Prop)
}
##---------------------------------------------------------------------------------------
## Parameters for alternative ## Here we just use one simulation parameter setting as an example
##---------------------------------------------------------------------------------------
coverage_rate <- 0.1 
alpha1 <- -2
alpha2 <- 1.0

## power plot results
powerDE.plot <- function(zingeR_res, MAST_res, edgeR_res, TrueGene){
METHOD = c("zingeR.DESeq2","MAST","edgeR")
Data = NULL
for(method in METHOD){
if(method == "zingeR.DESeq2"){
DE = zingeR_res$zingeR_pval;
iDEA_A = zingeR_res$iDEA_A;
iDEA_noA = zingeR_res$iDEA_noA;
}else if(method == "MAST"){
DE = MAST_res$zingeR_pval;
iDEA_A = MAST_res$iDEA_A;
iDEA_noA = MAST_res$iDEA_noA;
}else if(method == "edgeR"){
DE = edgeR_res$zingeR_pval;
iDEA_A =edgeR_res$iDEA_A;
iDEA_noA = edgeR_res$iDEA_noA;
}
FDR_DE = getFDR(DE,TrueGene,method)
FDR_iDEA = getFDR(iDEA_A,TrueGene,"iDEA")
FDR_iDEA_noA = getFDR(iDEA_noA,TrueGene,"iDEA")
Prop_DE = sum(names(FDR_DE)[which(FDR_DE<0.05)] %in% TrueGene)/length(TrueGene)
Prop_iDEA = sum(names(FDR_iDEA)[which(FDR_iDEA<0.05)] %in% TrueGene)/length(TrueGene)
Prop_iDEA_noA = sum(names(FDR_iDEA_noA)[which(FDR_iDEA_noA<0.05)] %in% TrueGene)/length(TrueGene)
            
data=data.frame(Method = c("iDEA: gene set","iDEA: no gene set",method),Proportion = c(Prop_iDEA,Prop_iDEA_noA,Prop_other))
data$Method = factor(data$Method,levels = c("iDEA: gene set","iDEA: no gene set",method))
data$summary = method
Data = rbind(Data,data)
}
library(ggplot2)
p=ggplot(Data, aes(fill=Method, y=Proportion,color = Method,x = summary)) +
geom_bar(position="dodge", stat="identity",size=0.75,aes(alpha = Method))+geom_text(aes(label=paste(round(Proportion*100,digits= 0),"%",sep="")),position=position_dodge(width=0.9), vjust="inward",size=5.5,colour = "black")+
theme(plot.margin = margin(1, 1, 1, 1, "cm"),
panel.background = element_blank(),
plot.background = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1),
axis.text = element_text(size = 18,color = "black"),
axis.title = element_text(size = 24,face="bold"),
axis.text.x=element_blank(),
axis.ticks.x=element_blank(),
legend.title=element_text(size = 20,face="bold"),
legend.text=element_text(size = 22),
legend.key = element_rect(colour = "transparent", fill = "white"),
legend.key.size = unit(1.5, 'lines'))+
xlab("Summary Statistics")+
ylab("% Truth Detected")+
theme(axis.title.x = element_text(vjust=-0.5))+
theme(axis.title.y = element_text(vjust=2.5))+ scale_color_manual(values=c("#F4A582","#F4A582","#80B1D3","#66C2A5","#BC80BD"),labels = c("iDEA: gene set","iDEA: no gene set","zingeR","MAST","edgeR"))+scale_fill_manual(values=c("#F4A582","#F4A582","#80B1D3","#66C2A5","#BC80BD"),labels = c("iDEA: gene set","iDEA: no gene set","zingeR","MAST","edgeR"))+ scale_alpha_manual(values=c(1,0,1,1,1))+ theme(legend.position="bottom")+scale_y_continuous(labels = scales::percent_format(accuracy = 2),limits = c(0,1))+ facet_wrap(~group,labeller=label_parsed,ncol = 4)+theme(strip.text = element_text(size = 22,face = "bold"))+guides(alpha = FALSE)+guides(color = FALSE)

return(p)
}

  






            



library(ggplot2)
p=ggplot(Data, aes(fill=Method, y=Proportion,color = Method,x = summary)) +
geom_bar(position="dodge", stat="identity",size=0.75,aes(alpha = Method))+geom_text(aes(label=paste(round(Proportion*100,digits= 0),"%",sep="")),position=position_dodge(width=0.9), vjust="inward",size=5.5,colour = "black")+
theme(plot.margin = margin(1, 1, 1, 1, "cm"),
panel.background = element_blank(),
plot.background = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1),
axis.text = element_text(size = 18,color = "black"),
axis.title = element_text(size = 24,face="bold"),
axis.text.x=element_blank(),
axis.ticks.x=element_blank(),
legend.title=element_text(size = 20,face="bold"),
legend.text=element_text(size = 22),
legend.key = element_rect(colour = "transparent", fill = "white"),
legend.key.size = unit(1.5, 'lines'))+
xlab("Summary Statistics")+
ylab("% Truth Detected")+
theme(axis.title.x = element_text(vjust=-0.5))+
theme(axis.title.y = element_text(vjust=2.5))+ scale_color_manual(values=c("#F4A582","#F4A582","#80B1D3","#66C2A5","#BC80BD"),labels = c("iDEA: gene set","iDEA: no gene set","zingeR","MAST","edgeR"))+scale_fill_manual(values=c("#F4A582","#F4A582","#80B1D3","#66C2A5","#BC80BD"),labels = c("iDEA: gene set","iDEA: no gene set","zingeR","MAST","edgeR"))+ scale_alpha_manual(values=c(1,0,1,1,1))+ theme(legend.position="bottom")+scale_y_continuous(labels = scales::percent_format(accuracy = 2),limits = c(0,1))+ facet_wrap(~group,labeller=label_parsed,ncol = 4)+theme(strip.text = element_text(size = 22,face = "bold"))+guides(alpha = FALSE)+guides(color = FALSE)
return(p)
}

      
