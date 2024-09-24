library(ggplot2)
#### combine the four groups of interest
alldata <- rbind(da,da1)

ggplot(alldata,aes(x=time,y=score,group=Group,color=Group))+geom_point()+geom_smooth(method='lm',se=FALSE)+
  theme_classic()+labs(y= "Ratio", x = "Cycle")+scale_color_manual(values=c("sienna3", "orange","royalblue",'turquoise1'))+
  theme(axis.title.x = element_text(size=18, face="bold"),
        axis.title.y = element_text(size=18, face="bold"),
        axis.text.x=element_text(size=18,face="bold"),
        axis.text.y=element_text(size=18,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14,face = "bold"))