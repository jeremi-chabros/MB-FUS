library(ggplot2)
library(forcats)
library(nlme) 
library(lme4)
library(readxl)
library(ggpubr)
library(zoo)

multiplesheets <- function(fname) {
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # return data frame
  data_frame
}

setwd("/Users/haoyuren/Documents/Documents - Bardiche/research/LxBiopsy")
data_raw<-multiplesheets("BT008_LongEarly_LxBiopsy_062724.xlsx")
names(data)
data<-data_raw$`BT008NA N=34 Dec2023`
head(data)
data.used<-data[,c("F1","F2","F3",
                   "F4","F5","F6",
                   "PFS since Dx as of progression or censor date (months)",
                   "OS since Dx as of death or censor date (months)",
                   "Has subject died as of censor date? Y/N"
)]

# data.used$Group<-ifelse(data.used$`PFS since Dx as of progression or censor date (months)`<=12 &
#                           data.used$`OS since Dx as of death or censor date (months)`<=24,"Short",
#                         ifelse(data.used$`PFS since Dx as of progression or censor date (months)`<=12 &
#                                  data.used$`OS since Dx as of death or censor date (months)`<=24,"Long","Others")
#                         "Long")
data.used$`PFS since Dx as of progression or censor date (months)` <- sapply(data.used$`PFS since Dx as of progression or censor date (months)`,as.numeric)
data.used$`OS since Dx as of death or censor date (months)`<-sapply(data.used$`OS since Dx as of death or censor date (months)`,as.numeric)
data.used$Group1<-ifelse(
  data.used$`PFS since Dx as of progression or censor date (months)`<=10,"Early progressors",
  ifelse(data.used$`OS since Dx as of death or censor date (months)`>=24,"Long-term Survivors","Others")
)

# data.used$Group1<-ifelse(
#   data.used$`PFS since Dx as of progression or censor date (months)`<=10,"Early progressors","Long-term Survivor" 
# )

# data.used$Group2<-ifelse(
#   data.used$`OS since Dx as of death or censor date (months)`<=24,"<=24mo",">24mo"
# )
#table(data.used$Group1)
# table(data.used$Group2)

# View(data.used)

data.used<-data.used[! data.used$Group1=="Others",]
# View(data.used[data.used$Group1=="Early progressors",])





#longitudinal data
data.long<-data.frame(
  Fragment_Ratio=unlist(data.used[,c("F1","F2","F3",
                                  "F4","F5","F6")]),
  Group1=as.character(rep(data.used$Group1,times=6)),
  #Group2=as.character(rep(data.used$Group2,times=6)),
  Cycle=rep(c("1","2","3",
              "4","5","6"),each=length(data.used$Group1)),
  Id=as.character(rep(1:length(data.used$Group1),times=6))
)
data.long<-data.long[complete.cases(data.long$Fragment_Ratio),]
# View(data.long)
#length(unique(data.long$Id))

#group1
#p <- ggplot(data = data.long, aes(x = Cycle, y = Fragment_Ratio, 
#                                  group = Group1,col=Group1))
#p + 
  #aes(x = fct_inorder(Cycle))+
  # stat_summary(aes(group = Group1),
  #              geom = "line",
  #              fun = mean, na.rm = T,
  #              #shape = 17,
  #              size = 2) +
  # stat_summary(aes(group = Id),
  #              geom = "point",
  #              fun = mean, na.rm = T,
  #              #shape = 17,
  #              size = 1)+
 # geom_line(linetype="dashed",alpha=.5,size=0.6,aes(group = Id))+
#  geom_point(alpha=.6,size=2,aes(group = Id))+
  # geom_smooth(size=1.5,method = "lm", aes(fill=Group1)) +
 # geom_smooth(size=1.5,method = "lm", fill=NA) +
#  geom_vline(xintercept=4,linetype="dashed") +
#  scale_color_manual(values=c( "#E69F00", "skyblue")) +
#  xlab("Cycle")+
#  ylab("Fragment Ratio")+
#  guides(color = guide_legend(title = "Group")) +
#  theme_classic() +
#  theme(axis.text=element_text(size=12),
#        axis.title=element_text(size=14,face="bold"),
#        legend.title = element_text(size=14))

# #group2
# p <- ggplot(data = data.long, aes(x = Cycle, y = cfDNA_Ratio, group = Id, color=Group2))
# p + aes(x = fct_inorder(Cycle))+
#   stat_summary(aes(group = Id),
#                geom = "line",
#                fun = mean, na.rm = T,
#                #shape = 17,
#                size = 1) +
#   stat_summary(aes(group = Id),
#                geom = "point",
#                fun = mean, na.rm = T,
#                #shape = 17,
#                size = 2)+
#   geom_line(alpha=.3,size=0.5)+
#   xlab("Time")+
#   guides(color = guide_legend(title = "OS")) +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title = element_text(size=14))

#statistical test
#data.long.smaller<-data.long[ data.long$Group1=="Early progressors",]
#data.long.smaller$Cycle<-as.numeric(data.long.smaller$Cycle)
#fit<-lme(Fragment_Ratio ~ Cycle,
#         random=~1|Id,
#         data=data.long.smaller)
#result<-summary(fit)
#result


#data.long.smaller<-data.long[ data.long$Group1=="Long-term Survivors",]
#data.long.smaller$Cycle<-as.numeric(data.long.smaller$Cycle)
#fit<-lme(Fragment_Ratio ~ Cycle,
#         random=~1|Id,
#         data=data.long.smaller)
#result<-summary(fit)
#result

#trend analysis
# data.long.smaller$Time<-ifelse(data.long.smaller$Cycle=="c1",1,
#                                ifelse(data.long.smaller$Cycle=="c2",2,
#                                       ifelse(data.long.smaller$Cycle=="c3",3,
#                                              ifelse(data.long.smaller$Cycle=="c4",4,
#                                                     ifelse(data.long.smaller$Cycle=="c5",5,6))))
# )
# 
# fit<-lme(cfDNA_Ratio ~Time,
#          random=~1|Id,
#          data=data.long.smaller[data.long.smaller$Group1=="Early progressors",])
# result<-summary(fit)
# result

# fit<-lme(cfDNA_Ratio ~Time,
#          random=~1|Id,
#          data=data.long.smaller[data.long.smaller$Group1=="Long-term Survivor",])
# result<-summary(fit)
# result

# ci<-intervals(fit,which = "fixed")
# multi_full_3l<-cbind(coef(result),ci$fixed)
# multi_full_3l<-multi_full_3l[,c(1,2,4,5,6,8)]
# multi_full_3l<-data.frame(multi_full_3l)
# colnames(multi_full_3l)<-c("Effect size",	"Std.Error",	"t-statistics",
#                            "P-value",	"95%LL",	"95%UL")
# round(multi_full_3l,3)







# #peak data
# data.peak<-data.frame(
#   cfDNA_Ratio=apply(data.used[,c("c1","c2","c3",
#                                  "c4","c5","c6")],1,function (x)
#                                  {
#                                    max(x,na.rm = T)
#                                  }),
#   Group1=as.character(data.used$Group1),
#   #Group2=as.character(data.used$Group2),
#   Id=as.character(1:length(data.used$Group1))
# )
# 
# data.peak<-data.peak[! is.infinite(data.peak$cfDNA_Ratio),]
# data.peak<-data.peak[complete.cases(data.peak$cfDNA_Ratio),]
# #group1
# p <- ggplot(data.peak, aes(x=Group1, y=cfDNA_Ratio,fill=Group1)) + 
#   geom_violin() +
#   ylab("Peak cfDNA Ratio") +
#   xlab("PFS") +
#   scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
#   geom_boxplot(width=0.1) +
#   guides(color = guide_legend(title = "Subgroups")) +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title = element_text(size=14))
# p

# #group2
# p <- ggplot(data.peak, aes(x=Group2, y=cfDNA_Ratio,fill=Group2)) + 
#   geom_violin() +
#   ylab("Peak cfDNA Ratio") +
#   xlab("OS") +
#   scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
#   geom_boxplot(width=0.1) +
#   guides(color = guide_legend(title = "Subgroups")) +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         legend.title = element_text(size=14))
# p



# #statistical test
# table(data.peak$Group1)
# t.test(
#   data.peak$cfDNA_Ratio[data.peak$Group1=="Early progressors"],
#   data.peak$cfDNA_Ratio[data.peak$Group1=="Long-term Survivor"],
#   var.equal=F
# ) 


# #peak vs last time (cycle 4)
# data.used.filter<-data.used[!is.na(data.used$F4),]
# data.peak.last<-data.frame(
#   Fragment_Ratio=c(as.vector(apply(data.used.filter[,c("F1","F2","F3",
#                                  "F4")],1,function (x)
#                                  {
#                                    max(x,na.rm = T)
#                                  })), unlist(data.used.filter[,"F4"])),
#   Group1=as.character(rep(data.used.filter$Group1,times=2)),
#   Group2=rep(c("Peak cycle","Last cycle"),each=length(data.used.filter$Group1)),
#   #Group2=as.character(data.used$Group2),
#   Id=as.character(rep(1:length(data.used.filter$Group1),times=2))
# )
# data.peak.last<-data.peak.last[! is.infinite(data.peak.last$Fragment_Ratio),]
# head(data.peak.last)
# 
# data.peak.last<-data.peak.last[complete.cases(data.peak.last$Fragment_Ratio),]
# 
# data.peak.last$Group2<-relevel(as.factor(data.peak.last$Group2),ref="Peak cycle")
# 
# fit<-lme(Fragment_Ratio ~Group1*Group2,
#          random=~1|Id,
#          data=data.peak.last)
# result<-summary(fit)
# result



#peak vs last time (whatever last)
##first select patients with at least two observations
index1<-apply(data.used[,c("F1","F2","F3",
                          "F4","F5","F6")], 1, function (x) ifelse(sum(is.na(x))<=4,T,F)) #at least two observations among 6 cycles
index2<-apply(data.used[,c(
                           "F4","F5","F6")], 1, function (x) ifelse(sum(is.na(x))<=2,T,F)) #at least one observations after cycle 3
#data.used.filter<-data.used[index1 & index2,]
data.used.filter<-data.used[index1,]
##create a vector containing the last observations
data.used.filter$last_obs<-apply(data.used.filter[,c("F1","F2","F3",
                                                     "F4","F5","F6")], 1, 
                                 function (x) 
                                 {
                                   na.locf(as.vector(x),na.rm = F)[6]
                                 }
)

data.peak.last<-data.frame(
  Fragment_Ratio=c(as.vector(apply(data.used.filter[,c("F1","F2","F3",
                                                       "F4","F5","F6")],1,function (x)
                                                       {
                                                         max(x,na.rm = T)
                                                       })), unlist(data.used.filter$last_obs)),
  Group1=as.character(rep(data.used.filter$Group1,times=2)),
  Group2=rep(c("Peak cycle","Last cycle"),each=length(data.used.filter$Group1)),
  #Group2=as.character(data.used$Group2),
  Id=as.character(rep(1:length(data.used.filter$Group1),times=2))
)
data.peak.last<-data.peak.last[! is.infinite(data.peak.last$Fragment_Ratio),]
head(data.peak.last)

data.peak.last<-data.peak.last[complete.cases(data.peak.last$Fragment_Ratio),]

data.peak.last$Group2<-relevel(as.factor(data.peak.last$Group2),ref="Peak cycle")

#fit<-lme(Fragment_Ratio ~Group1*Group2,
#         random=~1|Id,
#         data=data.peak.last)
#result<-summary(fit)
#result

#stat.test <- tibble::tribble(
#  ~group1, ~group2,   ~p.adj,
  #"1",     "2", "ns",
  #"1",     "3", "ns",
#  "1",     "2", "ns",
#  "1",     "3", "*",
#  "2",     "3", "ns"
#)
#p <- ggplot(data.peak.last, aes(x=Group2, y=Fragment_Ratio,col=Group1,group=Id)) + 
  #geom_violin(draw_quantiles = c(.25, .50, .75)) +
 # geom_line(size=0.8,aes(x=Group2,y=Fragment_Ratio, group = Id)) +
  # geom_line(aes(x=date,y=value, color = adjuster, group = 1)) +
  #geom_point(size=1.3,aes(x=Group2,y=Fragment_Ratio, group = Id)) +
  #ylab("Fragment Ratio") +
  #xlab("Cycle") +
  #scale_color_manual(values=c( "#E69F00", "skyblue")) +
  #guides(color = guide_legend(title = "Group1")) +
  #theme_classic() +
  #theme(axis.text=element_text(size=12),
  #      axis.title=element_text(size=14,face="bold"),
  #      legend.title = element_text(size=14))
#p


#normalized to last cycle
names(data.used.filter)
data.peak.last<-data.frame(
  Fragment_Ratio=c(unlist(data.used.filter[,"last_obs"])/as.vector(apply(data.used.filter[,c("F1","F2","F3",
                                             "F4","F5","F6")],1,function (x)
                                             {
                                               max(x,na.rm = T)
                                             }))),
  Group1=as.character(rep(data.used.filter$Group1,times=1)),
  Group2=rep(c("Peak-last ratio"),each=length(data.used.filter$Group1)),
  #Group2=as.character(data.used.filter$Group2),
  Id=as.character(rep(1:length(data.used.filter$Group1),times=1))
)
data.peak.last<-data.peak.last[! is.infinite(data.peak.last$Fragment_Ratio),]
head(data.peak.last)

data.peak.last<-data.peak.last[complete.cases(data.peak.last$Fragment_Ratio),]


 
#stat.test <- tibble::tribble(
#  ~group1, ~group2,   ~p.adj,
  #"1",     "2", "ns",
  #"1",     "3", "ns",
#  "1",     "2", "ns",
#  "1",     "3", "*",
#  "2",     "3", "ns"
#)
#p <- ggplot(data.peak.last, aes(x=Group1, y=Fragment_Ratio,fill=Group1)) + 
  # geom_violin(draw_quantiles = c(.25, .50, .75)) +
 # geom_violin() +
#  geom_boxplot(width=0.1) +  # Adjust width and position as needed
  #facet_wrap(~ Group2, scales = "free", ncol = 1) +
  # stat_pvalue_manual(
  #   stat.test,
  #   size=10,
  #   y.position = 105, step.increase = 0.1,
  #   label = "p.adj"
  # ) +
 # ylab("Normalized Peak Fragment Ratio") +
#xlab("Group") +
#  scale_fill_manual(values=c( "#E69F00", "skyblue")) +
  # guides(color = guide_legend(title = "Group1")) +
 # theme_classic() +
  #theme(axis.text=element_text(size=12),
   #     axis.title=element_text(size=14,face="bold"),
    #    legend.position = "none")
        
#p

#t.test(data.peak.last$Fragment_Ratio[data.peak.last$Group1=="Long-term Survivors"], mu=1)
 # t.test(
#   data.peak.last$Fragment_Ratio[data.peak.last$Group1=="Early progressors"],
#   data.peak.last$Fragment_Ratio[data.peak.last$Group1=="Long-term Survivor"],
#   var.equal=F
# )
#wilcox.test(data.peak.last$Fragment_Ratio[data.peak.last$Group1=="Long-term Survivors"], mu=1)
# wilcox.test(data.peak.last$cfDNA_Ratio ~ Group1, data = data.peak.last, alternative = "two.sided")

# fisher.test(matrix(c(4, 2,1,3), ncol = 2))



#make violin plot for three groups
active_tumor<-c(2.49,	1.65,	(1.71+1.12)/2,	1.74)
data_three_groups<-data.frame(Outcome=as.numeric(c(active_tumor,data.used.filter$last_obs)),
                              Group=as.factor(c(rep("Active Tumors",4),data.used.filter$Group1)))


stat.test <- tibble::tribble(
  ~group1, ~group2,   ~p.adj,
  #"1",     "2", "ns",
  #"1",     "3", "ns",
  "Active Tumors",     "Early progressors", "ns",
  "Active Tumors",     "Long-term Survivors", "**",
  "Early progressors",     "Long-term Survivors", "**"
)
data_three_groups$Group <- factor(data_three_groups$Group,levels = c("Active Tumors","Early progressors","Long-term Survivors"))
p <- ggboxplot(data_three_groups,x="Group", y="Outcome",fill="Group",
               palette =c( "gray","#E69F00", "skyblue"),
               bxp.errorbar = TRUE, bxp.errorbar.width = 0.4,lwd=0.7) + 
  # geom_violin(draw_quantiles = c(.25, .50, .75)) +
  #facet_wrap(~ Group2, scales = "free", ncol = 1) +
  stat_pvalue_manual(
    stat.test,
    size=7,
    y.position = 3.5, step.increase = 0.1,
    label = "p.adj"
  ) +
  # geom_violin() +
  #geom_boxplot(width=1) +  # Adjust width and position as needed
  
  ylab("Last Observed Fragment Length Ratio") +
  xlab("Group") +
  # scale_fill_manual(breaks = data_three_groups$Group,
  #   values=c( "gray","#E69F00", "skyblue")) +
  # guides(color = guide_legend(title = "Group1")) +
  theme_classic() +
  theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position = "none")

p


#data_three_groups$Group<-relevel(as.factor(data_three_groups$Group),ref="Long-term Survivors")
#fit<-lm(
#  Outcome~Group,
#  data=data_three_groups
#)
#summary(fit)












# last cycle to the first

#data.first.last<-data.frame(
#  cfDNA_Ratio=c(unlist(data.used.filter[,"C4"])/unlist(data.used.filter[,"c1"])),
#  Group1=as.character(rep(data.used.filter$Group1,times=1)),
#  Group2=rep(c("Peak-last ratio"),each=length(data.used.filter$Group1)),
  #Group2=as.character(data.used.filter$Group2),
#  Id=as.character(rep(1:length(data.used.filter$Group1),times=1))
#)
#data.first.last<-data.first.last[! is.infinite(data.first.last$cfDNA_Ratio),]

#data.first.last<-data.first.last[complete.cases(data.first.last$cfDNA_Ratio),]
#head(data.first.last)

#data.first.last<-data.first.last[data.first.last$cfDNA_Ratio<3,]



#stat.test <- tibble::tribble(
#  ~group1, ~group2,   ~p.adj,
  #"1",     "2", "ns",
  #"1",     "3", "ns",
#  "1",     "2", "ns",
#  "1",     "3", "*",
#  "2",     "3", "ns"
#)
#p <- ggplot(data.first.last, aes(x=Group1, y=cfDNA_Ratio,fill=Group1)) + 
#  # geom_violin(draw_quantiles = c(.25, .50, .75)) +
#  geom_violin() +
#  geom_boxplot(width=0.1) +  # Adjust width and position as needed
#  #facet_wrap(~ Group2, scales = "free", ncol = 1) +
#  # stat_pvalue_manual(
#  #   stat.test,
#  #   size=10,
#  #   y.position = 105, step.increase = 0.1,
#  #   label = "p.adj"
#  # ) +
#  ylab("Last Cycle Normalized to First Cycle") +
#  xlab("Group") +
#  scale_fill_manual(values=c( "#E69F00", "lightblue")) +
#  # guides(color = guide_legend(title = "Group1")) +
#  theme_classic() +
#  theme(axis.text=element_text(size=12),
#        axis.title=element_text(size=14,face="bold"),
#        legend.position = "none")#

#p

#t.test(
#  data.first.last$cfDNA_Ratio[data.first.last$Group1=="Early progressors"],
#  data.first.last$cfDNA_Ratio[data.first.last$Group1=="Long-term Survivors"],
#  var.equal=F
#)#

#wilcox.test(data.first.last$cfDNA_Ratio ~ Group1, data = data.first.last, alternative = "two.sided")


