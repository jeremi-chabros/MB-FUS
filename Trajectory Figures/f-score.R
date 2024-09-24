library(ggplot2)
library(forcats)
library(nlme) 
library(lme4)
library(readxl)
library(geepack)
library(ggpubr)
library(zoo)
library(dplyr)
library(purrr)
library(broom)


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

#load data
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

#creat disease severity based on PFS and OS, seperately
data.used$`PFS since Dx as of progression or censor date (months)` <- sapply(data.used$`PFS since Dx as of progression or censor date (months)`,as.numeric)
data.used$Group1<-ifelse(
  data.used$`PFS since Dx as of progression or censor date (months)`<=10,"Bad", "Good")


data.used$Group2<-ifelse(data.used$`OS since Dx as of death or censor date (months)`<=24,"Bad","Good")


table(data.used$Group1)
table(data.used$Group2)

# View(data.used)

#filtering
## select patients with at least two observations in total and has at least one observation after time 3
index1<-apply(data.used[,c("F1","F2","F3",
                           "F4","F5","F6")], 1, function (x) ifelse(sum(is.na(x))<=4,T,F)) #at least two observations among 6 cycles
index2<-apply(data.used[,c(
  "F4","F5","F6")], 1, function (x) ifelse(sum(is.na(x))<=2,T,F)) #at least one observations after cycle 3
data.used.filter<-data.used[index1,
                            c("F1","F2","F3",
                              "F4","F5","F6","Group1","Group2")]
View(data.used.filter)

#use data.used.filter and start your code here based on group 1 and group 2 separately
### change data format
long_dat1 <- reshape(data = data.used.filter,
                    varying = c(1:6),
                    sep = "",
                    direction = "long")
outcome1 <- sapply(long_dat1[,4], as.numeric)
long_dat1$F <- outcome1
dat1 <- split(long_dat1,long_dat1$id)


### getting slope for each individual
fit_lm_to_individual <- function(data, subject_col, time_col, value_col) {
  data %>%
    group_by(!!sym(subject_col)) %>%
    group_map(~ {
      model <- lm(reformulate(time_col, response = value_col), data = .x)
      tidy(model)
    }) %>%
    bind_rows(.id = "subject")
}
result1 <- fit_lm_to_individual(long_dat1, "id", "time", "F")
temp1 <- subset(result1,term=='time')
slope1 <- temp1$estimate

### define the cutoff point
groups1 <- ifelse(slope1>0.3,1,-1)

new_d1 <- cbind(data.used.filter,groups1)


##### making 2*2 table
##### define four groups for `PFS since Dx as of progression or censor date (months)
good_negative <- sum(ifelse(new_d1$Group1=="Good" & new_d1$groups=="-1", 1, 0),na.rm = TRUE)
good_positive <- sum(ifelse(new_d1$Group1=="Good" & new_d1$groups=="1", 1, 0),na.rm = TRUE)

bad_positive <- sum(ifelse(new_d1$Group1=="Bad" & new_d1$groups=="1", 1, 0),na.rm = TRUE)
bad_negative <- sum(ifelse(new_d1$Group1=="Bad" & new_d1$groups=="-1", 1, 0),na.rm = TRUE)


t.dat <- data.frame(
  "Good" = c(good_positive, good_negative),
  "Bad" = c(bad_positive, bad_negative),
  row.names = c("Positive", "Negative"),
  stringsAsFactors = FALSE
)
colnames(t.dat) <- c("Good", "Bad")


##### define gour groups for OS since Dx as of death or censor date (months)
good_negative2 <- sum(ifelse(new_d1$Group2=="Good" & new_d1$groups=="-1", 1, 0),na.rm = TRUE)
good_positive2 <- sum(ifelse(new_d1$Group2=="Good" & new_d1$groups=="1", 1, 0),na.rm = TRUE)

bad_positive2 <- sum(ifelse(new_d1$Group2=="Bad" & new_d1$groups=="1", 1, 0),na.rm = TRUE)
bad_negative2 <- sum(ifelse(new_d1$Group2=="Bad" & new_d1$groups=="-1", 1, 0),na.rm = TRUE)


t.dat2 <- data.frame(
  "Good" = c(good_positive2, good_negative2),
  "Bad" = c(bad_positive2, bad_negative2),
  row.names = c("Positive", "Negative"),
  stringsAsFactors = FALSE
)
colnames(t.dat2) <- c("Good", "Bad")

#### binomial test for both 2*2 table
binom.test(good_negative+bad_positive,10,0.5)
binom.test(good_negative2+bad_positive2,10,0.5)

### group1 bad, group 2 good
#### change the format for ggplot of each group of interest
d11 <- subset(new_d1,Group1=="Bad")
long11 <- reshape(d11,
                 varying = c(1:6),
                 sep = "",
                 direction = "long")
long11 <- na.omit(long11)


d21 <- subset(new_d1,Group2=="Good")
long21 <- reshape(d21,
                 varying = c(1:6),
                 sep = "",
                 direction = "long")
long21 <- na.omit(long21)

#ggplot(long11,aes(x=time,y=F,group=id))+geom_point(color='pink')+geom_smooth(method='lm',se=FALSE,color='orange')+
#  theme_classic()+labs(y= "Fragment Ratio", x = "Cycle")+
#  theme(axis.title.x = element_text(size=15, face="bold"),
#        axis.title.y = element_text(size=15, face="bold"),
#        axis.text.x=element_text(size=14,face="bold"),
#        axis.text.y=element_text(size=14,face="bold"))

#ggplot(long21,aes(x=time,y=F,group=id))+geom_point(color='pink')+geom_smooth(method='lm',se=FALSE,color='royalblue')+
#  theme_classic()+labs(y= "Fragment Ratio", x = "Cycle")+
#  theme(axis.title.x = element_text(size=15, face="bold"),
#        axis.title.y = element_text(size=15, face="bold"),
#        axis.text.x=element_text(size=14,face="bold"),
#        axis.text.y=element_text(size=14,face="bold"))

#### summary of these two sets of data
da1 <- rbind(long11,long21)

names(da1) <- c('Group1','Group2','Group','time','score','id')
da1$Group <- ifelse(da1$Group==1,'EP-FLRatio','LTS-FLRatio')


