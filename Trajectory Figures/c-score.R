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
data.used<-data[,c("c1","c2","c3",
                   "C4","c5","c6",
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
index1<-apply(data.used[,c("c1","c2","c3",
                           "C4","c5","c6")], 1, function (x) ifelse(sum(is.na(x))<=4,T,F)) #at least two observations among 6 cycles
index2<-apply(data.used[,c(
  "C4","c5","c6")], 1, function (x) ifelse(sum(is.na(x))<=2,T,F)) #at least one observations after cycle 3
data.used.filter<-data.used[index1,
                            c("c1","c2","c3",
                              "C4","c5","c6","Group1","Group2")]
View(data.used.filter)
names(data.used.filter) <- c("c1","c2","c3",
                             "c4","c5","c6","Group1","Group2")

#use data.used.filter and start your code here based on group 1 and group 2 separately
### change data format
long_dat <- reshape(data = data.used.filter,
                          varying = c(1:6),
                          sep = "",
                          direction = "long")
outcome <- sapply(long_dat[,4], as.numeric)
long_dat$c <- outcome

### getting slope for each individual
dat <- split(long_dat,long_dat$id)
fit_lm_to_individual <- function(data, subject_col, time_col, value_col) {
  data %>%
    group_by(!!sym(subject_col)) %>%
    group_map(~ {
      model <- lm(reformulate(time_col, response = value_col), data = .x)
      tidy(model)
    }) %>%
    bind_rows(.id = "subject")
}
result <- fit_lm_to_individual(long_dat, "id", "time", "c")
temp <- subset(result,term=='time')
slope <- temp$estimate

### define the cutoff point
groups <- ifelse(slope>0,1,-1)

new_d <- cbind(data.used.filter,groups)

##### making 2*2 table
##### define four groups for `PFS since Dx as of progression or censor date (months)
good_negative <- sum(ifelse(new_d$Group1=="Good" & new_d$groups=="-1", 1, 0),na.rm = TRUE)
good_positive <- sum(ifelse(new_d$Group1=="Good" & new_d$groups=="1", 1, 0),na.rm = TRUE)

bad_positive <- sum(ifelse(new_d$Group1=="Bad" & new_d$groups=="1", 1, 0),na.rm = TRUE)
bad_negative <- sum(ifelse(new_d$Group1=="Bad" & new_d$groups=="-1", 1, 0),na.rm = TRUE)


t.dat <- data.frame(
  "Good" = c(good_positive, good_negative),
  "Bad" = c(bad_positive, bad_negative),
  row.names = c("Positive", "Negative"),
  stringsAsFactors = FALSE
)
colnames(t.dat) <- c("Good", "Bad")


##### define gour groups for OS since Dx as of death or censor date (months)
good_negative2 <- sum(ifelse(new_d$Group2=="Good" & new_d$groups=="-1", 1, 0),na.rm = TRUE)
good_positive2 <- sum(ifelse(new_d$Group2=="Good" & new_d$groups=="1", 1, 0),na.rm = TRUE)

bad_positive2 <- sum(ifelse(new_d$Group2=="Bad" & new_d$groups=="1", 1, 0),na.rm = TRUE)
bad_negative2 <- sum(ifelse(new_d$Group2=="Bad" & new_d$groups=="-1", 1, 0),na.rm = TRUE)


t.dat2 <- data.frame(
  "Good" = c(good_positive2, good_negative2),
  "Bad" = c(bad_positive2, bad_negative2),
  row.names = c("Positive", "Negative"),
  stringsAsFactors = FALSE
)
colnames(t.dat2) <- c("Good", "Bad")

#### binomial test for both 2*2 table
binom.test(good_negative+bad_positive,21,0.5)
binom.test(good_negative2+bad_positive2,21,0.5)

#### change the format for ggplot of each group of interest
d1 <- subset(new_d,Group1=="Bad")
long1 <- reshape(d1,
                 varying = c(1:6),
                 sep = "",
                 direction = "long")
long1 <- na.omit(long1)
long1$c <- sapply(long1$c,as.numeric)


d2 <- subset(new_d,Group2=="Good")
long2 <- reshape(d2,
                 varying = c(1:6),
                 sep = "",
                 direction = "long")
long2$c <- sapply(long2$c,as.numeric)

#ggplot(long1,aes(x=time,y=c,group=id))+geom_point(color='pink')+geom_smooth(method='lm',se=FALSE,color='sienna3')+
#  theme_classic()+labs(y= "cfDNA Ratio", x = "Cycle")+
#  theme(axis.title.x = element_text(size=15, face="bold"),
#        axis.title.y = element_text(size=15, face="bold"),
#        axis.text.x=element_text(size=14,face="bold"),
#        axis.text.y=element_text(size=14,face="bold"))


#ggplot(long2,aes(x=time,y=c,group=id))+geom_point(color='pink')+geom_smooth(method='lm',se=FALSE,color='turquoise1')+
#  theme_classic()+labs(y= "cfDNA Ratio", x = "Cycle")+
#  theme(axis.title.x = element_text(size=15, face="bold"),
#        axis.title.y = element_text(size=15, face="bold"),
#        axis.text.x=element_text(size=14,face="bold"),
#        axis.text.y=element_text(size=14,face="bold"))


#### summary of these two sets of data
da <- rbind(long1,long2)

names(da) <- c('Group1','Group2','Group','time','score','id')
da$Group <- ifelse(da$Group==1,'EP-cfDNA','LTS-cfDNA')




