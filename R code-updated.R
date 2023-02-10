library(reader)
library(rms)
data <- read.csv("d:/data/IN.csv", header=TRUE, sep=",")
dd=datadist(data)
options(datadist="dd")
fit <- lrm(GBC~Age+Size+Course+CEA+CA199,data=data,x=T,y=T)
nom <- nomogram(fit, fun=plogis, lp=F, funlabel="Risk")

plot(nom,
     
     label.every = 1,
     
     col.conf = c('red','green'),
     
     conf.space = c(0.1,0.5),
     
     col.grid = gray(c(0.8,0.95)),
     
     which='shock')  ------#note Figure.2

fit
cal<-calibrate(fit,method='boot',B=1000)
plot(cal,xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability",ylab="Observed Probability",subtitles=FALSE)   ------#note Figure.5A CPF-C
cal

library(nomogramFormula)
option <- options(datadist = "dd")
options(option)
results <- formula_rd(nomogram = nom)
data$points <- points_cal(formula = results$formula,rd=data)
points <- points_cal(formula = results$formula,rd=data)
points

write.table(points, file='data.csv',row.names=FALSE,col.names=FALSE,sep=',')

getwd()   -----#note get CPF-C in Internal(IN) samples


library(nomogramFormula)
option <- options(datadist = "dd")
options(option)
results <- formula_rd(nomogram = nom)
dataex <- read.csv("d:/data/EX.csv", header=TRUE, sep=",")
dataex$points <- points_cal(formula = results$formula,rd=dataex)
points <- points_cal(formula = results$formula,rd=dataex)
points

write.table(points, file='dataex.csv',row.names=FALSE,col.names=FALSE,sep=',')

getwd()  -----#note get CPF-C in external(EX) samples

library(reader)
library(rms)
data1 <- read.csv("d:/data/shuju1.csv", header=TRUE, sep=",")
dd=datadist(data)
options(datadist="dd")
fit1 <- lrm(GBC~CPFA,data=data1,x=T,y=T)
cal1 <- calibrate(fit1,method='boot',B=1000)
plot(cal1,xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability",ylab="Observed Probability",subtitles=FALSE)   ------#note Figure.5A CPF-A


fit2 <- lrm(GBC~CPFB,data=data1,x=T,y=T)
cal2 <- calibrate(fit2,method='boot',B=1000)
plot(cal2,xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability",ylab="Observed Probability",subtitles=FALSE)   ------#note Figure.5A CPF-B

library(reader)
library(rms)
data2 <- read.csv("d:/data/shuju2.csv", header=TRUE, sep=",")
fita <- lrm(GBC~CPFA,data=data2,x=T,y=T)
cala <- calibrate(fita,method='boot',B=10000)
plot(cala,xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability",ylab="Observed Probability",subtitles=FALSE)   ------#note Figure.7A CPF-A

fitb <- lrm(GBC~CPFB,data=data2,x=T,y=T)
calb <- calibrate(fitb,method='boot',B=10000)
plot(calb,xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability",ylab="Observed Probability",subtitles=FALSE)   ------#note Figure.7A CPF-B

fitc <- lrm(GBC~CPFC,data=data2,x=T,y=T)
calc <- calibrate(fitc,method='boot',B=10000)
plot(calc,xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability",ylab="Observed Probability",subtitles=FALSE)   ------#note Figure.7A CPF-C


library(ResourceSelection)
data1 <- read.csv("d:/data/shuju1.csv", header=TRUE, sep=",")
model <- glm(GBC~CPFC,data=data1,family=binomial(link=logit))
hl <- hoslem.test(model$y,fitted(model),g=4)
hl    ------#note get Ph-l of CPF-C in Internal(IN) samples

model <- glm(GBC~CPFA,data=data1,family=binomial(link=logit))
hla <- hoslem.test(model$y,fitted(model),g=4)
hla    ------#note get Ph-l of CPF-A in Internal(IN) samples

model <- glm(GBC~CPFB,data=data1,family=binomial(link=logit))
hlb <- hoslem.test(model$y,fitted(model),g=4)
hlb    ------#note get Ph-l of CPF-B in Internal(IN) samples

library(ResourceSelection)
data2 <- read.csv("d:/data/shuju2.csv", header=TRUE, sep=",")
model <- glm(GBC~CPFC,data=data2,family=binomial(link=logit))
hl <- hoslem.test(model$y,fitted(model),g=4)
hlc    ------#note get Ph-l of CPF-C in EX samples

model <- glm(GBC~CPFA,data=data2,family=binomial(link=logit))
hla <- hoslem.test(model$y,fitted(model),g=4)
hla    ------#note get Ph-l of CPF-A in EX samples

model <- glm(GBC~CPFB,data=data2,family=binomial(link=logit))
hlb <- hoslem.test(model$y,fitted(model),g=4)
hlb    ------#note get Ph-l of CPF-B in EX samples



library(Hmisc) 
library(rms)
library(boot)
library(caret)
library(readr)
data1 <- read.csv("d:/data/shuju1.csv", header=TRUE, sep=",")    ------#note  IN samples/change 'shuju1.csv' to 'shuju2.csv' for EX samples
dd1=datadist(data1)
options(datadist="dd1")
formula1<-as.formula(GBC~CPFC)            ------#note  CPFA/CPFB/CPFC --> CPF-A/CPF-B/CPF-C, change if needed
fit1 <- lrm(formula1,data=data1,x=T,y=T)
data1$predvalue <- predict(fit1)
library(pROC)
modelROC <- roc(data1$GBC,data1$predvalue)
auc(modelROC)
ci(auc(modelROC))

v<-validate(fit1,method="boot",B=1000,dxy=T)
Dxy = v[rownames(v)=="Dxy",colnames(v)=="index.corrected"]
orig_Dxy =  v[rownames(v)=="Dxy",colnames(v)=="index.orig"]
bias_corrected_c_index <- abs(Dxy)/2+0.5
orig_c_index <- abs(orig_Dxy)/2+0.5
orig_c_index
bias_corrected_c_index        ------#note C-index of CPF-C in IN samples. 

library(rmda)
Data <- read.csv("d:/data/shuju1.csv", header=TRUE, sep=",")
CPFA <- decision_curve(GBC~CPFA,data = Data, family = binomial(link ='logit'),
                        thresholds= seq(0,0.9, by = 0.01),
                        confidence.intervals =0.95,study.design = 'case-control',
                        population.prevalence = 0.25)
CPFB <- decision_curve(GBC~CPFB,data = Data, family = binomial(link ='logit'),
                        thresholds= seq(0,0.9, by = 0.01),
                        confidence.intervals =0.95,study.design = 'case-control',
                        population.prevalence = 0.25)
CPFC <- decision_curve(GBC~CPFC,data = Data, family = binomial(link ='logit'),
                        thresholds= seq(0,0.9, by = 0.01),
                        confidence.intervals =0.95,study.design = 'case-control',
                        population.prevalence = 0.25)
List<- list(CPFA,CPFB,CPFC)
plot_decision_curve(List,curve.names= c('CPF-A','CPF-B','CPF-C'),
                    cost.benefit.axis =FALSE,col = c('red','blue','green'),
                    confidence.intervals =FALSE,standardize = FALSE)      -----#note  Figure.5B 

plot_clinical_impact(CPFA,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)    -----#note  Suppl. Figure.1 CPF-A

plot_clinical_impact(CPFB,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)    -----#note  Suppl. Figure.1 CPF-B

plot_clinical_impact(CPFC,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)    -----#note  Suppl. Figure.1 CPF-C

library(rmda)
Data <- read.csv("d:/data/shuju2.csv", header=TRUE, sep=",")
CPFA <- decision_curve(GBC~CPFA,data = Data, family = binomial(link ='logit'),
                        thresholds= seq(0,0.9, by = 0.01),
                        confidence.intervals =0.95,study.design = 'case-control',
                        population.prevalence = 0.25)
CPFB <- decision_curve(GBC~CPFB,data = Data, family = binomial(link ='logit'),
                        thresholds= seq(0,0.9, by = 0.01),
                        confidence.intervals =0.95,study.design = 'case-control',
                        population.prevalence = 0.25)
CPFC <- decision_curve(GBC~CPFC,data = Data, family = binomial(link ='logit'),
                        thresholds= seq(0,0.9, by = 0.01),
                        confidence.intervals =0.95,study.design = 'case-control',
                        population.prevalence = 0.25)
List<- list(CPFA,CPFB,CPFC)
plot_decision_curve(List,curve.names= c('CPF-A','CPF-B','CPF-C'),
                    cost.benefit.axis =FALSE,col = c('red','blue','green'),
                    confidence.intervals =FALSE,standardize = FALSE)      -----#note  Figure.7B 

plot_clinical_impact(CPFA,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)    -----#note  Suppl. Figure.2 CPF-A

plot_clinical_impact(CPFB,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)    -----#note  Suppl. Figure.2 CPF-B

plot_clinical_impact(CPFC,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)    -----#note  Suppl. Figure.2 CPF-C






