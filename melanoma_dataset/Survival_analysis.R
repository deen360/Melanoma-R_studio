library(boot)
library(survival)
library(ggplot2)

data(melanoma)
summary(melanoma)

data <- melanoma

#View(data)


#converting sex to factors 
data$sex <- as.factor(data$sex)
#levels(data$sex) <- c("female","male")


#converting  ulcer to factors 
data$ulcer <- as.factor(data$ulcer)
#levels(data$ulcer) <- c("absent","present")


#converting  status to factors 
data$status <- as.factor(data$status)
#levels(data$status) <- c("died from melonema","still alive", "died from other causes")

#box plot showing the age and the sex 
Sex_Age_plot <- ggplot(data, aes(x=sex, y=age, fill=sex)) + stat_boxplot(geom ='errorbar') +
  geom_boxplot()
Sex_Age_plot

#locating the outlier of the box plot 
out <-boxplot.stats(data$age)$out
out
out_ind <- which(data$age %in% c(out))
data[out_ind,] #there is a female aged 4 which is an outlier


#proportion of status on pie chart
ggplot(data, aes(x="", y=status, color=status)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + theme_void() +
  scale_color_manual(labels = c("died from melanoma","still alive", "died from other causes "),
                     values = c("red", "green", "blue"))  

#factors sex status ulcer
ggplot(data, aes(x=sex, y=ulcer, color=ulcer)) + 
  geom_bar(stat="identity", width=0.5,) +
  scale_color_manual(labels = c("abscent", "present"),
                     values = c("orange", "blue")) 


# A basic scatterplot showing the age and time, between male and females 
ggplot(data, aes(x=time, y=age, color=sex)) + 
  geom_point(size=2) + scale_color_manual(labels = c("female", "male"),
                                          values = c("orange", "blue")) +
   geom_smooth(method = "loess", se = FALSE )


# A basic scatterplot of the thickness and time showing the different status
ggplot(data, aes(x=time, y=thickness, color=status)) + 
  geom_point(size=2) + 
  scale_color_manual(labels = c("died from melanoma","still alive", "died from other causes "),
                     values = c("red", "green", "blue"))+
  geom_smooth(method = "loess", se = FALSE )
 

#creating a new column where 1 = event and 0 = censored
# still alive and died form other causes are taken as censored observations 
recodestatus<-function(x){
  if(x==1){rs=1} ## event happened 
  if(x==2){rs=0} ## no event / censored
  if(x==3){rs=0} ## no event / censored
  return(rs)
}
for(i in 1:length(data$status)){
  data$recodedStatus[i]<-recodestatus(data$status[i])
}

summary(data)

#loading of data into surv function 
mySurvival<-Surv(time=data$time, event = data$recodedStatus)


#fitting into the survfit
## single survival curve: no comparisons
fit_surv<-survfit(mySurvival~1)

plot(fit_surv,xlab="Time", ylab="Probability of survival",mark=3)
title("Kaplan meier  survival (General)")


#checking for the median survival time
fit_surv

#showing the median on the chart where p(0.5)
abline(h=0.5)


#sex
## survival curve: for sex
fit_surv_sex<-survfit(mySurvival~ data$sex)

#showing the median survival time of sex each sex
fit_surv_sex

#with confidence interval
plot(fit_surv_sex,xlab="Time (days)", ylab="Probability of survival", conf.int =TRUE, col=c("orange","blue"))
legend("bottomleft", c("female","male"), col=c("orange","blue"), lty=1, title ="Sex")
title("Kaplan meier  survival (sex)")


#without confidence interval
plot(fit_surv_sex,xlab="Time (days)", ylab="Probability of survival", conf.int = "none", col=c("orange","blue"),mark=3)
legend("bottomleft", c("female","male"), col=c("orange","blue"), lty=1,title ="Sex")
title("Kaplan meier  survival (sex)")

#showing the median on the chart for sex 
abline(h=0.5)

## now to check the survival difference in sex  (log rank test)
## Q: Is it better by chance, or statistically significant?
survdiff(mySurvival~data$sex) #p < 5% , reject null hypothesis and say sex has an effect on survival


#ulcer
## single survival curve: for ulcer
fit_surv_ulcer<-survfit(mySurvival~ data$ulcer)

#checking for median survival time of ulcer 
summary(fit_surv_ulcer)

#with confidence interval
plot(fit_surv_ulcer,xlab="Time (days)", ylab="Probability of survival", conf.int =TRUE, col=c("green","red"))
legend("bottomleft", c("absent","present"), col=c("green","red"), lty=1,title ="Ulcer")

#without confidence interval and the mark time=T
plot(fit_surv_ulcer,xlab="Time (days)", ylab="Probability of survival", conf.int ="none", col=c("green","red"),mark=3)
legend("bottomleft", c("absent","present"), col=c("green","red"), lty=1,title ="Ulcer")


#showing the median on the chart for ulcer
#only people with ulcer have median survival time 
abline(h=0.5)
abline(v=3042) 


## now to check the survival difference in ulcer  (log rank test)
survdiff(mySurvival~data$ulcer) #p < 5% , reject null hypothesis and say ulcer has an effect on survival



#cox proportional hazard model 

#model 1
cox_mod_1<-coxph(Surv(time, recodedStatus) ~ data$sex+data$age+data$year+data$thickness+data$ulcer, data = data)
summary(cox_mod_1)

#model 2 (- year)
cox_mod_2<-coxph(Surv(time, recodedStatus) ~ data$sex+data$age+data$thickness+data$ulcer, data = data)
summary(cox_mod_2)

#likely hood ratio test to see if we can drop year 
anova(cox_mod_1,cox_mod_2) #looks similar, no significant difference between model 1 and 2

#model 3 (- year,age)
cox_mod_3<-coxph(Surv(time, recodedStatus) ~ sex + thickness + ulcer, data = data)
summary(cox_mod_3)

#likely hood ratio test to see if we can drop age
#p values are close, so we drop age
anova(cox_mod_2,cox_mod_3) #looks similar, no significant difference #looks similar, no significant difference between model 2 and 3

#model 4 (- year,age,sex)
#p values are close, so we drop sex 
#sex is a factor but its just 2 levels, so we can use the wild test to check
cox_mod_4<-coxph(Surv(time, recodedStatus) ~ data$thickness+data$ulcer, data = data)
summary(cox_mod_4)
anova(cox_mod_3,cox_mod_4) #looks similar, no significant difference between model 3 and 4

#model 5 (- year,age,sex,thickness)
cox_mod_5<-coxph(Surv(time, recodedStatus) ~ ulcer, data = data)
summary(cox_mod_5)
anova(cox_mod_4,cox_mod_5)#THERE IS A DIFFERENCE BETWEN MODELS


#model 6
# the p values are far an thickness has an impact 
cox_mod_6<-coxph(Surv(time, recodedStatus) ~ data$thickness, data = data)
summary(cox_mod_6)#THERE IS A DIFFERENCE BETWEN MODELS
anova(cox_mod_4,cox_mod_6)




#model 4 was chosen to be the best (the model the thickness and ulcer)
# we do the verification 
bestfit<- cox_mod_4
summary(bestfit)
#the residuals for deviance and martingale
mresid_deviance <- resid(bestfit, type="deviance")
mresid_martingale <- resid(bestfit, type = "martingale")

summary(mresid_deviance)
summary(mresid_martingale) 

#plotting the deviance residuals with time 
plot(data$time, mresid_deviance, main = "lowess(deviance residuals)", pch = 1, col=c("green","red"))
lines(lowess(data$time, mresid_deviance))
legend("topright", c("censored","event"), col=c("green","red"), lty=1)

#deviance check
#all values are with the 2.5 range
idx.outliers <- which(residuals(bestfit,"deviance") > 2.5)
data[idx.outliers,]


#plotting the martingale residuals with time
plot(data$time, mresid_martingale, main = "lowess(martingale residuals)", pch = 1, col=c("green","red"))
lines(lowess(data$time, mresid_martingale))
legend("topright", c("censored","event"), col=c("green","red"), lty=1)

#checking for outliers in martigale 
#martingale looks skewed [−∞, 1]
idx.outliers <- which(residuals(bestfit,"martingale") < -1.1)
data[idx.outliers,]


#plotting th deviance residuals with thickness 
plot(data$thickness, mresid_deviance, main = "lowess(deviance residuals thickness)", pch = 3, col=c("green","red"))
lines(lowess(data$thickness, mresid_deviance))
legend("topright", c("censored","event"), col=c("green","red"), lty=1)

#plotting the martingale residuals with thickness
plot(data$thickness, mresid_martingale, main = "lowess(martingale residuals thickness)", pch = 3, col=c("green","red"))
lines(lowess(data$thickness, mresid_martingale))
legend("bottomleft", c("censored","event"), col=c("green","red"), lty=1)


#check for the proportional hazard assumption

#ulcer
cox.zph(cox_mod_5)
plot(cox.zph(cox_mod_5))
abline(h=0,col="red")

#thickness
cox.zph(cox_mod_6)
plot(cox.zph(cox_mod_6))
abline(h=0,col="red")



#The proportional hazard ratio does not seem to be completely met,
#indicating that the proportional hazard changes some time, 

#conclusion
#For 1 unit additional thickness of tumour,
#the risk of death is increased by 12.1% or (1.121), 
#and Having ulcer can increase the risk of death by 238% or (3.380).



