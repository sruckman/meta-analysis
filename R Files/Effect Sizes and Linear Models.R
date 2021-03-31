#Meta-Analysis Effect Size Calculations and Linear Models
setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis")

library(nlme)
library(DescTools)

####################### Effect Size Calculations ########################
#Effect Sizes will be in rho and Fisher Z values
metadat <- read.csv("Excel Sheets/meta data.csv", header=T)

#Add and rename column to hold effect sizes
metadat <- cbind(metadat,rep(0,length(metadat$Study)),rep(0,length(metadat$Study)),rep(0,length(metadat$Study)))
names(metadat)[32] <- "rho"
names(metadat)[33] <- "Fisher_Z"
names(metadat)[34] <- "SE"

#Calculate rho effect size
for (i in 1:length(metadat$Study)){
  if (metadat$Stat.Test[i] == "t"){
    metadat$rho[i] <- sqrt((metadat$Test.Statistic[i]^2)/(metadat$Test.Statistic[i]^2 + metadat$df1[i]))
  } else if (metadat$Stat.Test[i] == "mean"){
    sdp <- sqrt((metadat$sd1[i]^2+metadat$sd2[i]^2)/2)
    d <- (metadat$mean1[i]-metadat$mean2[i])/sdp
    metadat$rho[i] <- d/sqrt(d^2+4)
  } else if (metadat$Stat.Test[i] == "X2"){
    metadat$rho[i] <- sqrt(metadat$Test.Statistic[i]/metadat$Sample.Size[i])
  } else if (metadat$Stat.Test[i] == "F"){
    SumSq <- metadat$df1[i]*metadat$Test.Statistic[i]/(metadat$df2[i])
    metadat$rho[i] <- SumSq/(1+SumSq)
  } else if (metadat$Stat.Test[i] == "r"){
    metadat$rho[i] <- metadat$r[i]
  } else {
    print[i]
  }
}

#Calculate Fisher Z and SE
for (i in 1:length(metadat$Study)){
  if (metadat$rho[i] < 1){
  metadat$Fisher_Z[i] <- FisherZ(metadat$rho[i])
  metadat$SE[i] <- 1/sqrt(metadat$Sample.Size[i]-3)
  } else {
    metadat$Fisher_Z[i] <- FisherZ(0.99)
    metadat$SE[i] <- 1/sqrt(metadat$Sample.Size[i]-3)
  }
}

#Check for any missing values in r2, Fisher Z, and SE
missingrow <- rep(0,length(metadat$Study))
for (i in 1:length(metadat$Study)){
  if (metadat$Fisher_Z == "NaN"){
    missingrow[i] <- i
  } else {
    print("No missing value")
  } 
}
missingrow <- subset(missingrow, missingrow != 0)
missingrow

####################### Mixed Effects Model #########################
#Fit our entire model: Fisher Z values

fitZ <- lme(Fisher_Z ~ 1, random = list(~1|Study, ~1|Species, ~1|Pattern, ~1|Age, 
                                        ~1|Sex, ~1|Location, ~1|Season), 
            data=metadat)
summary(fitZ)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error  DF  t-value p-value
#(Intercept) 0.2918458 0.04649168 145 6.277377       0

#QREML
sum(fitZ$residuals^2)

#k=148, df = k-1 = 147, X^2 with df = 147 = 176.294
#Since 141.642 is much less than 176.294, we fail to reject our null
#There is no significant heterogeneity

#Subset the data based on class of pigments:

#Melanocortin:
rmel <- subset(metadat, metadat$Classification == "melanocortin")

fitZmel <- lme(Fisher_Z ~ 1,  
               random = list(~1|Study, ~1|Species, ~1|Age, ~1|Sex, ~1|Location, 
                             ~1|Season), data=rmel)
summary(fitZmel)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  0.245422 0.04373504 68 5.611563       0

#QREML
sum(fitZmel$residuals^2)

#k=70, df = k-1 = 69, X^2 with df = 69 = 89.391
#Since 31.17634 is much less than 89.391, we fail to reject our null
#There is no significant heterogeneity

#Carotenoid:
rcarot <- subset(metadat, metadat$Classification == "carotenoid")

fitZcarot <- lme(Fisher_Z ~ 1, 
                 random = list(~1|Study, ~1|Species, ~1|Age, ~1|Sex, ~1|Location, 
                               ~1|Season), data=rcarot)
summary(fitZcarot)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept) 0.09518496 0.05302605 57 1.79506  0.0779

#QREML
sum(fitZcarot$residuals^2)

#k=57, df = k-1 = 56, X^2 with df = 56 = 74.468
#Since 23.09871 is much less than 74.468, we fail to reject our null
#There is no significant heterogeneity

#Unknown:
runknown <- subset(metadat, metadat$Classification == "unknown")

fitZunknown <- lme(Fisher_Z ~ 1, 
                 random = list(~1|Study, ~1|Species, ~1|Age, ~1|Sex, ~1|Location, 
                               ~1|Season), data=runknown)
summary(fitZunknown)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  0.475266  0.1159627 10 4.098439  0.0021

#QREML
sum(fitZunknown$residuals^2)

#k=10, df = k-1 = 9, X^2 with df = 9 = 16.919
#Since 3.114768 is much less than 16.919, we fail to reject our null
#There is no significant heterogeneity

#None:
rnone <- subset(metadat, metadat$Classification == "none")

fitZnone <- lme(Fisher_Z ~ 1, 
                   random = list(~1|Study, ~1|Species, ~1|Age, ~1|Sex, ~1|Location, 
                                 ~1|Season), data=rnone)
summary(fitZnone)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  1.539107  0.3104424 10 4.957786   6e-04

#QREML
sum(fitZnone$residuals^2)

#k=10, df = k-1 = 9, X^2 with df = 9 = 16.919
#Since 22.32293 is greater than 16.919, we reject our null
#There is significant heterogeneity

####################### Confidence Intervals #######################
#Entire Model:
Zmean <- 0.2918458 
uci <- 0.2918458 + 1.96*0.04649168
lci <- 0.2918458 - 1.96*0.04649168
print(c(lci,uci))

#Melanocortin:
Zmel <- 0.245422 
ucimel <- 0.245422 + 1.96*0.04373504
lcimel <- 0.245422 - 1.96*0.04373504
print(c(lcimel,ucimel))

#Carotenoid:
Zcar <- 0.09518496 
ucicar <- 0.09518496 + 1.96*0.05302605
lcicar <- 0.09518496 - 1.96*0.05302605
print(c(lcicar,ucicar))

#Unknown:
Zun <- 0.475266   
uciun <- 0.475266 + 1.96*0.1159627
lciun <- 0.475266 - 1.96*0.1159627
print(c(lciun,uciun))

#None:
Znone <- 1.539107   
ucinone <- 1.539107 + 1.96*0.3104424
lcinone <- 1.539107 - 1.96*0.3104424
print(c(lcinone,ucinone))

####################### Plots #######################
#Funnel Plot:
Zmean <- mean(metadat$Fisher_Z)
Zse <- mean(metadat$SE)
uciZ <- Zmean + 1.96*Zse 
lciZ <- Zmean - 1.96*Zse
print(c(lciZ,uciZ))

cols <- rep(0,length(metadat$Classification))
for (i in 1:length(metadat$Classification)){
  if (metadat$Classification[i] == "carotenoid"){
    cols[i] <- "orange"
  } else if (metadat$Classification[i] == "melanocortin"){
    cols[i] <- "black"
  } else if (metadat$Classification[i] == "unknown"){
    cols[i] <- "purple"
  } else if (metadat$Classification[i] == "none"){
    cols[i] <- "blue"
  } 
}

plot(metadat$Fisher_Z, metadat$SE, pch = 16, xlab = "Fisher Z",
     ylab = "Standard Error", xlim=c(-3,3),
     col=cols)
abline(v=Zmean, lty = 3)
segments(Zmean,1.05,uciZ,0)
segments(Zmean,1.05,lciZ,0)
abline(v=0,col="red")

legend(x = -3.15, y = 1.02, legend = c("Carotenoid", "Melanocortin", "Unknown", "Eye"),
       fill = c("orange", "black", "purple", "blue"), cex = 0.7,
       border = "white", box.col = "white")

#Export 4x4.5

#Fisher z plot for each model:
plot(NA,xlim=c(-1,2.5),ylim=c(-1,1),axes=F,ann=F)
axis(1)
#Fisher Z
#Melanin
segments(lcimel,0.8,ucimel,0.8);
points(Zmel,0.8,pch=16,xpd=NA)
#Carotenoid
segments(lcicar,0.6,ucicar,0.6);
points(Zcar,0.6,pch=16,col = "orange",xpd=NA)
#Unknown
segments(lciun,0.4,uciun,0.4);
points(Zun,0.4,pch=16,col = "purple",xpd=NA)
#None
segments(lcinone,0.2,ucinone,0.2);
points(Znone,0.2,pch=16,col = "blue",xpd=NA)
#Overall
segments(lci,-0.5,uci,-0.5);
points(Zmean,-0.5,pch=16,xpd=NA)
#Add dashed line at 0
abline(v = 0, lty =3)
#Add axis labels
title(xlab = "Fisher Z")
text(-1,-0.5,"Overall", cex = 0.8, adj = c(0,0))
text(-1,0.8,"Melanocortin", cex = 0.8, adj = c(0,0))
text(-1,0.6,"Carotenoid", cex = 0.8, adj = c(0,0))
text(-1,0.4,"Unknown", cex = 0.8, adj = c(0,0))
text(-1,0.2,"Eye", cex = 0.8, adj = c(0,0))

#Export 5x5
