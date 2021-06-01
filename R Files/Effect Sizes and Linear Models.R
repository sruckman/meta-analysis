#Meta-Analysis Effect Size Calculations and Linear Models
setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis")

library(nlme)
library(DescTools)
library(scales)

####################### Effect Size Calculations ########################
#Effect Sizes will be in rho and Fisher Z values
metadat <- read.csv("Excel Sheets/meta data.csv", header=T)

#Add and rename column to hold effect sizes
metadat <- cbind(metadat,rep(0,length(metadat$Study)),rep(0,length(metadat$Study)),
                 rep(0,length(metadat$Study)),rep(0,length(metadat$Study)))
names(metadat)[35] <- "rho"
names(metadat)[36] <- "Fisher_Z"
names(metadat)[37] <- "SE"
names(metadat)[38] <- "Weight"

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

#Calculate Fisher Z, SE, and Weight
for (i in 1:length(metadat$Study)){
  if (metadat$rho[i] < 1){
  metadat$Fisher_Z[i] <- FisherZ(metadat$rho[i])
  metadat$SE[i] <- 1/sqrt(metadat$Sample.Size[i]-3)
  metadat$Weight[i] <- sqrt(metadat$Sample.Size[i]-3)
  } else {
    metadat$Fisher_Z[i] <- FisherZ(0.99)
    metadat$SE[i] <- 1/sqrt(metadat$Sample.Size[i]-3)
    metadat$Weight[i] <- sqrt(metadat$Sample.Size[i]-3)
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

write.csv(metadat, "metafull.csv")
####################### Mixed Effects Model #########################
### Fit our entire model: Fisher Z values ###

fitZ <- lme(Fisher_Z ~ 1, 
            random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age, 
                          ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity), 
            weights = varFixed(~Weight), data=metadat)
summary(fitZ)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error  DF  t-value p-value
#(Intercept) 0.3059841 0.05830351  80 5.248125       0

#QREML
sum(fitZ$residuals^2)

#k=149, df = k-1 = 148, X^2 with df = 148 = 177.390
#Since 179.6119 is more than 177.390, we reject our null
#There is significant heterogeneity

### Subset the data based on class of pigments: ###

#### Melanocortin: ####
rmel <- subset(metadat, metadat$Classification == "melanocortin")

fitZmel <- lme(Fisher_Z ~ 1,  
               random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age,
                             ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity), 
               weights = varFixed(~Weight), data=rmel)
summary(fitZmel)

#Fixed effects:  Fisher_Z ~ 1 
#                 Value Std.Error DF  t-value p-value
#(Intercept)  0.2851931 0.0678219 40 4.205029   1e-04

#QREML
sum(fitZmel$residuals^2)

#k=71, df = k-1 = 70, X^2 with df = 70 = 90.531
#Since 27.273 is much less than 90.531, we fail to reject our null
#There is no significant heterogeneity

#### Eumelanin: ####
reu <- subset(metadat, metadat$Eu_Pheomelanin == "eumelanin")

fitZeu <- lme(Fisher_Z ~ 1,  
              random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age,
                            ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity), 
              weights = varFixed(~Weight), data=reu)
summary(fitZeu)

#Fixed effects:  Fisher_Z ~ 1 
#                Value   Std.Error DF  t-value p-value
#(Intercept)  0.2805119 0.08520859 29 3.292061  0.0026

#QREML
sum(fitZeu$residuals^2)

#k=40, df = k-1 = 39, X^2 with df = 39 = 54.572
#Since 15.19769 is much less than 54.572, we fail to reject our null
#There is no significant heterogeneity

#### Pheomelanin: ####
rph <- subset(metadat, metadat$Eu_Pheomelanin == "pheomelanin")

fitZph <- lme(Fisher_Z ~ 1,  
              random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age,
                            ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity), 
              weights = varFixed(~Weight), data=rph)
summary(fitZph)

#Fixed effects:  Fisher_Z ~ 1 
#                Value   Std.Error DF  t-value p-value
#(Intercept)  0.2786944  0.1011602 19 2.754981  0.0126

#QREML
sum(fitZph$residuals^2)

#k=32, df = k-1 = 31, X^2 with df = 31 = 44.985
#Since 13.37308 is much less than 44.985, we fail to reject our null
#There is no significant heterogeneity

#### Carotenoid: ####
rcarot <- subset(metadat, metadat$Classification == "carotenoid")

fitZcarot <- lme(Fisher_Z ~ 1, 
                 random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age,
                               ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity), 
                 weights = varFixed(~Weight), data=rcarot)
summary(fitZcarot)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF t-value p-value
#(Intercept) 0.1611488 0.06641673 32 2.426328  0.0211

#QREML
sum(fitZcarot$residuals^2)

#k=59, df = k-1 = 58, X^2 with df = 58 = 76.778
#Since 46.50761 is much less than 76.778, we fail to reject our null
#There is no significant heterogeneity

#### Unknown: ####
runknown <- subset(metadat, metadat$Classification == "unknown")

fitZunknown <- lme(Fisher_Z ~ 1, 
                   random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age,
                                 ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity), 
                   weights = varFixed(~Weight), data=runknown)
summary(fitZunknown)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  0.5394419 0.1921933  5 2.806768  0.0377

#QREML
sum(fitZunknown$residuals^2)

#k=8, df = k-1 = 7, X^2 with df = 7 = 14.067
#Since 3.878083 is much less than 14.067, we fail to reject our null
#There is no significant heterogeneity

#### Eye: ####
rnone <- subset(metadat, metadat$Classification == "none")

fitZnone <- lme(Fisher_Z ~ 1, 
                random = list(~1|Authors, ~1|Species), 
                weights = varFixed(~Weight), data=rnone)
summary(fitZnone)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  1.574762 0.5345479  8  2.945971  0.0185

#QREML
sum(fitZnone$residuals^2)

#k=10, df = k-1 = 9, X^2 with df = 9 = 16.919
#Since 26.36726 is greater than 16.919, we reject our null
#There is significant heterogeneity

#### Vertebrates: ####
rvert <- subset(metadat, metadat$Vert_Invert == "vertebrate")

fitZvert <- lme(Fisher_Z ~ 1, 
                random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age,
                              ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity), 
                weights = varFixed(~Weight), data=rvert)
summary(fitZvert)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  0.292554 0.06176129 74 4.736851       0

#QREML
sum(fitZvert$residuals^2)

#k=141, df = k-1 = 140, X^2 with df = 140 = 168.613
#Since 176.9077 is greater than 168.613, we reject our null
#There is significant heterogeneity

#### Invertebrates: ####
rin <- subset(metadat, metadat$Vert_Invert == "invertebrate")

fitZin <- lme(Fisher_Z ~ 1, 
                random = list(~1|Authors, ~1|Species, ~1|Pattern, 
                              ~1|Sex, ~1|Location, ~1|Plasticity), 
                weights = varFixed(~Weight), data=rin)
summary(fitZin)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  0.4672572 0.1305277  6 3.579754  0.0116

#QREML
sum(fitZin$residuals^2)

#k=8, df = k-1 = 7, X^2 with df = 7 = 14.067
#Since 2.234939 is greater than 14.067, we fail to reject our null
#There is no significant heterogeneity


####################### Wald Tests #######################

####################### Confidence Intervals #######################
#Entire Model:
Zmean <- 0.3059841   
uci <- Zmean + 1.96*0.05830351
lci <- Zmean - 1.96*0.05830351
print(c(lci,uci))

#Melanocortin:
Zmel <- 0.2851931  
ucimel <- Zmel + 1.96*0.0678219
lcimel <- Zmel - 1.96*0.0678219
print(c(lcimel,ucimel))

#Eumelanin:
Zeu <- 0.2805119  
ucieu <- Zeu + 1.96*0.08520859
lcieu <- Zeu - 1.96*0.08520859
print(c(lcieu,ucieu))

#Pheomelanin:
Zph <- 0.2786944  
uciph <- Zph + 1.96*0.1011602
lciph <- Zph - 1.96*0.1011602
print(c(lciph,uciph))


#Carotenoid:
Zcar <- 0.1611488  
ucicar <- Zcar + 1.96*0.06641673
lcicar <- Zcar - 1.96*0.06641673
print(c(lcicar,ucicar))

#Unknown:
Zun <- 0.5394419   
uciun <- Zun + 1.96*0.1921933
lciun <- Zun - 1.96*0.1921933
print(c(lciun,uciun))

#Eye:
Znone <- 1.574762   
ucinone <- Znone + 1.96*0.5345479
lcinone <- Znone - 1.96*0.5345479
print(c(lcinone,ucinone))

#Vertebrates
Zvert <- 0.292554  
ucivert <- Zvert + 1.96*0.06176129
lcivert <- Zvert - 1.96*0.06176129
print(c(lcivert,ucivert))

#Invertebrates
Zin <- 0.4672572    
uciin <- Zin + 1.96*0.1305277
lciin <- Zin - 1.96*0.1305277
print(c(lciin,uciin))

####################### Plots #######################
#Funnel Plot:
Zmeans <- mean(metadat$Fisher_Z)
Zse <- mean(metadat$SE)
uciZ <- Zmeans + 1.96*Zse 
lciZ <- Zmeans - 1.96*Zse
print(c(lciZ,uciZ))

cols <- rep(0,length(metadat$Classification))
for (i in 1:length(metadat$Classification)){
  if (metadat$Classification[i] == "carotenoid"){
    cols[i] <- "orange"
  } else if (metadat$Eu_Pheomelanin[i] == "eumelanin"){
    cols[i] <- "black"
  } else if (metadat$Eu_Pheomelanin[i] == "pheomelanin"){
    cols[i] <- "orangered3"
  } else if (metadat$Classification[i] == "unknown"){
    cols[i] <- "darkorchid4"
  } else if (metadat$Classification[i] == "none"){
    cols[i] <- "dodgerblue4"
  } 
}

plot(metadat$Fisher_Z, metadat$SE, pch = 20, xlab = "Fisher Z",
     ylab = "Standard Error", xlim=c(-3,3),
     col=alpha(cols,0.75))
abline(v=Zmeans, lty = 3)
segments(Zmeans,1.05,uciZ,0)
segments(Zmeans,1.05,lciZ,0)
#abline(v=0,col="red")

legend(x = -3.15, y = 1.02, 
       legend = c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown", "Eye"),
       fill = c("orange", "black","orangered3", "darkorchid4", "dodgerblue4"), cex = 0.7,
       border = "white", box.col = "white")

#Export 5x5

#Fisher z plot for each model:
plot(NA,xlim=c(-1,3),ylim=c(-1,1),axes=F,ann=F)
axis(1)
#Fisher Z
#Carotenoid
segments(lcicar,1,ucicar,1);
points(Zcar,1,pch=16,col = "orange",xpd=NA)
#Eumelanin
segments(lcieu,0.8,ucieu,0.8);
points(Zeu,0.8,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lciph,0.6,uciph,0.6);
points(Zph,0.6,pch=16,col = "orangered3",xpd=NA)
#Unknown
segments(lciun,0.4,uciun,0.4);
points(Zun,0.4,pch=16,col = "darkorchid4",xpd=NA)
#Eye
segments(lcinone,0.2,ucinone,0.2);
points(Znone,0.2,pch=16,col = "dodgerblue4",xpd=NA)
#Vertebrates
segments(lcivert,0,ucivert,0);
points(Zvert,0,pch=16,col = "darkturquoise",xpd=NA)
#Invertebrates
segments(lciin,-0.2,uciin,-0.2);
points(Zin,-0.2,pch=16,col = "darkslategray4",xpd=NA)

#Overall
segments(lci,-0.6,uci,-0.6);
points(Zmean,-0.6,pch=16,xpd=NA)
#Melanin
segments(lcimel,-0.8,ucimel,-0.8);
points(Zmel,-0.8,pch=16,xpd=NA)
#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-1,-0.6,"Overall", cex = 0.9, adj = c(0,0))
text(-1,-0.8,"Melanocortin", cex = 0.9, adj = c(0,0))
text(-1,1,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-1,0.8,"Eumelanin", cex = 0.9, adj = c(0,0))
text(-1,0.6,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-1,0.4,"Unknown", cex = 0.9, adj = c(0,0))
text(-1,0.2,"Eye", cex = 0.9, adj = c(0,0))
text(-1,0,"Vertebrates", cex = 0.9, adj = c(0,0))
text(-1,-0.2,"Invertebrates", cex = 0.9, adj = c(0,0))

#Export 6x6
