#Test Meta-Analysis for Lab Meeting R Code


###################### Effect Size Calculator ##########################
#Final results will be in terms of r^2 (effect size correlation)

metadat <- read.csv("C:\\Users\\sarah\\Documents\\test meta data.csv", header=T)

#Add and rename column to hold effect sizes
metadat <- cbind(metadat,rep(0,length(metadat$Study)),rep(0,length(metadat$Study)),rep(0,length(metadat$Study)))
names(metadat)[30] <- "r2"
names(metadat)[31] <- "Fisher_Z"
names(metadat)[32] <- "Weight"

#Calculate r effect size
for (i in 1:length(metadat$Study)){
  if (metadat$Stat.Test[i] == "t"){
    metadat$r2[i] <- sqrt((metadat$Test.Statistic[i]^2)/(metadat$Test.Statistic[i]^2 + metadat$df1[i]))
  } else if (metadat$Stat.Test[i] == "mean"){
    sdp <- sqrt((metadat$sd1[i]^2+metadat$sd2[i]^2)/2)
    d <- (metadat$mean1[i]-metadat$mean2[i])/sdp
    metadat$r2[i] <- d/sqrt(d^2+4)
  } else if (metadat$Stat.Test[i] == "X2"){
    metadat$r2[i] <- sqrt(metadat$Test.Statistic[i]/metadat$Sample.Size[i])
  } else if (metadat$Stat.Test[i] == "F"){
    SumSq <- metadat$df1[i]*metadat$Test.Statistic[i]/(metadat$df2[i])
    metadat$r2[i] <- SumSq/(1+SumSq)
  } else {
    print[i]
  }
}

#Calculate Fisher Z and weights
for (i in 1:length(metadat$Study)){
  metadat$Fisher_Z[i] <- 0.5*(log((1+metadat$r2[i])/(1-metadat$r2[i])))
  metadat$Weight[i] <- metadat$Sample.Size[i]-3
}

#Check for any missing values in r2, Fisher Z, and Weights
missingrow <- rep(0,length(metadat$Study))
for (i in 1:length(metadat$Study)){
  if (metadat$r2 == 0){
    missingrow[i] <- i
  } else {
    print("No missing value")
  } 
}
missingrow <- subset(missingrow, missingrow != 0)
missingrow


############################ Mixed Effects Model #########################
library(nlme)
fitr2 <- lme(r2 ~ 1, random = list(~1|Study, ~1|Species, ~1|Pattern, ~1|Age, 
                                 ~1|Sex, ~1|Location, ~1|Season), data=metadat)
summary(fitr2)

#Fixed effects:  r2 ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept) 0.1638881 0.05271579 48 3.108899  0.0032


#QREML
sum(fitr2$residuals^2)

#k=48, df = k-1 = 47, X^2 with df = 47 = 64.001
#Since 18.2015 is much less than 64.001, we fail to reject our null
#There is no significant heterogeneity

fitZ <- lme(Fisher_Z ~ 1, random = list(~1|Study, ~1|Species, ~1|Pattern, ~1|Age, 
                                 ~1|Sex, ~1|Location, ~1|Season), 
                                  weights = varFixed(~Weight), data=metadat)
summary(fitZ)

#Variance function:
#  Structure: fixed weights
#Formula: ~Weight 
#Fixed effects:  Fisher_Z ~ 1 
#                Value Std.Error DF  t-value p-value
#(Intercept) 0.173721 0.06669966 48 2.604526  0.0122

#QREML
sum(fitZ$residuals^2)

#k=48, df = k-1 = 47, X^2 with df = 47 = 64.001
#Since 28.6743 is much less than 64.001, we fail to reject our null
#There is no significant heterogeneity

rmel <- subset(metadat, metadat$Classification == "melanocortin")
fitrmel <- lme(r2 ~ 1, 
               random = list(~1|Study, ~1|Species, ~1|Age, ~1|Sex, ~1|Location, 
                             ~1|Season), data=rmel)
summary(fitrmel)

#Fixed effects:  r2 ~ 1 
#                Value Std.Error DF  t-value p-value
#(Intercept) 0.1825197 0.0966214 17 1.889019  0.0761

rcarot <- subset(metadat, metadat$Classification == "carotenoid")
fitrcarot <- lme(r2 ~ 1, 
               random = list(~1|Study, ~1|Species, ~1|Age, ~1|Sex, ~1|Location, 
                             ~1|Season), data=rcarot)
summary(fitrcarot)

#Fixed effects:  r2 ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept) 0.1693981 0.05644887 16 3.000913  0.0085

fitZmel <- lme(Fisher_Z ~ 1, weights = varFixed(~Weight), 
               random = list(~1|Study, ~1|Species, ~1|Age, ~1|Sex, ~1|Location, 
                             ~1|Season), data=rmel)
summary(fitZmel)

#Fixed effects:  Fisher_Z ~ 1 
#                Value Std.Error DF  t-value p-value
#(Intercept) 0.1570945 0.1149621 17 1.366489  0.1896

fitZcarot <- lme(Fisher_Z ~ 1, weights = varFixed(~Weight),
                 random = list(~1|Study, ~1|Species, ~1|Age, ~1|Sex, ~1|Location, 
                               ~1|Season), data=rcarot)
summary(fitZcarot)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept) 0.1879457 0.06683081 16 2.812261  0.0125

########################## Plots ##############################
#Effect Size vs Sample Size
plot(metadat$r2~metadat$Sample.Size, pch=16, xlab = "Sample Size",
     ylab = "Effect Size", 
     col=ifelse(metadat$Classification == "carotenoid","orange", "black"))
abline(h=0)
abline(h=mean(metadat$r2), lty=3)

#Export 4x4

#Confidence Intervals for Fisher Z and r plot
Zmean <- mean(metadat$Fisher_Z)
Zse <- 1/sqrt(sum(metadat$Sample.Size)-3)
uciZ <- Zmean + 1.96*Zse 
lciZ <- Zmean - 1.96*Zse
print(c(lciZ,uciZ))

library(DescTools)
ucir <- FisherZInv(uciZ)
lcir <- FisherZInv(lciZ)
rmean <- FisherZInv(Zmean)
print(c(lcir,ucir))

plot(NA,xlim=c(-0.5,0.5),ylim=c(-1,1),axes=F,ann=F)
axis(1)
#r Effect Size
segments(lcir,-0.5,ucir,-0.5);
points(rmean,-0.5,pch=16,xpd=NA)
#Fisher Z
segments(lciZ,0.5,uciZ,0.5);
points(Zmean,0.5,pch=16,xpd=NA)
#Add dashed line at 0
abline(v = 0, lty =3)
#Add axis labels
title(xlab = "Effect Size")
text(-0.4,-0.5,"r")
text(-0.4,0.5,"Fisher Z")

#Export 4x4

#Confidence Intervals for Fisher Z and r plot Melanocortin Data
Zmean <- mean(rmel$Fisher_Z)
Zse <- 1/sqrt(sum(rmel$Sample.Size)-3)
uciZ <- Zmean + 1.96*Zse 
lciZ <- Zmean - 1.96*Zse
print(c(lciZ,uciZ))

ucir <- FisherZInv(uciZ)
lcir <- FisherZInv(lciZ)
rmean <- FisherZInv(Zmean)
print(c(lcir,ucir))


#Confidence Intervals for Fisher Z and r plot Carotenoid Data
Zmeanc <- mean(rcarot$Fisher_Z)
Zsec <- 1/sqrt(sum(rcarot$Sample.Size)-3)
uciZc <- Zmeanc + 1.96*Zsec 
lciZc <- Zmeanc - 1.96*Zsec
print(c(lciZc,uciZc))

ucirc <- FisherZInv(uciZc)
lcirc <- FisherZInv(lciZc)
rmeanc <- FisherZInv(Zmeanc)
print(c(lcirc,ucirc))

plot(NA,xlim=c(-0.5,0.5),ylim=c(-1,1),axes=F,ann=F)
axis(1)
#r Effect Size
segments(lcir,-0.4,ucir,-0.4);
points(rmean,-0.4,pch=16,xpd=NA)
segments(lcirc,-0.6,ucirc,-0.6);
points(rmeanc,-0.6,pch=16,col = "orange",xpd=NA)
#Fisher Z
segments(lciZ,0.6,uciZ,0.6);
points(Zmean,0.6,pch=16,xpd=NA)
segments(lciZc,0.4,uciZc,0.4);
points(Zmeanc,0.4,pch=16,col = "orange",xpd=NA)
#Add dashed line at 0
abline(v = 0, lty =3)
#Add axis labels
title(xlab = "Effect Size")
text(-0.4,-0.5,"r")
text(-0.4,0.5,"Fisher Z")

#Export 4x4

#Funnel Plot
Zmean <- mean(metadat$Fisher_Z)
Zse <- 1/sqrt(sum(metadat$Sample.Size)-3)
uciZ <- Zmean + 1.96*Zse 
lciZ <- Zmean - 1.96*Zse
print(c(lciZ,uciZ))


plot(metadat$Fisher_Z, metadat$Weight, pch = 16, xlab = "Fisher Z",
     ylab = "Standard Deviation", xlim=c(-2,2),
     col=ifelse(metadat$Classification == "carotenoid","orange", "black"))
abline(v=Zmean, lty = 3)
segments(Zmean,245,Zse*sqrt(sum(metadat$Sample.Size))+Zmean,-1)
segments(Zmean,245,-1*Zse*sqrt(sum(metadat$Sample.Size))+Zmean,-1)

#Export 4x4

############################ Phylogenetic MCMC ##############################
library(MCMCglmm)
library(dplyr)

phylo <- ape::read.tree("C:\\Users\\sarah\\Documents\\species.nwk")
plot(phylo)

specieslist<-phylo$tip
metacopy <- data.frame(metadat$Species,metadat$r2)
r2 <- c()
for (i in 1:length(metacopy)){
  if (metacopy$metadat.Species[i] == specieslist[1]){
    r2[i] <- metacopy$metadat.r2
  }
}

testmcmc <- MCMCglmm(metadat$r2 ~ 1, random = ~ metadat$Species,
                     data=metadat)
