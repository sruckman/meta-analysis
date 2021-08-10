#Meta-Analysis Effect Size Calculations and Linear Models
setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis")

library(nlme)
library(DescTools)
library(scales)
library(metafor)

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
                          ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity, 
                          ~1|Classification), 
            weights = varFixed(~Weight), data=metadat)
summary(fitZ)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error  DF  t-value p-value
#(Intercept) 0.3059841 0.05830351  80 5.248125       0

#QREML for Heterogeneity Test
sum(fitZ$residuals^2)

#k=149, df = k-1 = 148, X^2 with df = 148 = 177.390
#Since 196.262 is more than 177.390, we reject our null
#There is significant heterogeneity
 
#Publication Bias (Egger's Test):
ranktest(x = metadat$Fisher_Z, sei = metadat$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.1824, p = 0.0011

#There is a significant publication bias. 

#### Only Using Known Color Genes ####
known <- subset(metadat, metadat$Classification != "unknown")

fitkn <- lme(Fisher_Z ~ 1, 
            random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age, 
                          ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity, 
                          ~1|Classification), 
            weights = varFixed(~Weight), data=known)
summary(fitkn)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error  DF  t-value p-value
#(Intercept) 0.2365914 0.04686318 73 5.048556       0

#QREML
sum(fitkn$residuals^2)

#k=131, df = k-1 = 130, X^2 with df = 130 = 157.610
#Since 79.2127 is less than 157.610, we do not reject our null
#There is no significant heterogeneity

ranktest(x = known$Fisher_Z, sei = known$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.1044, p = 0.0801

#There is no significant publication bias

#### Subset the data based on class of pigments: ####
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

#Publication Bias (Egger's Test):
ranktest(x = rmel$Fisher_Z, sei = rmel$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.1834, p = 0.0267

#There is a significant publication bias.

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

#Publication Bias (Egger's Test):
ranktest(x = reu$Fisher_Z, sei = reu$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.2859, p = 0.0115

#There is a significant publication bias.

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

#Publication Bias (Egger's Test):
ranktest(x = rph$Fisher_Z, sei = rph$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = -0.0550, p = 0.6695

#There is no significant publication bias.

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

#Publication Bias (Egger's Test):
ranktest(x = rcarot$Fisher_Z, sei = rcarot$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.0380, p = 0.6750

#There is no significant publication bias.

#### Unknown: ####
runknown <- subset(metadat, metadat$Classification == "unknown")

fitZunknown <- lme(Fisher_Z ~ 1, 
                   random = list(~1|Authors, ~1|Species, ~1|Pattern, ~1|Age,
                                 ~1|Sex, ~1|Location, ~1|Season, ~1|Plasticity), 
                   weights = varFixed(~Weight), data=runknown)
summary(fitZunknown)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  0.9531082  0.316791 11 3.008634  0.0119

#QREML
sum(fitZunknown$residuals^2)

#k=18, df = k-1 = 17, X^2 with df = 17 = 27.587
#Since 100.8544 is much greater than 27.587, we reject our null
#There is significant heterogeneity

#Publication Bias (Egger's Test):
ranktest(x = runknown$Fisher_Z, sei = runknown$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.4303, p = 0.0155

#There is a significant publication bias.

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

#Publication Bias (Egger's Test):
ranktest(x = rvert$Fisher_Z, sei = rvert$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.1698, p = 0.0031

#There is a significant publication bias.

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
#Since 2.234939 is less than 14.067, we fail to reject our null
#There is no significant heterogeneity

#Publication Bias (Egger's Test):
ranktest(x = rin$Fisher_Z, sei = rin$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.6910, p = 0.0178

#There is a significant publication bias.

#### Plastic: ####
rpl <- subset(metadat, metadat$Plasticity == "Plastic")

fitZpl <- lme(Fisher_Z ~ 1, 
              random = list(~1|Authors, ~1|Species, ~1|Pattern, 
                            ~1|Sex, ~1|Location, ~1|Plasticity), 
              weights = varFixed(~Weight), data=rpl)
summary(fitZpl)

#Fixed effects:  Fisher_Z ~ 1 
#                Value  Std.Error DF  t-value p-value
#(Intercept)  0.3703848 0.1340427 33 2.763185  0.0093

#QREML
sum(fitZpl$residuals^2)

#k=49, df = k-1 = 48, X^2 with df = 48 = 65.171
#Since 107.8682 is greater than 65.171, we reject our null
#There is significant heterogeneity

#Publication Bias (Egger's Test):
ranktest(x = rpl$Fisher_Z, sei = rpl$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.2443, p = 0.0142

#There is a significant publication bias.

#### Non-Plastic: ####
rnpl <- subset(metadat, metadat$Plasticity == "No")

fitZnpl <- lme(Fisher_Z ~ 1, 
              random = list(~1|Authors, ~1|Species, ~1|Pattern, 
                            ~1|Sex, ~1|Location, ~1|Plasticity), 
              weights = varFixed(~Weight), data=rnpl)
summary(fitZnpl)

#Fixed effects:  Fisher_Z ~ 1 
#                 Value  Std.Error DF  t-value p-value
#(Intercept)  0.2704223 0.04637794 48 5.830839       0

#QREML
sum(fitZnpl$residuals^2)

#k=100, df = k-1 = 99, X^2 with df = 99 = 123.225	
#Since 46.37348 is less than 123.225, we fail to reject our null
#There is no significant heterogeneity

#Publication Bias (Egger's Test):
ranktest(x = rnpl$Fisher_Z, sei = rnpl$SE)

#Rank Correlation Test for Funnel Plot Asymmetry
#Kendall's tau = 0.1132, p = 0.0997

#There is no significant publication bias.
####################### Confidence Intervals #######################
#Entire Model:
Zmean <- 0.3059841   
uci <- Zmean + 1.96*0.05830351
lci <- Zmean - 1.96*0.05830351
print(c(lci,uci))

#Only Known Color Genetics:
Zknown <- 0.2365914 
ucikn <- Zknown + 1.96*0.04686318
lcikn <- Zknown - 1.96*0.04686318
print(c(lcikn,ucikn))

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
Zun <- 0.9531082     
uciun <- Zun + 1.96*0.316791
lciun <- Zun - 1.96*0.316791
print(c(lciun,uciun))

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

#Plastic
Zpl <- 0.3703848 
ucipl <- Zpl + 1.96*0.1340427
lcipl <- Zpl - 1.96*0.1340427
print(c(lcipl, ucipl))

#Non-Plastic
Znpl <- 0.2704223
ucinpl <- Znpl + 1.96*0.04637794
lcinpl <- Znpl - 1.96*0.04637794
print(c(lcinpl, ucinpl))

####################### Plots #######################
#### Funnel Plot: ####
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
  } 
}

funnel(x = metadat$Fisher_Z, sei = metadat$SE, yaxis = "sei", 
       xlab = "Fisher Z", col = alpha(cols, 0.75), back = "white",
       xlim = c(-3.5,3), pch=shapes)
legend(x = -3.73, y = -0.05, 
       legend = c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown"),
       fill = c("orange", "black","orangered3", "darkorchid4"), cex = 0.7,
       border = "white", box.col = "white", )

#Export 5x5

#Funnel Plot for Plasticity and Verts/Inverts
colsv <- rep(0,length(metadat$Vert_Invert))
for (i in 1:length(metadat$Vert_Invert)){
  if (metadat$Vert_Invert[i] == "vertebrate"){
    colsv[i] <- "dodgerblue2"
  } else if (metadat$Vert_Invert[i] == "invertebrate"){
    colsv[i] <- "deeppink"
  } 
}

shapes <- rep(0,length(metadat$Plasticity))
for (i in 1:length(metadat$Vert_Invert)){
  if (metadat$Plasticity[i] == "Plastic"){
    shapes[i] <- 20
  } else if (metadat$Plasticity[i] == "No"){
    shapes[i] <- 18
  } 
}


funnel(x = metadat$Fisher_Z, sei = metadat$SE, yaxis = "sei", 
       xlab = "Fisher Z", col = alpha(colsv, 0.75), back = "white",
       xlim = c(-3.5,3), pch = shapes)
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown", "Plastic", "Non-Plastic")
Cols <- c("orange", "black","orangered3", "darkorchid4", "black", "black")
points <- c(15,15,15,15,20,18)
ys <- c(0,0.05,0.1,0.15,0.2,0.25)
for(i in 1:6){
  points(x=-3.6, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-3.6,y=ys[i], labels=words[i], pos=4,cex=.75)
}


#Export 5x5

#### Fisher z plot for each model: ####
plot(NA,xlim=c(-1,3),ylim=c(-1,1.6),axes=F,ann=F)
axis(1)
#Fisher Z
#Carotenoid
segments(lcicar,1.2,ucicar,1.2);
points(Zcar,1.2,pch=16,col = "orange",xpd=NA)
#Eumelanin
segments(lcieu,1,ucieu,1);
points(Zeu,1,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lciph,0.8,uciph,0.8);
points(Zph,0.8,pch=16,col = "orangered3",xpd=NA)
#Unknown
segments(lciun,0.6,uciun,0.6);
points(Zun,0.6,pch=16,col = "darkorchid4",xpd=NA)
#Vertebrates
segments(lcivert,0.4,ucivert,0.4);
points(Zvert,0.4,pch=16,col = "dodgerblue4",xpd=NA)
#Invertebrates
segments(lciin,0.2,uciin,0.2);
points(Zin,0.2,pch=16,col = "dodgerblue2",xpd=NA)
#Plastic
segments(lcipl,0,ucipl,0);
points(Zpl,0,pch=16,col = "darkgreen",xpd=NA)
#Non-Plastic
segments(lcinpl,-0.2,ucinpl,-0.2);
points(Znpl,-0.2,pch=16,col = "mediumseagreen",xpd=NA)


#Overall
segments(lci,-0.6,uci,-0.6);
points(Zmean,-0.6,pch=16,xpd=NA)
#Color Genes Only
segments(lciin,-0.8,uciin,-0.8);
points(Zin,-0.8,pch=16,xpd=NA)
#Melanin
segments(lcimel,-1,ucimel,-1);
points(Zmel,-1,pch=16,xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-1,-0.6,"Overall", cex = 0.9, adj = c(0,0))
text(-1,-0.8,"Only Known", cex = 0.9, adj = c(0,0))
text(-1,-1,"Melanocortin", cex = 0.9, adj = c(0,0))
text(-1,1.2,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-1,1,"Eumelanin", cex = 0.9, adj = c(0,0))
text(-1,0.8,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-1,0.6,"Unknown", cex = 0.9, adj = c(0,0))
text(-1,0.4,"Vertebrates", cex = 0.9, adj = c(0,0))
text(-1,0.2,"Invertebrates", cex = 0.9, adj = c(0,0))
text(-1,0,"Plastic", cex = 0.9, adj = c(0,0))
text(-1,-0.2,"Non-Plastic", cex = 0.9, adj = c(0,0))

#Export 6x6
