#Meta-Analysis Effect Size Calculations
setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis")

library(DescTools)

#Effect Sizes will be in rho and Fisher Z values
metadat <- read.csv("Excel Sheets/meta_data.csv", header=T)

#Add and rename column to hold effect sizes calculations and standard errors
metadat <- cbind(metadat,rep(0,length(metadat$Study)),rep(0,length(metadat$Study)),
                 rep(0,length(metadat$Study)),rep(0,length(metadat$Study)),
                 rep(0,length(metadat$Study)))
names(metadat)[42] <- "rho"
names(metadat)[43] <- "SE_r"
names(metadat)[44] <- "Fisher_Z"
names(metadat)[45] <- "SE_Z"
names(metadat)[46] <- "Weight"

#Calculate rho effect size based on test statsitics given
for (i in 1:length(metadat$Study)){
  if (metadat$Stat.Test[i] == "t"){
    rpb <- sqrt((metadat$Test.Statistic[i]^2)/(metadat$Test.Statistic[i]^2 + metadat$df1[i]))
    p <- metadat$Sample.Size[i]/(metadat$Sample.Size[i]+metadat$Sample.Size[i])
    metadat$rho[i] <- ((sqrt(p*(1-p)))/dnorm(qnorm(p)))*rpb
  } else if (metadat$Stat.Test[i] == "mean"){
    m <- metadat$n1[i]+metadat$n2[i]-2
    d <- (metadat$mean2[i] - metadat$mean1[i])/(sqrt(((metadat$n1[i]-1)*metadat$sd1[i]^2)+((metadat$n2[i]-1)*metadat$sd2[i]^2)/(m)))
    h <- (m/metadat$n1[i]) + (m/metadat$n2[i])
    rpb <- d/sqrt((d^2)+h)
    p <- metadat$n1[i]/(metadat$n1[i]+metadat$n2[i])
    metadat$rho[i] <- ((sqrt(p*(1-p)))/dnorm(qnorm(p)))*rpb
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

#Calculate Fisher Z, SE for both r and Fisher Z, and Weight
#if values are greater than 1 then we will use 0.99
#if valuse are less than -1 then we will use -0.99
for (i in 1:length(metadat$Study)){
  if (metadat$rho[i] < 1 && metadat$rho[i] > -1){
    metadat$SE_r[i] <- (sqrt(1-metadat$rho[i]^2))/(sqrt(metadat$Sample.Size[i]-2))
    metadat$Fisher_Z[i] <- FisherZ(metadat$rho[i])
    metadat$SE_Z[i] <- 1/sqrt(metadat$Sample.Size[i]-3)
    metadat$Weight[i] <- sqrt(metadat$Sample.Size[i]-3)
  } else if (metadat$rho[i] > 1){
    metadat$SE_r[i] <- (sqrt(1-0.99^2))/(sqrt(metadat$Sample.Size[i]-2))
    metadat$Fisher_Z[i] <- FisherZ(0.99)
    metadat$SE_Z[i] <- 1/sqrt(metadat$Sample.Size[i]-3)
    metadat$Weight[i] <- sqrt(metadat$Sample.Size[i]-3)
  } else {
    metadat$SE_r[i] <- (sqrt(1-(-0.99)^2))/(sqrt(metadat$Sample.Size[i]-2))
    metadat$Fisher_Z[i] <- FisherZ(-0.99)
    metadat$SE_Z[i] <- 1/sqrt(metadat$Sample.Size[i]-3)
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

write.csv(metadat, "meta_complete_data2.csv")
