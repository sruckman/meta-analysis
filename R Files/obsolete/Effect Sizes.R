#Effect Size Calculator 
#Final results will be in terms of r^2 (effect size correlation)

metadat <- read.csv("C:\\Users\\sarah\\Documents\\test meta data.csv", header=T)
r2 <- rep(NA,length(metadat$Study))
metadat <- cbind(metadat,r2)

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


#Means and Standard Deviations
msd <- function(m1,sd1,m2,sd2){
  sdp <- sqrt((sd1^2+sd^2)/2)
  d <- (m1-m2)/sdp
  r <- d/sqrt(d^2+4)
  print(r)
}

#T values and df
tdf <- function(t,df){
  r <- sqrt((t^2)/(t^2+df))
  print(r)
}

#Wald and sample sizes
wald <- function(w,n1,n2){
  d <- w*sqrt((1/n1)+1/n2)
  print(d)
}

#Chi-square with df = 1
x1 <- function(x,n){
  r <- sqrt(x/n)
  print(r)
}

#One way ANOVA
oneano <- function(f,n1,n2){
  r <- sqrt(f)/sqrt(f+n1+n2-2)
  print(r)
}

#glmm 
glmm <- function(FValue,NumDF,DenDF){
  SumSq = NumDF*FValue/(DenDF)
  r = SumSq/(1+SumSq)
}