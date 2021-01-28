############################## Practice Output ##########################
#Set your working directory using setwd()
setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis\\Learning R")

mcmc1 <- read.csv("mcmc.log.csv")
mcmc2 <- read.csv("mcmc.log2.csv")

#par(mfcol = c(1,2))
plot(mcmc1$likelihood, type = "l")
plot(mcmc2$likelihood, type = "l")

#MCMC 1 shows a good run
par(mfcol = c(1,1))
plot(density(mcmc1$codon2[5001:10000]), xlim=c(-2,14), ylim=c(0,2))
polygon(density(mcmc1$codon2[5001:10000]), col=rgb(1,0,0,0.2))
polygon(density(mcmc1$codon3[5001:10000]), col=rgb(0,0,1,0.2))
polygon(density(mcmc1$codon1[5001:10000]), col=rgb(0,1,0,0.2))

mean(mcmc1$codon2[5001:10000])
mean(mcmc1$codon3[5001:10000])

#The rate parameters are around 0.5281994 for codon 2 and is pretty concentrated. Codon 3 is more 
#spread out and is centered around 9.44909.

########################### Your Turn #################################
log <- read.csv("test.log.csv")
plot(log$p, type = "l")

plot(density(log$asc1[1500:1800]), xlim=c(0,12), ylim=c(0,1), main = "")
polygon(density(log$asc1[1500:1800]), col=rgb(1,0,0,0.2))
polygon(density(log$desc1[1500:1800]), col=rgb(0,0,1,0.2))

mean(log$asc1[1500:1800])
mean(log$desc1[1500:1800])

#Fission mean = 5.457657 (red)
#Fusion mean = 2.130224 (blue)