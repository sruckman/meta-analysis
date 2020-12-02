#KEY R Review Before Holidays

#You measure height of students at king and the gym. 
#Are the heights you measured significantly different?
king <- c(126,164,148,120,178,183)
gym <- c(151,109,151,174,118,136)
t.test(king,gym)
#Paired = T if paired data
#Use x and mu = for one sample t test

#What test did you use for this question?
#two-sample t test

#What p-value was associated with this test?
#p-value = 0.3811

#What do you infer from your test? 
#You can infer that the heights measured from king and gym are 
#not significantly different


#You grow plants with three different potting soils and measure 
#height at 21 days does your data support any difference in the 
#growth with these soils.
soil1 <- c(23, 12, 45, 23, 21, 45, 21)
soil2 <- c(35, 45, 21, 34, 67, 23, 16)
soil3 <- c(16, 21, 18, 33, 16, 21, 19)
cond <- as.factor(rep(c("s1","s2","s3"),each=7))
pheight <- c(soil1,soil2,soil3)

fit <- lm(pheight~cond)
summary(aov(fit))
#p-value = 0.162
#No, since the p-value is greater than alpha, we can conclude that 
#the growths are not significantly different in the different soils

#Stickleback fish occur in deep water and shallow water populations. 
#These populations rarely interbreed. It has been hypothesized that 
#these fish have genetic adaptations to their habitat. To test this 
#you grow fish from both strains in both deep and shallow water. 
#Does the data below support the hypothesis that these fish are 
#adapted to their natural habitat? The values in the table are 
#fitnesses for fish in your experiment

dat <- as.data.frame(matrix(,24,3))
colnames(dat) <- c("fit","geno","hab")
dat$fit <- c(.97, .78, .99, .87, .91, .89,
             .61, .87,.88, .78, .80, .37,
             .56, .95,.73, .81, .89, .64,
             .77, .95,.93, .95, .89, .94)
dat$geno <- rep(c("deep","shallow"),each=12)
dat$hab <- rep(c("deep","shallow","deep","shallow"),each=6)  
conds <- rep(1:4,each=6)
fit <- lm(dat$fit~conds)
summary(aov(fit))

#P-value = 0.846
#No, since the p-value is greater than alpha, we cannot conclude
#that the fish are adapted to their natural habitats. 


############################## MCMC Data: #####################################
#Download the two mcmc log files from the course website. Choose the MCMC that represents a "good" run? 
#Provide a description of the rate parameter for codon2 and codon3.
mcmc1 <- read.csv("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis\\Learning R\\mcmc.log.csv")
mcmc2 <- read.csv("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis\\Learning R\\mcmc.log2.csv")

#par(mfcol = c(1,2))
plot(mcmc1$likelihood, type = "l")
plot(mcmc2$likelihood, type = "l")

#MCMC 1 shows a good run
par(mfcol = c(1,1))
plot(density(mcmc1$codon2[5001:10000]), xlim=c(-2,14), ylim=c(0,2))
polygon(density(mcmc1$codon2[5001:10000]), col=rgb(1,0,0,0.2))
polygon(density(mcmc1$codon3[5001:10000]), col=rgb(0,0,1,0.2))

mean(mcmc1$codon2[5001:10000])
mean(mcmc1$codon3[5001:10000])

#The rate parameters are around 0.5281994 for codon 2 and is pretty concentrated. Codon 3 is more 
#spread out and is centered around 9.44909.
