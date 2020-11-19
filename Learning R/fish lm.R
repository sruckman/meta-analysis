
#fixed effects model
fish <- read.csv("C:\\Users\\sarah\\Documents\\Exp Design\\Scripts\\fish2.csv")
#You have different sample sizes and this means that you lose power with the smaller sample sizes
#This is wide data and we need to make it long so that each value has one column
fitness <- unlist(fish)
env <- c(rep("deep",24), rep("shallow",24))
geno <- c(rep("deep",12), rep("shallow",12), rep("deep",12),rep("shallow",12))
dat2 <- data.frame(fitness,env,geno)              
dat2 <- dat2[!is.na(dat2$fitness),]  #also dat2[complete.cases(dat2),]

library(ggplot2)
#library(ggraptR)
#ggraptR(dat2)
ggplot(dat2, aes(y=fitness, x=as.factor(env))) + 
  geom_boxplot(aes(fill=as.factor(geno)), stat="boxplot", 
               position=position_dodge(0.4), width=0.2) +
  theme_bw() + 
  theme(text=element_text(color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="geno")) + 
  xlab("Environment") + 
  ylab("Fitness")

fit <- glm(dat2$fitness ~ dat2$env*dat2$geno)
summary(fit)
#All terms are significant
#Exp for dat2$envshallow: The mean for shallow = 17.475-12.125 = 5.35
#Exp for dat2$genoshallow: The mean shallow genotype fitness mean 17-11
#Intercept = deep expectation for both environment and genotype
#Exp for interaction: how different are shallow shallow to deep deep = 17-12-11+11
#these values minimize the residuals


