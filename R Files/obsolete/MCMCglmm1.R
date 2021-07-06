#MCMCglmm
setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis")

library(ape)
library(MCMCglmm)
library(phytools)
library(TreeTools)

#### Phylogeny ####
tree <- read.tree("Excel Sheets/list.nwk")
root(tree, outgroup = "Phymactis_clematis")
plotTree(tree, fsize = 0.7, lwd =1, ftype = "i")

#Export 6x6
#To clear the tree:
dev.off()


#### Data ####
metadata <- read.csv("Excel Sheets/metafull.csv", header=T)
names(metadata)[5] <- "animal"
#Check to see if the labels match
tree$tip.label<-gsub("_"," ",tree$tip.label)
tips <- TipLabels(tree)
tips %in% metadata$animal

################# MCMC ##############
#for MCMCglmm to run we need to change our edge lengths from 0 to 0.000001
for (i in 1:length(tree$edge.length)){
  if (tree$edge.length[i] == 0){
    tree$edge.length[i] <- 0.000001
  } else {
    print(i)
  }
}
print(tree$edge.length)

#priors:
prior.ex<- list(G = list(G1 = list(V = 1, nu = 0.02), 
                         G2 = list(V = 1, nu = 0.02), 
                         G3 = list(V = 1, nu = 0.02),
                         G4 = list(V = 1, nu = 0.02), 
                         G5 = list(V = 1, nu = 0.02), 
                         G6 = list(V = 1, nu = 0.02), 
                         G7 = list(V = 1, nu = 0.02), 
                         G8 = list(V = 1, nu = 0.02),
                         G9 = list(V = 1, nu = 0.02),
                         G10 = list(V = 1, nu = 0.02)), 
                R = list(V=1, nu=0.02))

#Run the model:
full.mcmc <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + Pattern + 
                        Age + Sex + Location + Season + Plasticity +
                        Classification + us(Weight):units, data=metadata, pedigree = tree, 
                      nitt = 200000, thin = 50, burnin = 190000, 
                      prior = prior.ex)

summary(full.mcmc)

plot(full.mcmc$Sol)

# Proportion of variance explained by random factors
rand <- full.mcmc$VCV/apply(full.mcmc$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#Units = Effect Size Specific Effects
