#MCMCglmm
setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis")

library(ape)
library(MCMCglmm)
library(TreeTools)
library(phyr)
library(phytools)

#### Phylogeny ####
tree <- read.tree("Excel Sheets/obsolete/list.nwk")
root(tree, outgroup = "Phymactis_clematis")
plotTree(tree, fsize = 0.7, lwd =1, ftype = "i")

#Export 6x6

#tree <- read.nexus("Excel Sheets/MyNexusTreefile.nex")
#plot(tree)

#### Data ####
metadata <- read.csv("Excel Sheets/metafull.csv", header=T)
names(metadata)[5] <- "animal"
labels <- read.csv("Excel Sheets/species.csv", header=F)
hit <- sample(which(metadata$animal == labels[1,1]), 1)
dat.prune <- metadata[hit, ]
for(i in 2:57){
  hit <- which(metadata$animal == labels[i,1])
  if(length(hit)>1) {hit <- sample(hit, 1)
  dat.prune <- rbind(dat.prune, metadata[hit, ])
  } else {
    dat.prune <- rbind(dat.prune, metadata[hit, ])
  }
}

#fit <- pglmm(Fisher_Z ~ 1, random.effects = list(~1|Authors, ~1|Pattern, ~1|Age, 
#                                         ~1|Sex, ~1|Location, ~1|Season), 
#             data=dat.prune, REML = T, tree = tree, sp = species, site = site)

################# MCMC ##############
#for MCMCglmm to run we need to change our edge lengths from 0 to 0.00001
print(tree$edge.length)
for (i in 1:length(tree$edge.length)){
  if (tree$edge.length[i] == 0){
    tree$edge.length[i] <- 0.00001
  } else {
    print(i)
  }
}

#Check to see if the labels match
tree$tip.label<-gsub("_"," ",tree$tip.label)
tips <- TipLabels(tree)
tips %in% dat.prune$animal 

#Run the model!
full.mcmc <- MCMCglmm(Fisher_Z ~ 1, random = ~animal,
                      data=dat.prune, pedigree = tree, nitt = 20000, thin = 50, 
                      burnin = 10000)
#random = Authors + Species + Pattern + Age + 
#  Sex + Location + Season
summary(full.mcmc)
plot(full.mcmc$Sol)
