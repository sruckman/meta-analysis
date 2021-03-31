library(MCMCglmm)
library(ape)

Adams.phylo <- read.nexus("evo_314_sm_tables2.nex",tree.names=T)
Adams.data <- read.csv("evo_314_sm_tables1.csv",header=T)

mammaltree<-read.nexus("evo_314_sm_tables2.nex",tree.names=T)
specieslist<-mammaltree$tip
plot(mammaltree)
sigma<-vcv.phylo(mammaltree,cor=T)
data<-read.csv("evo_314_sm_tables1.csv",header=T)
attach(data)
corr<-as.matrix(data[,(2)])                           
sampsize<-as.matrix(data[,(3)]) 
X<-as.matrix(data[,(4)])
N<-length(corr)
names(corr)<-data[,1]
names(X)<-data[,1]
names(sampsize)<-data[,1]

#Re-order data to match tree order
corr<-corr[match(rownames(sigma),names(corr))]  
X<-X[match(rownames(sigma),names(X))]  
sampsize<-sampsize[match(rownames(sigma),names(sampsize))]  

## effect size calculations for correlation coefficients
FisherZ<- 0.5*log((1+corr)/(1-corr))
effect_var<-1/(sampsize-3)
weight<-1/effect_var
plot(Adams.phylo)

IA <- inverseA(Adams.phylo, nodes = "TIPS")
IAasreml <- sm2asreml(IA$Ainv, IA$node.names)

m10.mcmc <- MCMCglmm(FisherZ ~ 1, random = ~animal, mev = weight,
                     data=Adams.data, pedigree = Adams.phylo)
