#MCMCglmm using Fisher Z not rho

setwd("C:\\Users\\sarah\\Documents\\meta-analysis")

library(ape)
library(MCMCglmm)
library(phytools)
library(TreeTools)
library(ggplot2)
library(emmeans)
library(bayestestR)
library(meta)
library(DescTools)

#### Phylogeny ####
#First, we need to read in the tree and root it to convert it from an 
#ultrametric tree to a rooted tree for the analysis
#Our outgroup is the sea anemone, Phymactis clematis 
tree <- read.tree("Excel Sheets/list.nwk")
tree <- root(tree, "Phymactis_clematis")
#Force ultrametric because tree due to rounding issues
tree <- force.ultrametric(tree)

#check if the tree is ultrametric and rooted for analysis
is.ultrametric(tree)
# TRUE
is.rooted(tree)
# TRUE
#plot the tree
plotTree(tree, fsize = 0.7, lwd =1, ftype = "i")

#Export 6x6
#To clear the tree so we can plot the rest of our figures:
dev.off()

#for MCMCglmm to run we need to change our edge lengths from 0 to 0.000001
for (i in 1:length(tree$edge.length)){
  if (tree$edge.length[i] == 0){
    tree$edge.length[i] <- 0.000001
  } else {
    print(i)
  }
}
print(tree$edge.length)

#### Data ####
#read in the dataset with the calculated effect sizes
metadata <- read.csv("Excel Sheets/meta_complete_data2.csv", header=T)

#Change Eu_Phae to represent all pigments
metadata$Eu_Pheomelanin <- ifelse(metadata$Eu_Pheomelanin=="N/A", metadata$Classification, metadata$Eu_Pheomelanin)
unique(metadata$Eu_Pheomelanin)

#we must change species to animal for the analysis to run
names(metadata)[5] <- "animal"
#Check to see if the labels on our tree match the animal column
tree$tip.label<-gsub("_"," ",tree$tip.label)
tips <- TipLabels(tree)
tips %in% metadata$animal
metadata$animal %in% tips

#store for later figures
metadatas <- metadata
#remove pterinidine point (only one)
metadata <- metadata[c(1:134, 136:170),]

#Change variables to be factors not characters
metadata[,c(2:7,12:23)] <- lapply(metadata[,c(2:7,12:23)], factor)

################ RANDOM EFFECTS MODEL ####
load("R Files/allRnd_Z.RDATA")

#priors:
prior.ex <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model:
allRnd.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                        data=metadatas, pedigree = tree, 
                        nitt = 5000000, thin = 1000, burnin = 2500000, 
                        prior = prior.ex)


#Save the model for later
#save(allRnd.Z, file = "allRnd_Z.RDATA")

#Summary of Results
summary(allRnd.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -283.4245  

#Location effects: Fisher_Z ~ 1  
#             post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)    0.24843  0.04359  0.47720     2500 0.0256 *


plot(allRnd.Z$Sol)

#proportion of samples above 0
p_significance(allRnd.Z$Sol,threshold = 0)
#Practical Significance (threshold: 0.00)
#Parameter   |   ps
#(Intercept) | 0.99


# Proportion of variance explained by random factors
rand <- allRnd.Z$VCV/apply(allRnd.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

###### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.007306089     0.025687274     0.965300234     0.001706402

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 2.568727  

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.7306089

## total heterogeneity percent
(I2s*100)+(I2u*100)
#3.299336

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2105515

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(allRnd.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.24842743 0.04358969 0.47719709  

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$Weight
MR<-metadata$Fisher_Z-pred_matrix[1:169]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Call:
#glm(formula = zMR ~ Precision, family = "gaussian", data = metadata)
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.70269    0.37245   1.887   0.0609 .
#Precision   -0.11253    0.04491  -2.506   0.0132 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Trim_Fill <- meta::trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 169 (with 0 added studies)

#                                         95%-CI     z  p-value
#Random effects model   0.0147 [-0.0599; 0.0894] 0.39  0.6992

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.248)
FisherZInv(0.044)
FisherZInv(0.477)


FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = "black", pch = 16,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Data Points", "Trim and Fill Points")
points <- c(16,8)
ys <- c(20.1,19.3)
for(i in 1:2){
  points(x=-2.9, y=ys[i], pch=points[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Run the model without tree ####
priors <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                        G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
                        G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
               R = list(V=1, nu=0.02))

tree.Z <- MCMCglmm(Fisher_Z ~ 1, 
                      random = ~ animal + Authors + us(SE_Z):units, 
                      data=metadata, 
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = priors)
tree.Z$DIC
#-282.069

#Save the model for later
#save(tree.Z, file = "tree_Z.RDATA")

#### Run the model without tree and species ####
priors <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))

species.Z <- MCMCglmm(Fisher_Z ~ 1, 
                    random = ~ Authors + us(SE_Z):units, 
                    data=metadata, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = priors)
species.Z$DIC
#-282.1048

#Save the model for later
#save(species.Z, file = "species_Z.RDATA")

#### Run the model without Authors/Study ####
au.Z <- MCMCglmm(Fisher_Z ~ 1, 
                 random = ~ animal + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = priors)
au.Z$DIC
#-265.5481

#Save the model for later
#save(au.Z, file = "au_Z.RDATA")

#### Run the model without weight ####
weight.Z <- MCMCglmm(Fisher_Z ~ 1, 
                 random = ~ animal + Authors, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = priors)
weight.Z$DIC
#175.2044

################ RANDOM EFFECTS SUBSET BY ####
#### Run the model subset by social rank uncontrolled removed ####
#Subset data by removing papers that did not control for social rank
rankcon <- metadatas[metadatas$Social_Rank_Controlled != "uncontrolled",]

#load the model
load("R Files/rsocial_Z.RDATA")

rsocial.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                     data=rankcon, pedigree = tree, 
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex)
rsocial.Z$DIC
#-221.3196

#Save the model
save(rsocial.Z, file = "rsocial_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rsocial.Z$VCV/apply(rsocial.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.007249351     0.033621569     0.954778062     0.004351018

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 3.362157   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.7249351 

## total heterogeneity percent
(I2s*100)+(I2u*100)
#4.087092 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.1603061

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rsocial.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.1949752 0.0472759 0.3356330   

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-rankcon$Weight
MR<-rankcon$Fisher_Z-pred_matrix[1:133]
zMR<-MR*Precision
rankcon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=rankcon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept) -0.119553   0.305016  -0.392    0.696
#Precision    0.001078   0.034706   0.031    0.975
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Trim_Fill <- meta::trimfill(MR, rankcon$SE_Z)
Trim_Fill

#Number of studies: k = 145 (with 12 added studies)

#                                         95%-CI     z  p-value
#Random effects model  0.0267 [-0.0252; 0.0785] 1.01  0.3133

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.1949752)
FisherZInv(0.0472759)
FisherZInv(0.3356330)

cols <- c(rep("black", 133), rep("blue", 12))
p <- c(rep(16, 133), rep(15,12))
FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = p,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Data Points", "Trim and Fill Points")
points <- c(16,15)
cols <- c("black", "blue")
ys <- c(20.1,19.3)
for(i in 1:2){
  points(x=-2.9, y=ys[i], pch=points[i], col =  cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Run the model subset by social rank uncontrolled removed and color class ####
#load the model
load("R Files/rsocolor_Z.RDATA")

prior.ex2 <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
                           G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))

rsocolor.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, 
                       random = ~animal + Authors + us(SE_Z):units,
                      data=rankcon, pedigree = tree, 
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex2)
rsocolor.Z$DIC
#-214.4816

#Save the model
save(rsocolor.Z, file = "rsocolor_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rsocolor.Z$VCV/apply(rsocolor.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.013228267     0.046839832     0.934703746     0.005228155

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
#4.683983     

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.322827

## total heterogeneity percent
(I2s*100)+(I2u*100)
#6.00681 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2025885

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rsocolor.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.20562785 -0.04057176  0.45300421  

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-rankcon$Weight
MR<-rankcon$Fisher_Z-pred_matrix[1:133]
zMR<-MR*Precision
rankcon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=rankcon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept) -0.06733    0.30131  -0.223    0.824
#Precision   -0.01685    0.03428  -0.491    0.624
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Trim and Fill is only used for funnel plot creation since Egger's regression is non-sign

Trim_Fill <- meta::trimfill(MR, rankcon$SE_Z)
Trim_Fill

#Number of studies combined: k = 144 (with 11 added studies)

#                                         95%-CI     z  p-value
#Random effects model  0.0133 [-0.0381; 0.0647] 0.51  0.6123

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.20562785)
FisherZInv(-0.04057176)
FisherZInv(0.45300421)


mycol <- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
cols <- c(rep("black", 124), rep(mycol, 9))
p <- c(rep(16, 124), rep(15,9))
FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = p,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#Export 8x8


#### Run the model subset by condition uncontrolled removed ####
condcon <- metadatas[metadatas$Condition != "None measured",]

#load the model
load("R Files/rcond_Z.RDATA")

rcond.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                      data=condcon, pedigree = tree, 
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex)
rcond.Z$DIC
#-127.0492

#Save the model
#save(rcond.Z, file = "rcond_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rcond.Z$VCV/apply(rcond.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.03908488      0.14313502      0.80315777      0.01462234

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 14.3135    

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
#3.908488 

## total heterogeneity percent
(I2s*100)+(I2u*100)
#18.22199 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.1985594 

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rcond.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.3096058 0.0149827 0.5863683  

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-condcon$Weight
MR<-condcon$Fisher_Z-pred_matrix[1:112]
zMR<-MR*Precision
condcon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=condcon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   1.1802     0.5162   2.286 0.024153 *  
#Precision    -0.2755     0.0712  -3.869 0.000186 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Trim_Fill <- meta::trimfill(MR, condcon$SE_Z)
Trim_Fill

#Number of studies combined: k = 122 (with 10 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.1535 [-0.2477; -0.0593] -3.19  0.0014

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.3096058)
FisherZInv(0.0149827)
FisherZInv(0.5863683)
FisherZInv(0.3096058-0.1535)
FisherZInv(0.0149827-0.1535)
FisherZInv(0.5863683-0.1535)

cols <- c(rep("black", 101), rep("blue", 10))
p <- c(rep(16, 101), rep(15, 10))
FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = p,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Data Points", "Trim and Fill Points")
points <- c(16,15)
cols <- c("black", "blue")
ys <- c(20.1,19.3)
for(i in 1:2){
  points(x=-2.9, y=ys[i], pch=points[i], col =  cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Run the model subset by condition uncontrolled removed and color class ####
#load the model
load("R Files/rcondcolor_Z.RDATA")

rcondcolor.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, random = ~animal + Authors + us(SE_Z):units,
                    data=condcon, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
rcondcolor.Z$DIC
#-123.7725

#Save the model
#save(rcondcolor.Z, file = "rcondcolor_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rcondcolor.Z$VCV/apply(rcondcolor.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.04998185      0.14717948      0.78790836      0.01493032

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 14.71795   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 4.998185 

## total heterogeneity percent
(I2s*100)+(I2u*100)
#19.71613  

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.23566165

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rcondcolor.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.32098198 -0.04149762  0.71388893 

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-condcon$Weight
MR<-condcon$Fisher_Z-pred_matrix[1:112]
zMR<-MR*Precision
condcon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=condcon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.73132    0.51321   1.425  0.15699   
#Precision   -0.20735    0.07079  -2.929  0.00413 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Trim_Fill <- meta::trimfill(MR, condcon$SE_Z)
Trim_Fill

#Number of studies combined: k = 120 (with 8 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.1460 [-0.2364; -0.0557] -3.17  0.0015

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.32098198)
FisherZInv(-0.04149762)
FisherZInv(0.71388893)
FisherZInv(0.32098198-0.1460)
FisherZInv(-0.04149762-0.1460)
FisherZInv(0.71388893-0.1460)

mycol <- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
cols <- c(rep("black", 101), rep(mycol, 7))
p <- c(rep(16, 101), rep(15, 7))
FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = p,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

#Export 8x8


#### Run the model subset by social rank and condition uncontrolled removed ####
randccon <- rankcon[rankcon$Condition != "None measured",]

#load the model
load("R Files/rsandc_Z.RDATA")

rsandc.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                      data=randccon, pedigree = tree, 
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex)
rsandc.Z$DIC
#-114.3896

#Save the model
#save(rsandc.Z, file = "rsandc_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rsandc.Z$VCV/apply(rsandc.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.009821188     0.029387298     0.952688840     0.008102675 

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
#2.93873    

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.9821188 

## total heterogeneity percent
(I2s*100)+(I2u*100)
#3.920849   

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2075871 

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rsandc.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.23155547 0.07943347 0.41159614 

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-randccon$Weight
MR<-randccon$Fisher_Z-pred_matrix[1:87]
zMR<-MR*Precision
randccon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=randccon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   1.2909     0.4632   2.787 0.006569 ** 
#Precision    -0.2163     0.0619  -3.494 0.000757 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Trim_Fill <- meta::trimfill(MR, randccon$SE_Z)
Trim_Fill

#Number of studies combined: k = 94 (with 7 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.0673 [-0.1426; 0.0080] -1.75  0.0798

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.23155547)
FisherZInv(0.07943347)
FisherZInv(0.41159614)
FisherZInv(0.23155547-0.0673)
FisherZInv(0.07943347-0.0673)
FisherZInv(0.41159614-0.0673)

cols <- c(rep("black", 78), rep("blue", 9))
p <- c(rep(16, 78), rep(15, 9))
FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = p,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Data Points", "Trim and Fill Points")
points <- c(16,15)
cols <- c("black", "blue")
ys <- c(20.1,19.3)
for(i in 1:2){
  points(x=-2.9, y=ys[i], pch=points[i], col =  cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Run the model subset by social rank and condition uncontrolled removed and color####
#load the model
load("R Files/rsandcolor_Z.RDATA")

rsandcolor.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, random = ~animal + Authors + us(SE_Z):units,
                     data=randccon, pedigree = tree, 
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)
rsandcolor.Z$DIC
#-111.6156

#Save the model
save(rsandcolor.Z, file = "rsandcolor_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rsandcolor.Z$VCV/apply(rsandcolor.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.011952785     0.034955841     0.944591379     0.008499994

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 3.495584   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
#1.195279  

## total heterogeneity percent
(I2s*100)+(I2u*100)
#4.690863  

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2157207 

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rsandcolor.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.232823226 -0.007874676  0.483511010

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-randccon$Weight
MR<-randccon$Fisher_Z-pred_matrix[1:87]
zMR<-MR*Precision
randccon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=randccon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   1.2478     0.4550   2.742 0.007437 ** 
#Precision    -0.2093     0.0608  -3.442 0.000897 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Trim_Fill <- meta::trimfill(MR, randccon$SE_Z)
Trim_Fill

#Number of studies combined: k = 95 (with 8 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.0737 [-0.1478; 0.0005] -1.95  0.0516

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.232823226)
FisherZInv(-0.007874676)
FisherZInv(0.483511010)
FisherZInv(0.232823226-0.0737)
FisherZInv(-0.007874676-0.0737)
FisherZInv(0.483511010-0.0737)

cols <- c(rep("black", 78), rep("blue", 8))
p <- c(rep(16, 78), rep(15, 8))
FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = p,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Data Points", "Trim and Fill Points")
points <- c(16,15)
cols <- c("black", "blue")
ys <- c(20.1,19.3)
for(i in 1:2){
  points(x=-2.9, y=ys[i], pch=points[i], col =  cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Run the model subset by age uncontrolled removed ####
#Subset data by removing papers that did not control for age
agecon <- metadatas[metadatas$Age_Controlled != "uncontrolled",]

#load the model
load("R Files/rage_Z.RDATA")

rage.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                      data=agecon, pedigree = tree, 
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex)
rage.Z$DIC
#14.1832

#Save the model
#save(rage.Z, file = "rage_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rage.Z$VCV/apply(rage.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.1584632       0.3573036       0.3483839       0.1358493 

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 35.73036    

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 15.8463

## total heterogeneity percent
(I2s*100)+(I2u*100)
#51.57668 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2431849

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rage.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2598076 -0.4301027  1.0783987  

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-agecon$Weight
MR<-agecon$Fisher_Z-pred_matrix[1:25]
zMR<-MR*Precision
agecon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=agecon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   0.5188     1.6742   0.310    0.759
#Precision    -0.1730     0.2526  -0.685    0.500
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Trim_Fill <- meta::trimfill(MR, agecon$SE_Z)
Trim_Fill

#Number of studies combined: k = 25 (with 0 added studies)

#                                         95%-CI     z  p-value
#Random effects model  -0.0871 [-0.3202; 0.1459] -0.73  0.4638



FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = "black", pch = 16,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#Export 8x8


#### Run the model subset by age uncontrolled removed and color class####

#load the model
load("R Files/ragec_Z.RDATA")

ragec.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, random = ~animal + Authors + us(SE_Z):units,
                   data=agecon, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex)
ragec.Z$DIC
#12.20796

#Save the model
#save(ragec.Z, file = "ragec_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- ragec.Z$VCV/apply(ragec.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.2907266       0.3401662       0.2801957       0.0889115 

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 34.01662    

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 29.07266

## total heterogeneity percent
(I2s*100)+(I2u*100)
#51.57668 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.4038967

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(ragec.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2891451 -1.2252288  1.8440485  

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-agecon$Weight
MR<-agecon$Fisher_Z-pred_matrix[1:25]
zMR<-MR*Precision
agecon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=agecon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   1.5899     1.6885   0.942    0.356
#Precision    -0.3957     0.2548  -1.553    0.134
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Trim_Fill <- meta::trimfill(MR, agecon$SE_Z)
Trim_Fill

#Number of studies combined: k = 25 (with 0 added studies)

#                                         95%-CI     z  p-value
#Random effects model  -0.1161 [-0.3516; 0.1195] -0.97  0.3342



FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = "black", pch = 16,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#Export 8x8


#### Run the model subset by age, social rank, and condition uncontrolled removed and color class####
raccon <- rankcon[rankcon$Age_Controlled != "uncontrolled",]
raccon <- raccon[raccon$Condition != "None measured",]

#load the model
load("R Files/rankagecond_Z.RDATA")

rankagecond.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                    data=raccon, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex)
rankagecond.Z$DIC
#-4.435873

#Save the model
save(rankagecond.Z, file = "rankagecond_Z.RDATA")

###### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rankagecond.Z$VCV/apply(rankagecond.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.10442912      0.05032449      0.57606732      0.26917907

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samlng error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 5.032449     

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 10.44291

## total heterogeneity percent
(I2s*100)+(I2u*100)
#15.47536 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2463342

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rankagecond.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.1769460 -0.1807661  0.5232827  

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-raccon$Weight
MR<-raccon$Fisher_Z-pred_matrix[1:15]
zMR<-MR*Precision
raccon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=raccon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   1.6627     0.7944   2.093   0.0565 .
#Precision    -0.3286     0.1379  -2.383   0.0331 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Trim_Fill <- meta::trimfill(MR, raccon$SE_Z)
Trim_Fill

#Number of studies combined: k = 18 (with 3 added studies)

#                                         95%-CI     z  p-value
#Random effects model  -0.0926 [-0.2012; 0.0159] -1.67  0.0945



FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = "black", pch = 16,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#Export 8x8


### MIXED EFFECTS MODELS ####
#### Run the model with Class only ####
prior.ex2 <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
                           G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))

class.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree,
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)
class.Z$DIC
#-279.3766

#Save the model for later
save(class.Z, file = "class_Z.RDATA")

load("R Files/class_Z.RDATA")

color.Z <- class.Z

summary(color.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -279.3766

#Location effects: Fisher_Z ~ Eu_Pheomelanin - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#Eu_Pheomelanincarotenoid    0.14825 -0.11504  0.41019     2500 0.2224  
#Eu_Pheomelanineumelanin     0.31717  0.07431  0.58482     2500 0.0160 *
#Eu_Pheomelaninpheomelanin   0.24738 -0.11430  0.63211     2500 0.1544  
#Eu_Pheomelaninstructural    0.26541 -0.15531  0.75116     2347 0.2304  
#Eu_Pheomelaninunknown       0.37169  0.01358  0.71101     2500 0.0448 *
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(color.Z$Sol)
autocorr(color.Z$Sol)
autocorr(color.Z$VCV)

xsim <- simulate(color.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- color.Z$VCV/apply(color.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

###### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.009755790     0.025980433     0.962570925     0.001692852

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samling error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogeneity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
# 2.598043 

## species level heterogeneity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
#0.975579  

## total heterogeneity
(I2s*100) + (I2u*100)
#3.573622 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2606474 

# Proportion of variance explained by random factors
Sol <- color.Z$Sol/apply(color.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

###### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(color.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2573985 -0.7637793  1.2765428 

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(color.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2612630 -0.0288979  0.5645946   

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$Weight
MR<-metadata$Fisher_Z-pred_matrix[1:169]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.49888    0.36251   1.376   0.1706  
#Precision   -0.08963    0.04371  -2.051   0.0419 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#### Since the intercept of Egger's regression is not significant 
#### there is no publication bias. Trim and fill is only used for 
#### funnel plot creation

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 169 (with 0 added studies)

#                                       95%-CI    z p-value
#Random effects model 0.0021 [-0.0690; 0.0732] 0.06  0.9534

FisherZInv(0.2612630)
FisherZInv(-0.0288979)
FisherZInv(0.5645946)

cols <- rep(0,158)
for (i in 1:158){
  if (i > 158){
    cols[i] <- "forestgreen"
  } else if (metadatas$Classification[i] == "carotenoid"){
    cols[i] <- "orange"
  } else if (metadatas$Eu_Pheomelanin[i] == "eumelanin"){
    cols[i] <- "black"
  } else if (metadatas$Eu_Pheomelanin[i] == "pheomelanin"){
    cols[i] <- "orangered3"
  } else if (metadatas$Classification[i] == "unknown"){
    cols[i] <- "darkorchid4"
  } else if (metadatas$Classification[i] == "structural"){
    cols[i] <- "cornflowerblue"
  } else if (metadatas$Classification[i] == "pteridine"){
    cols[i] <- "deeppink"
  } 
}

shapes <- rep(0,167)
for (i in 1:167){
  if (i >= 158){
    shapes[i] <- 8
  } else if (metadatas$Plasticity[i] == "Plastic"){
    shapes[i] <- 20
  } else if (metadatas$Plasticity[i] == "No"){
    shapes[i] <- 18
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown",
           "Structural", "Pteridine",  "Plastic", "Non-Plastic", 
           "Trim and Fill Points")
Cols <- c("orange","black", "orangered3", "darkorchid4", 
          "cornflowerblue", "deeppink", "grey50", "grey50", 
          "forestgreen")
points <- c(15,15,15,15,15,15,20,18,8)
ys <- c(20.1,19.3,18.5,17.7,16.9,16.1,15.3,14.5,13.7)
for(i in 1:8){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


###### Medians and 95% Credible Intervals ####
emmeans(color.Z, ~ Eu_Pheomelanin, data=metadata, level = 0.95)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.150   -0.1150     0.410
#eumelanin       0.309    0.0743     0.585
#pheomelanin     0.239   -0.1143     0.632
#structural      0.261   -0.1553     0.751
#unknown         0.375    0.0136     0.711
#Point estimate displayed: median 
#HPD interval probability: 0.95
Zcar <- 0.150; lcar <- -0.1150; ucar <- 0.410
Zeu <- 0.309;  leu <-   0.0743; ueu <- 0.585
Zph <- 0.239;  lph <-  -0.1143; uph <- 0.632
Zst <- 0.261;  lst <-  -0.1553; ust <- 0.751
Zun <- 0.375;  lun <-   0.0136; uun <- 0.711

#Probability above 0
p_significance(color.Z$Sol, threshold = 0)
#Parameter                 |   ps
#Eu_Pheomelanincarotenoid  | 0.89
#Eu_Pheomelanineumelanin   | 0.99
#Eu_Pheomelaninpheomelanin | 0.92
#Eu_Pheomelaninstructural  | 0.88
#Eu_Pheomelaninunknown     | 0.98


#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.7,1.9),axes=F,ann=F)
axis(1)
#Fisher Z
#Overall Model
segments(-0.03201199,1.8,0.58309319,1.8)
points(0.25903120,1.8,pch=16,col="black",xpd=NA)
#Carotenoid
segments(lcar,1.6,ucar,1.6);
points(Zcar,1.6,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu,1.4,ueu,1.4);
points(Zeu,1.4,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph,1.2,uph,1.2);
points(Zph,1.2,pch=16,col = "black",xpd=NA)
#Structural
segments(lst,1,ust,1);
points(Zst,1,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun,0.8,uun,0.8);
points(Zun,0.8,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.8,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,1.6,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Eumelanin*", cex = 0.9, adj = c(0,0), font=2)
text(-2.1,1.2,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Structural", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Unknown*", cex = 0.9, adj = c(0,0), font = 2)

#Export 6x6
#### Run the model with Class No Unknowns coded as Eumelanic ####
metadat <- metadata
metadat$Eu_Pheomelanin <- ifelse(metadata$Eu_Pheomelanin =="unknown", yes = "eumelanin", no = metadata$Eu_Pheomelanin)
unique(metadat$Eu_Pheomelanin)

cnounk.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadat, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
cnounk.Z$DIC
#-279.5898

#Save the model for later
#save(cnounk.Z, file = "cnounk_Z.RDATA")

# Proportion of variance explained by random factors
rand <- cnounk.Z$VCV/apply(cnounk.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

###### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.009404615     0.025238146     0.963618685     0.001738554

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samling error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogeneity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
#  2.523815 

## species level heterogeneity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.9404615

## total heterogeneity
(I2s*100) + (I2u*100)
#3.464276  

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2585012 

# Proportion of variance explained by random factors
Sol <- cnounk.Z$Sol/apply(cnounk.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

###### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(cnounk.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2631020 -0.7371434  1.2710833

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(cnounk.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.25997369 -0.01171932  0.54473681  

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadat$Weight
MR<-metadat$Fisher_Z-pred_matrix[1:169]
zMR<-MR*Precision
metadat[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadat)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.51585    0.36471   1.414   0.1591  
#Precision   -0.09137    0.04397  -2.078   0.0393 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Trim_Fill <- trimfill(MR, metadat$SE_Z)
Trim_Fill

#Number of studies combined: k = 169 (with 0 added studies)

#                                       95%-CI    z p-value
#Random effects model 0.0035 [-0.0684; 0.0753] 0.09  0.9247

FisherZInv(0.25997369)
FisherZInv(-0.01171932)
FisherZInv(0.54473681)
FisherZInv(0.25997369+0.0035)
FisherZInv(-0.01171932+0.0035)
FisherZInv(0.54473681+0.0035)

cols <- rep(0,158)
for (i in 1:158){
  if (i > 158){
    cols[i] <- "forestgreen"
  } else if (metadatas$Classification[i] == "carotenoid"){
    cols[i] <- "orange"
  } else if (metadatas$Eu_Pheomelanin[i] == "eumelanin"){
    cols[i] <- "black"
  } else if (metadatas$Eu_Pheomelanin[i] == "pheomelanin"){
    cols[i] <- "orangered3"
  } else if (metadatas$Classification[i] == "unknown"){
    cols[i] <- "darkorchid4"
  } else if (metadatas$Classification[i] == "structural"){
    cols[i] <- "cornflowerblue"
  } else if (metadatas$Classification[i] == "pteridine"){
    cols[i] <- "deeppink"
  } 
}

shapes <- rep(0,167)
for (i in 1:167){
  if (i >= 158){
    shapes[i] <- 8
  } else if (metadatas$Plasticity[i] == "Plastic"){
    shapes[i] <- 20
  } else if (metadatas$Plasticity[i] == "No"){
    shapes[i] <- 18
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown",
           "Structural", "Pteridine",  "Plastic", "Non-Plastic", 
           "Trim and Fill Points")
Cols <- c("orange","black", "orangered3", "darkorchid4", 
          "cornflowerblue", "deeppink", "grey50", "grey50", 
          "forestgreen")
points <- c(15,15,15,15,15,15,20,18,8)
ys <- c(20.1,19.3,18.5,17.7,16.9,16.1,15.3,14.5,13.7)
for(i in 1:8){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


###### Medians and 95% Credible Intervals ####
emmeans(cnounk.Z, ~ Eu_Pheomelanin, data=metadat, level = 0.95)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.144   -0.0999     0.429
#eumelanin       0.319    0.0868     0.579
#pheomelanin     0.241   -0.1105     0.605
#structural      0.258   -0.1962     0.701
#Point estimate displayed: median 
#HPD interval probability: 0.95
Zcar <- 0.144; lcar <- -0.0999; ucar <- 0.429
Zeu <- 0.319;  leu <-   0.0868; ueu <- 0.579
Zph <- 0.241;  lph <-  -0.1105; uph <- 0.605
Zst <- 0.258;  lst <-  -0.1962; ust <- 0.701


#Probability above 0
p_significance(cnounk.Z$Sol, threshold = 0)
#Parameter                 |   ps
#Eu_Pheomelanincarotenoid  | 0.89
#Eu_Pheomelanineumelanin   | 0.99
#Eu_Pheomelaninpheomelanin | 0.92
#Eu_Pheomelaninstructural  | 0.89



#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.9,1.9),axes=F,ann=F)
axis(1)
#Fisher Z
#Overall Model
segments(-0.03201199,1.8,0.58309319,1.8)
points(0.25903120,1.8,pch=16,col="black",xpd=NA)
#Carotenoid
segments(lcar,1.6,ucar,1.6);
points(Zcar,1.6,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu,1.4,ueu,1.4);
points(Zeu,1.4,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph,1.2,uph,1.2);
points(Zph,1.2,pch=16,col = "black",xpd=NA)
#Structural
segments(lst,1,ust,1);
points(Zst,1,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.8,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,1.6,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Eumelanin*", cex = 0.9, adj = c(0,0), font=2)
text(-2.1,1.2,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Structural", cex = 0.9, adj = c(0,0))


#Export 6x6
#### Run the model with Aggression Measures only ####
aggr.Z <- MCMCglmm(Fisher_Z ~ Aggression - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
aggr.Z$DIC
#-282.6515

#Save the model
save(aggr.Z, file = "agg_Z.RDATA")

load("R Files/agg_Z.RDATA")

summary(aggr.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC:  -282.6515
#Location effects: Fisher_Z ~ Aggression - 1 
#                   post.mean l-95% CI u-95% CI eff.samp  pMCMC   
#AggressionDirect     0.30374  0.10199  0.53342     2701 0.0176 *
#AggressionIndirect   0.19190 -0.01917  0.40682     2516 0.0640 .
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Model diagnostics
plot(aggr.Z$Sol)
autocorr(aggr.Z$Sol)
autocorr(aggr.Z$VCV)

xsim <- simulate(aggr.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- aggr.Z$VCV/apply(aggr.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

###### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.006068031     0.025983939     0.966256711     0.001691319 

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[4]

#varM = samling error effect
varM = Rand_Var[3]

##Total variance equation
varT = varA + varS + varE + varM 

## study level heterogeneity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
#  2.598394

## species level heterogeneity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
#0.6068031

## total heterogeneity
(I2s*100) + (I2u*100)
#3.205197  

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.1798293 

# Proportion of variance explained by random factors
Sol <- aggr.Z$Sol/apply(aggr.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

###### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(aggr.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2487491 -0.7198333  1.2158071

###### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(aggr.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.24682483 0.04033735 0.46899498 

###### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$Weight
MR<-metadata$Fisher_Z-pred_matrix[1:169]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.45546    0.36991   1.231    0.220
#Precision   -0.06855    0.04460  -1.537    0.126
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#### Since the intercept of Egger's regression is not significant 
#### there is no publication bias. Trim and fill is only used for 
#### funnel plot creation

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 169 (with 0 added studies)

#                                       95%-CI    z p-value
#Random effects model 0.0178 [-0.0552; 0.0907] 0.48  0.6331

FisherZInv(0.24682483)
FisherZInv(0.04033735)
FisherZInv(0.46899498)

shapes <- rep(0,159)
for (i in 1:167){
  if (metadatas$Aggression[i] == "Direct"){
    shapes[i] <- 16
  } else if (metadatas$Aggression[i] == "Indirect"){
    shapes[i] <- 18
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Direct Measure", "Indirect Measure")
points <- c(16,18)
ys <- c(20.1,19.3)
for(i in 1:2){
  points(x=-2.9, y=ys[i], pch=points[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


###### Medians and 95% Credible Intervals ####
emmeans(aggr.Z, ~ Aggression, data=metadata, level = 0.95)
#Aggression emmean lower.HPD upper.HPD
#Direct      0.304    0.1020     0.533
#Indirect    0.190   -0.0192     0.407                 

Zdir <- 0.304; ldir <- 0.1020; udir <- 0.533
Zindir <- 0.190; lindir <- -0.0192; uindir <- 0.407

FisherZInv(0.304)
FisherZInv(0.1020)
FisherZInv(0.533)

FisherZInv(0.190)
FisherZInv(-0.0192)
FisherZInv(0.407)

#Probability above 0
p_significance(aggr.Z$Sol, threshold = 0)
#Parameter          |   ps
#AggressionDirect   | 0.99
#AggressionIndirect | 0.97

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(1.3,1.9),axes=F,ann=F)
axis(1)
#Fisher Z
#Overall Model
segments(0.05915659,1.8,0.48136754,1.8)
points(0.25228065,1.8,pch=16,col="black",xpd=NA)
#Direct
segments(ldir,1.6,udir,1.6);
points(Zdir,1.6,pch=16,col = "black",xpd=NA)
#Indirect
segments(lindir,1.4,uindir,1.4);
points(Zindir,1.4,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.8,"Overall Model*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,1.6,"Direct Measure*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,1.4,"Indirect Measure", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Run the model with Age Controlled only ####
#Run the model
acon.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
acon.Z$DIC
#-281.2581


#### Run the model with Plasticity only ####
#Run the model
plasticity.Z <- MCMCglmm(Fisher_Z ~ Plasticity - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
plasticity.Z$DIC
#-281.7001

#Save the model for later
#save(plasticity.Z, file = "plasticity_Z.RDATA")

#### Run the model with Time Lag Analysis ####
#Run the model with linear term
year.Z <- MCMCglmm(Fisher_Z ~ as.numeric(Publication.Year) - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadatas, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
year.Z$DIC
#-283.2241

pred_matrix <- predict(year.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#0.24644261 0.04738101 0.45841392 

#Save the model for later
#save(year.Z, file = "yearlinear_Z.RDATA")

#Run the model with quadratic term
metadatas$Publication.Year2 <- metadatas$Publication.Year^2

prior.y   <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.2, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 1, alpha.mu = 0.2, alpha.V = 0.5),
                           G3 = list(V = 1, nu = 1, alpha.mu = 0.2, alpha.V = 0.5)), 
                  R = list(V=1, nu=1))

year2.Z <- MCMCglmm(as.numeric(Fisher_Z) ~ as.numeric(Publication.Year2) - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadatas, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.y)
year2.Z$DIC
#552.48691

#Save the model for later
save(year2.Z, file = "yearquad_Z.RDATA")

pred_matrix <- predict(year2.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#0.263393213 -0.007684623  0.574176250

#### Run the model with Sex only ####
#Run the model
sex.Z <- MCMCglmm(Fisher_Z ~ Sex - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
sex.Z$DIC
#-280.0899

#Save the model for later
save(sex.Z, file = "sex_Z.RDATA")

#### Run the model with Vert/Invert only ####
#Run the model
vert.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
vert.Z$DIC
#-280.8181

#Save the model for later
#save(vert.Z, file = "vert_Z.RDATA")

#### Run the model with Location only ####
#Run the model
loc.Z <- MCMCglmm(Fisher_Z ~ Location - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
loc.Z$DIC
#-280.4865

#Save the model for later
#save(loc.Z, file = "loc_Z.RDATA")

#### Run the model with Season only ####
#Run the model
season.Z <- MCMCglmm(Fisher_Z ~ Season - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
season.Z$DIC
#-280.0494

#### Run the model with Life Stage only ####
age.Z <- MCMCglmm(Fisher_Z ~ Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
age.Z$DIC
#-279.958

#### Run the model with Geography only ####
geo.Z <- MCMCglmm(Fisher_Z ~ Geographic - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
geo.Z$DIC
#-275.9462

#save model
#save(geo.Z, file = "geo_Z.RDATA")

#### Run the model with Social Rank only ####
social.Z <- MCMCglmm(Fisher_Z ~ Social_Rank_Controlled - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
social.Z$DIC
#-276.8918

#save model
#save(social.Z, file = "social_Z.RDATA")

#### Run the model with Obs vs Exp only ####
obsexp.Z <- MCMCglmm(Fisher_Z ~ Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
obsexp.Z$DIC
#-281.5224

#save model
#save(obsexp.Z, file = "obsexp_Z.RDATA")

#### Run the model with Condition only ####
condition.Z <- MCMCglmm(Fisher_Z ~ Condition - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
condition.Z$DIC
#-280.5764

#save model
#save(condition.Z, file = "condition_Z.RDATA")

######## COMBOS OF TWOS ####
###### Run the model with Class and Vert/Invert ####
#Run the model
cv.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree, 
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)
cv.Z$DIC
#-278.2739

#Save the model for later
#save(cv.Z, file = "cv_Z.RDATA")

###### Run the model with Class and Location ####
#Run the model
cl.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Location - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
cl.Z$DIC
#-278.6673

#Save the model for later
#save(cl.Z, file = "cl_Z.RDATA")

###### Run the model with Class and Season ####
#Run the model
csea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Season - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
csea.Z$DIC
#-278.9705

###### Run the model with Class and Life Stage ####
#Run the model
ca.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Age - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
ca.Z$DIC
#-277.8247

###### Run the model with Plasticity and Life Stage ####
#Run the model
pa.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
pa.Z$DIC
#-279.7642

###### Run the model with Sex and Life Stage ####
#Run the model
sa.Z <- MCMCglmm(Fisher_Z ~ Sex + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
sa.Z$DIC
#-278.7395

###### Run the model with Vert/Invert and Life Stage ####
#Run the model
verta.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
verta.Z$DIC
#-279.4668

###### Run the model with Location and Life Stage ####
#Run the model
la.Z <- MCMCglmm(Fisher_Z ~ Location + Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
la.Z$DIC
#-279.2101

###### Run the model with Season and Life Stage ####
#Run the model
seaa.Z <- MCMCglmm(Fisher_Z ~ Season + Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
seaa.Z$DIC
#-279.6151

###### Run the model with Class and Plasticity ####
#Run the model
pc.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
pc.Z$DIC
#-280.078

###### Run the model with Class and Sex ####
#Run the model
cs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
cs.Z$DIC
#-277.89

###### Run the model with Plasticity and Sex ####
#Run the model
ps.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
ps.Z$DIC
#-280.4701

###### Run the model with Plasticity and Vert/Invert ####
#Run the model
pv.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
pv.Z$DIC
#-281.9165

###### Run the model with Plasticity and Location ####
#Run the model
ploc.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Location - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
ploc.Z$DIC
#-281.8055

#Save the model for later
save(ploc.Z, file = "ploc_Z.RDATA")

###### Run the model with Plasticity and Season ####
#Run the model
psea.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Season - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
psea.Z$DIC
#-279.5517

###### Run the model with Sex and Vert/Invert ####
#Run the model
svert.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
svert.Z$DIC
#-278.8167

###### Run the model with Sex and Location ####
#Run the model
seloc.Z <- MCMCglmm(Fisher_Z ~ Sex + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
seloc.Z$DIC
#-279.6257

###### Run the model with Sex and Season ####
#Run the model
ss.Z <- MCMCglmm(Fisher_Z ~ Sex + Season - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
ss.Z$DIC
#-278.654

###### Run the model with Vert/Invert and Location ####
#Run the model
vl.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
vl.Z$DIC
#-280.3721

###### Run the model with Vert/Invert and Season ####
#Run the model
vsea.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Season - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
vsea.Z$DIC
#-279.9379

###### Run the model with Location and Season ####
#Run the model
sealoc.Z <- MCMCglmm(Fisher_Z ~ Season + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
sealoc.Z$DIC
#-280.3347


###### Run the model with Color and Geography ####
cg.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Geographic - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
cg.Z$DIC
#-274.2157

###### Run the model with Color and Social Rank ####
csoc.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Social_Rank_Controlled - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
csoc.Z$DIC
#-275.2933

###### Run the model with Color and Agg Measure ####
cagg.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Aggression - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
cagg.Z$DIC
#-280.3148

#Save the model for later
#save(cagg.Z, file = "cagg_Z.RDATA")

###### Run the model with Color and Condition ####
ccond.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Condition - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
ccond.Z$DIC
#-279.0012

###### Run the model with Color and Obs vs Exp ####
cobs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Obs_vs_Exp - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
cobs.Z$DIC
#-279.0433

###### Run the model with Agg Measures and Plasticity ####
ap.Z <- MCMCglmm(Fisher_Z ~ Aggression + Plasticity - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

ap.Z$DIC
#-281.9284

#Save the model
save(ap.Z, file = "ap_Z.RDATA")

###### Run the model with Agg Measures and Sex ####
as.Z <- MCMCglmm(Fisher_Z ~ Aggression + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

as.Z$DIC
#-280.2149

###### Run the model with Agg Measures and Vert/Invert ####
av.Z <- MCMCglmm(Fisher_Z ~ Aggression + Vert_Invert - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

av.Z$DIC
#-281.2818

###### Run the model with Agg Measures and Location ####
al.Z <- MCMCglmm(Fisher_Z ~ Aggression + Location - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

al.Z$DIC
#-281.4216

###### Run the model with Agg Measures and Season ####
asea.Z <- MCMCglmm(Fisher_Z ~ Aggression + Season - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

asea.Z$DIC
#-281.2245

###### Run the model with Agg Measures and Life Stage ####
aage.Z <- MCMCglmm(Fisher_Z ~ Aggression + Age - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

aage.Z$DIC
#-281.1896

###### Run the model with Agg Measures and Geography ####
aggeo.Z <- MCMCglmm(Fisher_Z ~ Aggression + Geographic - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

aggeo.Z$DIC
#-276.2135

#Save the model for later
save(aggeo.Z, file = "aggeo_Z.RDATA")

###### Run the model with Agg Measures and Social Rank ####
agsoc.Z <- MCMCglmm(Fisher_Z ~ Aggression + Social_Rank_Controlled - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

agsoc.Z$DIC
#-277.8547

###### Run the model with Agg Measures and Obs vs Exp ####
agob.Z <- MCMCglmm(Fisher_Z ~ Aggression + Obs_vs_Exp - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

agob.Z$DIC
#-281.965

#Save the model for later
#save(agob.Z, file = "agob_Z.RDATA")

###### Run the model with Agg Measures and Condition ####
agcond.Z <- MCMCglmm(Fisher_Z ~ Aggression + Condition - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

agcond.Z$DIC
#-281.0882

###### Run the model with Plasticity and Geography ####
pg.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Geographic - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree,
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)

pg.Z$DIC
#-275.7684

###### Run the model with Plasticity and Social Rank ####
psoc.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Social_Rank_Controlled - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

psoc.Z$DIC
#-278.233

###### Run the model with Plasticity and Obs vs Exp ####
pob.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Obs_vs_Exp - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

pob.Z$DIC
#-281.7183

###### Run the model with Plasticity and Condition ####
pcond.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Condition - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree,
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)

pcond.Z$DIC
#-281.3461

###### Run the model with Sex and Geography ####
sg.Z <- MCMCglmm(Fisher_Z ~ Sex + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

sg.Z$DIC
#-276.2906

###### Run the model with Sex and Social Rank ####
ssoc.Z <- MCMCglmm(Fisher_Z ~ Sex + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

ssoc.Z$DIC
#-277.4439

###### Run the model with Sex and Obs vs Exp ####
sob.Z <- MCMCglmm(Fisher_Z ~ Sex + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

sob.Z$DIC
#-279.73

###### Run the model with Sex and Condition ####
scond.Z <- MCMCglmm(Fisher_Z ~ Sex + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

scond.Z$DIC
#-280.1891

###### Run the model with Vert_Invert and Geography ####
vg.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

vg.Z$DIC
#-275.2178

###### Run the model with Vert_Invert and Social Rank ####
vsoc.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

vsoc.Z$DIC
#-276.8043

###### Run the model with Vert_Invert and Obs vs Exp ####
vob.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

vob.Z$DIC
#-280.5685

###### Run the model with Vert_Invert and Condition ####
vcond.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

vcond.Z$DIC
#-280.4669

###### Run the model with Location and Geography ####
lg.Z <- MCMCglmm(Fisher_Z ~ Location + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

lg.Z$DIC
#-274.3192

###### Run the model with Location and Social Rank ####
lsoc.Z <- MCMCglmm(Fisher_Z ~ Location + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

lsoc.Z$DIC
#-276.936

###### Run the model with Location and Obs vs Exp ####
lob.Z <- MCMCglmm(Fisher_Z ~ Location + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

lob.Z$DIC
#-280.2143

###### Run the model with Location and Condition ####
lcond.Z <- MCMCglmm(Fisher_Z ~ Location + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

lcond.Z$DIC
#-279.7334

###### Run the model with Season and Geography ####
seag.Z <- MCMCglmm(Fisher_Z ~ Season + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

seag.Z$DIC
#-274.6825

###### Run the model with Season and Social Rank ####
seasoc.Z <- MCMCglmm(Fisher_Z ~ Season + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

seasoc.Z$DIC
#-275.4664

###### Run the model with Season and Obs vs Exp ####
seaob.Z <- MCMCglmm(Fisher_Z ~ Season + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

seaob.Z$DIC
#-280.537

###### Run the model with Season and Condition ####
seacond.Z <- MCMCglmm(Fisher_Z ~ Season + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

seacond.Z$DIC
#-274.6474

###### Run the model with Life Stage and Geography ####
ag.Z <- MCMCglmm(Fisher_Z ~ Age + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

ag.Z$DIC
#-275.4199

###### Run the model with Life Stage and Social Rank ####
asoc.Z <- MCMCglmm(Fisher_Z ~ Age + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

asoc.Z$DIC
#-276.195

###### Run the model with Life Stage and Obs vs Exp ####
aob.Z <- MCMCglmm(Fisher_Z ~ Age + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

aob.Z$DIC
#-280.3408

###### Run the model with Life Stage and Condition ####
acond.Z <- MCMCglmm(Fisher_Z ~ Age + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

acond.Z$DIC
#-280.6246

###### Run the model with Geography and Social Rank ####
gsoc.Z <- MCMCglmm(Fisher_Z ~ Social_Rank_Controlled + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

gsoc.Z$DIC
#-273.7864

###### Run the model with Geography and Obs vs Exp ####
gob.Z <- MCMCglmm(Fisher_Z ~ Geographic + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

gob.Z$DIC
#-275.708

###### Run the model with Geography and Condition ####
gcond.Z <- MCMCglmm(Fisher_Z ~ Geographic + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

gcond.Z$DIC
#-275.4908

###### Run the model with Social Rank and Obs vs Exp ####
socob.Z <- MCMCglmm(Fisher_Z ~ Social_Rank_Controlled + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

socob.Z$DIC
#-277.3733

###### Run the model with Social Rank and Condition ####
soccond.Z <- MCMCglmm(Fisher_Z ~ Social_Rank_Controlled + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

soccond.Z$DIC
#-276.6564

###### Run the model with Obs_vs_Exp and Condition ####
obscond.Z <- MCMCglmm(Fisher_Z ~ Obs_vs_Exp + Condition - 1, 
                      random = ~animal + Authors + us(SE_Z):units, 
                      data=metadata, pedigree = tree,
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex2)

obscond.Z$DIC
#-280.5617

###### Run the model with Age controlled and Color ####
agc.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Eu_Pheomelanin - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
agc.Z$DIC
#-278.3787

###### Run the model with Age controlled and Agg Measure ####
agagg.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Aggression - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000,
                  prior = prior.ex2)
agagg.Z$DIC
#-281.6363

###### Run the model with Age controlled and Plasticity ####
agp.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Plasticity - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000,
                    prior = prior.ex2)
agp.Z$DIC
#-281.7649

###### Run the model with Age controlled and Sex ####
ags.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Sex - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000,
                    prior = prior.ex2)
ags.Z$DIC
#-278.3953

###### Run the model with Age controlled and Vert_Invert ####
agv.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Vert_Invert - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000,
                  prior = prior.ex2)
agv.Z$DIC
#-280.3873

###### Run the model with Age controlled and Location ####
agl.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Location - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000,
                  prior = prior.ex2)
agl.Z$DIC
#-279.0037

###### Run the model with Age controlled and Season ####
agsea.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Season - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000,
                  prior = prior.ex2)
agsea.Z$DIC
#-279.9674

###### Run the model with Age controlled and Life Stage ####
agage.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000,
                    prior = prior.ex2)
agage.Z$DIC
#-280.3441

###### Run the model with Age controlled and Geography ####
aggeo.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Geographic - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000,
                    prior = prior.ex2)
aggeo.Z$DIC
#-274.6264

###### Run the model with Age controlled and Social Rank ####
agsoc.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Social_Rank_Controlled - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000,
                    prior = prior.ex2)
agsoc.Z$DIC
#-276.8725

###### Run the model with Age controlled and Obs vs. Exp ####
agobs.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Obs_vs_Exp - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000,
                    prior = prior.ex2)
agobs.Z$DIC
#-280.8294

###### Run the model with Age controlled and Condition ####
agcon.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Obs_vs_Exp - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000,
                    prior = prior.ex2)
agcon.Z$DIC
#-280.034

#### INTERACTION TERMS ####
###### Mixed Effects Model Plasticity, Sex, and Plasticity x Sex ####
#Run the model
pxs.Z <- MCMCglmm(Fisher_Z ~ Sex + Plasticity + Sex*Plasticity - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000,
                  prior = prior.ex2)
pxs.Z$DIC
#-277.3966

###### Mixed Effects Model Class, Plasticity, and Class x Plasticity ####
#Run the model
ebyp.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity 
                   + Eu_Pheomelanin*Plasticity - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000,
                   prior = prior.ex2)
ebyp.Z$DIC
#-275.9657

###### Mixed Effects Model Class, Sex, and Class x Sex ####
#Run the model
ebys.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Sex*Eu_Pheomelanin - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000,
                   prior = prior.ex2)
ebys.Z$DIC
#-274.0648

#### GRAPHS ####
###### Fisher Z with Pub Bias correction subset combined ####
plot(NA,xlim=c(-2,2),ylim=c(-0.1,0.7),axes=F,ann=F)
axis(1)

#### Random Effects Model
#Mean
segments(0.04646076-0.0255,0.6,0.48347851-0.0255,0.6);
points(0.26002580-0.0255,0.6,pch=16,col = "black",xpd=NA)

#### Mixed: Plasticity
#Plastic
segments(lpl-0.0545,0.4,upl-0.0545,0.4);
points(Zpl-0.0545,0.4,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl-0.0545,0.2,unpl-0.0545,0.2);
points(Znpl-0.0545,0.2,pch=16,col = "black",xpd=NA)
#Mean Plasticity Model
segments(-0.051,0,0.467,0);
points(0.207,0,pch=16,col = "black",xpd=NA)

#### Mixed: Plasticity and Vert/Invert
#Plastic
segments(lpla-0.0488,-0.2,upla-0.0488,-0.2);
points(Zpla-0.0488,-0.2,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpla-0.0488,-0.4,unpla-0.0488,-0.4);
points(Znpla-0.0488,-0.4,pch=16,col = "black",xpd=NA)
#Vert
segments(lve-0.0488,-0.6,uve-0.0488,-0.6);
points(Zve-0.0488,-0.6,pch=16,col = "black",xpd=NA)
#Invert
segments(lin-0.0488,-0.8,uin-0.0488,-0.8);
points(Zin-0.0488,-0.8,pch=16,col = "black",xpd=NA)
#Mean Plasticity and Vert/Invert Model
segments(-0.190,-1,0.585,-1);
points(0.208,-1,pch=16,col = "black",xpd=NA)

#Add line at 0 and separate models
abline(h = 0.5, lty = 1)
abline(h = -0.1, lty = 1)
abline(v = 0, lty = 1)

#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Random Effects Model*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,0.4,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.2,"Non-Plastic*", cex = 0.9, adj = c(0,0))
text(-2.1,0,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,-0.2,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,-0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,-0.6,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,-0.8,"Invertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,-1,"Overall Model", cex = 0.9, adj = c(0,0))

#Export 8x8 (6x6 for plasticity and random only)
###### Fisher Z with Publication Bias correction all combined ####
plot(NA,xlim=c(-2,2),ylim=c(-5.1,1.7),axes=F,ann=F)
axis(1)
#### Mixed: Color Class
#Carotenoid
segments(lcar ,1.6,ucar,1.6);
points(Zcar,1.6,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu,1.4,ueu,1.4);
points(Zeu,1.4,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph,1.2,uph,1.2);
points(Zph,1.2,pch=16,col = "black",xpd=NA)
#Structural 
segments(lst,1,ust,1);
points(Zst,1,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun,0.8,uun,0.8);
points(Zun,0.8,pch=16,col = "black",xpd=NA)
#Mean Class Model
segments(-0.084,0.6,0.620,0.6);
points(0.273,0.6,pch=16,col = "black",xpd=NA)

#### Mixed: Plasticity
#Plastic
segments(lpl-0.0545,0.4,upl-0.0545,0.4);
points(Zpl-0.0545,0.4,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl-0.0545,0.2,unpl-0.0545,0.2);
points(Znpl-0.0545,0.2,pch=16,col = "black",xpd=NA)
#Mean Plasticity Model
segments(-0.051,0,0.467,0);
points(0.207,0,pch=16,col = "black",xpd=NA)

#### Mixed: Sex
#Both Sexes
segments(lb-0.0695,-0.2,ub-0.0695,-0.2);
points(Zb-0.0695,-0.2,pch=16,col = "black",xpd=NA)
#Females
segments(lf-0.0695,-0.4,uf-0.0695,-0.4);
points(Zf-0.0695,-0.4,pch=16,col = "black",xpd=NA)
#Males
segments(lm-0.0695,-0.6,um-0.0695,-0.6);
points(Zm-0.0695,-0.6,pch=16,col = "black",xpd=NA)
#Mean Sex Model
segments(-0.050,-0.8,0.472,-0.8);
points(0.193,-0.8,pch=16,col = "black",xpd=NA)

#### Mixed: Plasticity and Sex
#Plastic
segments(lp-0.0927,-1,up-0.0927,-1);
points(Zp-0.0927,-1,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnp-0.0927,-1.2,unp-0.0927,-1.2);
points(Znp-0.0927,-1.2,pch=16,col = "black",xpd=NA)
#Both Sexes
segments(lbs-0.0927,-1.4,ubs-0.0927,-1.4);
points(Zbs-0.0927,-1.4,pch=16,col = "black",xpd=NA)
#Females
segments(lfs-0.0927,-1.6,ufs-0.0927,-1.6);
points(Zfs-0.0927,-1.6,pch=16,col = "black",xpd=NA)
#Males
segments(lms-0.0927,-1.8,ums-0.0927,-1.8);
points(Zms-0.0927,-1.8,pch=16,col = "black",xpd=NA)
#Mean Plasticity and Sex Model
segments(-0.125,-2,0.477,-2);
points(0.172,-2,pch=16,col = "black",xpd=NA)

#### Mixed: Vert/Invert
#Vert
segments(lv-0.0165,-2.2,uv-0.0165,-2.2);
points(Zv-0.0165,-2.2,pch=16,col = "black",xpd=NA)
#Invert
segments(li-0.0165,-2.4,ui-0.0165,-2.4);
points(Zi-0.0165,-2.4,pch=16,col = "black",xpd=NA)
#Mean Vert/Invert Model
segments(-0.065,-2.6,0.595,-2.6);
points(0.2446,-2.6,pch=16,col = "black",xpd=NA)

#### Mixed: Plasticity and Vert/Invert
#Plastic
segments(lpla-0.0488,-2.8,upla-0.0488,-2.8);
points(Zpla-0.0488,-2.8,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpla-0.0488,-3,unpla-0.0488,-3);
points(Znpla-0.0488,-3,pch=16,col = "black",xpd=NA)
#Vert
segments(lve-0.0488,-3.2,uve-0.0488,-3.2);
points(Zve-0.0488,-3.2,pch=16,col = "black",xpd=NA)
#Invert
segments(lin-0.0488,-3.4,uin-0.0488,-3.4);
points(Zin-0.0488,-3.4,pch=16,col = "black",xpd=NA)
#Mean Plasticity and Vert/Invert Model
segments(-0.190,-3.6,0.585,-3.6);
points(0.208,-3.6,pch=16,col = "black",xpd=NA)

#### Plasticity and Location
#Plastic
segments(lplas,-3.8,uplas,-3.8);
points(Zplas,-3.8,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnplas,-4,unplas,-4);
points(Znplas,-4,pch=16,col = "black",xpd=NA)
#Captivity
segments(lcap,-4.2,ucap,-4.2);
points(Zcap,-4.2,pch=16,col = "black",xpd=NA)
#Domestic
segments(ldom,-4.4,udom,-4.4);
points(Zdom,-4.4,pch=16,col = "black",xpd=NA)
#Lab
segments(llab,-4.6,ulab,-4.6);
points(Zlab,-4.6,pch=16,col = "black",xpd=NA)
#Field
segments(lfi,-4.8,ufi,-4.8);
points(Zfi,-4.8,pch=16,col = "black",xpd=NA)
#Mean Plasticity and Vert/Invert Model
segments(-0.036,-5,0.519,-5);
points(0.237,-5,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
abline(h = 0.5, lty = 1)
abline(h = -0.1, lty = 1)
abline(h = -0.9, lty = 1)
abline(h = -2.1, lty = 1)
abline(h = -2.7, lty = 1)
abline(h = -3.7, lty = 1)

#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Eumelanin*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,1.2,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Structural", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Unknown*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,0.6,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.2,"Non-Plastic*", cex = 0.9, adj = c(0,0))
text(-2.1,0,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,-0.2,"Both Sexes*", cex = 0.9, adj = c(0,0))
text(-2.1,-0.4,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,-0.6,"Males*", cex = 0.9, adj = c(0,0))
text(-2.1,-0.8,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,-1,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,-1.2,"Non-Plastic*", cex = 0.9, adj = c(0,0))
text(-2.1,-1.4,"Both Sexes*", cex = 0.9, adj = c(0,0))
text(-2.1,-1.6,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,-1.8,"Males", cex = 0.9, adj = c(0,0))
text(-2.1,-2,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,-2.2,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,-2.4,"Invertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,-2.6,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,-2.8,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,-3,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,-3.2,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,-3.4,"Invertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,-3.6,"Overall Model", cex = 0.9, adj = c(0,0))
text(-2.1,-3.8,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,-4,"Non-Plastic*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-4.2,"Captivity*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-4.4,"Domestic", cex = 0.9, adj = c(0,0))
text(-2.1,-4.6,"Lab*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-4.8,"Field", cex = 0.9, adj = c(0,0))
text(-2.1,-5,"Overall Model", cex = 0.9, adj = c(0,0))

#Export 10x10

###### Fisher Z for each Study #####
order <- read.csv("Excel Sheets/species_order.csv")
order <- cbind(order,rep(0,length(order$Study)),rep(0,length(order$Study)))
names(order)[48] <- "lci"
names(order)[49] <- "uci"

for (i in 1:length(order$Study)){
  order$uci[i] <- order$Fisher_Z[i] + 1.96*order$SE_Z[i]
  order$lci[i] <- order$Fisher_Z[i] - 1.96*order$SE_Z[i]
}

#Fisher Z Plot
colsa <- rep(0,length(order$Classification))
for (i in 1:length(order$Classification)){
  if (order$Classification[i] == "carotenoid"){
    colsa[i] <- "darkorange"
  } else if (order$Eu_Pheomelanin[i] == "eumelanin"){
    colsa[i] <- "black"
  } else if (order$Eu_Pheomelanin[i] == "pheomelanin"){
    colsa[i] <- "tan4"
  } else if (order$Classification[i] == "unknown"){
    colsa[i] <- "slateblue4"
  }  else if (order$Classification[i] == "pteridine"){
    colsa[i] <- "turquoise4"
  } else if (order$Classification[i] == "structural"){
    colsa[i] <- "violet"
  } 
}

shapesa <- rep(0,length(order$Plasticity))
for (i in 1:length(order$Vert_Invert)){
  if (order$Plasticity[i] == "Plastic"){
    shapesa[i] <- 15
  } else if (order$Plasticity[i] == "No"){
    shapesa[i] <- 16
  } 
}

plot(NA,xlim=c(-4,4),ylim=c(0,350),axes=F,ann=F)
axis(1)
polygon(x = c(-4,-4,4,4), y = c(8,14,14,8), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-4,-4,4,4), y = c(16,40,40,16), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-4,-4,4,4), y = c(170,234,234,170), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-4,-4,4,4), y = c(258,342,342,258), col = alpha("gray", 0.8), density = NA)

abline(v=0)
for (i in 1:length(order$Study)){
  segments(order$lci[i],i*2,order$uci[i],i*2);
  points(order$Fisher_Z[i],i*2,pch=shapesa[i], col = colsa[i],xpd=NA)
  text(-4, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
}

#Random effects model no pub bias only need original point
points(x=0.24842743, y = 350,pch = 17)
segments(0.04358969, 350, 0.47719709, 350)

#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Structural", 
           "Pteridine", "Unknown", "Plastic", "Non-Plastic", "Before Correction",
           "After Correction")
Cols <- c("darkorange","black", "tan4", "violet", 
          "turquoise4", "slateblue4", "black", "black", "black","black")
points <- c(15,15,15,15,15,15,15,16, 17, 18)
ys <- c(140,135,130,125,120,115,105,100,95,90)
for(i in 1:10){
  points(x=3.5, y=ys[i], pch=points[i], col=Cols[i])
  text(x=3.5,y=ys[i], labels=words[i], pos=4, cex=.5, font = 2)
}
title(xlab = "Fisher Z")
text(3.5,145,cex = 0.5, "Colors:", font = 2)
text(3.5, 110, cex = 0.5, "Shapes:", font = 2)

#labels for Families
text(3.75,320,cex=0.5,"Actinopterygii", font=2)
text(3.75,235, cex=0.5, "Amphibia", font=2)
text(3.75,213,cex=0.5, "Reptilia", font=2)
text(3.75, 161, cex=0.5, "Aves", font=2)
text(3.75,39,cex=0.5, "Mammalia", font=2)
text(3.75,16,cex=0.5, "Malacostraca", font=2)
text(3.75,13,cex=0.5,"Insecta",font=2)
text(3.75,6,cex=0.5,"Anthoza",font=2)
text(-4.075,348,cex=0.5, pos=4,"Random Effects Model", font=2)


#Export 15x10 and portiate  

###### Fisher Z by Phylogeny ####
order <- read.csv("Excel Sheets/species_order.csv")
order <- cbind(order,rep(0,length(order$Study)),rep(0,length(order$Study)))
names(order)[48] <- "lci"
names(order)[49] <- "uci"

for (i in 1:length(order$Study)){
  order$uci[i] <- order$rho[i] + 1.96*order$SE_r[i]
  order$lci[i] <- order$rho[i] - 1.96*order$SE_r[i]
}

plot(NA,xlim=c(-2,2),ylim=c(1.9,300),axes=F,ann=F)
axis(1)
polygon(x = c(-2,-2,2,2), y = c(8,14,14,8), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-2,-2,2,2), y = c(16,40,40,16), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-2,-2,2,2), y = c(150,196,196,150), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-2,-2,2,2), y = c(220,302,302,220), col = alpha("gray", 0.8), density = NA)

abline(v=0)
for (i in 1:length(order$Study)){
  segments((order$lci[i]),i*2,(order$uci[i]),i*2);
  points((order$rho[i]),i*2,pch=shapesa[i],col = colsa[i],xpd=NA)
  text(-2, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
  text(1.7, i*2, order$Class[i], cex= 0.6, adj = c(0,0), font = 2)
}

#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Structural", 
           "Pteridine", "Unknown", "Plastic", "Non-Plastic")
Cols <- c("darkorange","black", "tan4", "violet", 
          "turquoise4", "slateblue4", "black", "black")
points <- c(15,15,15,15,15,15,15,16)
ys <- c(140,135,130,125,120,115,105,100)
for(i in 1:8){
  points(x=3.75, y=ys[i], pch=points[i], col=Cols[i])
  text(x=3.75,y=ys[i], labels=words[i], pos=4, cex=.5, font = 2)
}
title(xlab = "Fisher Z")
text(3.75,145,cex = 0.5, "Colors:", font = 2)
text(3.75, 110, cex = 0.5, "Shapes:", font = 2)
#Export 12x15

