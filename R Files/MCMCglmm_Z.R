#MCMCglmm using Fisher Z not rho

setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis")

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
#remove tip as this paper and species were removed because it did not 
#fit the inclusion criteria
tree <- drop.tip(tree, "Luscinia cyanura")

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
metadata <- metadata[c(1:134, 136:159),]


################ Random Effects Model with strict prior ####
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

#DIC: -263.3073 

#Location effects: Fisher_Z ~ 1  
#             post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)    0.25089  0.01804  0.45813     2500 0.0216 *


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

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.006422849     0.022519270     0.969416776     0.001641105

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
# 2.251927   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.6422849

## total heterogeneity percent
(I2s*100)+(I2u*100)
#2.894212

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2100122

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(allRnd.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.25088699 0.01803682 0.45812805  

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$Weight
MR<-metadata$Fisher_Z-pred_matrix[1:149]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Call:
#glm(formula = zMR ~ Precision, family = "gaussian", data = metadata)
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.77740    0.38150   2.038   0.0433 *  
#Precision   -0.11600    0.04523  -2.565   0.0113 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Trim_Fill <- meta::trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 158 (with 0 added studies)

#                                         95%-CI     z  p-value
#Random effects model   0.0233 [-0.0551; 0.1017] 0.58  0.5599

#### Adjusted mean and 95% confidence interval #####
#              fit               lwr              upr 
#0.25088699+0.0233; 0.01803682+0.0233; 0.45812805+0.0233;
#        0.274187        0.04133682         0.481428
0.25088699+0.0233; 0.01803682+0.0233; 0.45812805+0.0233;


#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.25088699)
FisherZInv(0.01803682)
FisherZInv(0.45812805)
FisherZInv(0.25088699+0.0233)
FisherZInv(0.01803682+0.0233)
FisherZInv(0.45812805+0.0233)

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


###### Run the model without tree ####
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
#-262.1931

#Save the model for later
#save(tree.Z, file = "tree_Z.RDATA")

###### Run the model without tree and species ####
priors <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))

species.Z <- MCMCglmm(Fisher_Z ~ 1, 
                    random = ~ Authors + us(SE_Z):units, 
                    data=metadata, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = priors)
species.Z$DIC
#-262.5397

#Save the model for later
#save(species.Z, file = "species_Z.RDATA")

###### Run the model without Authors/Study ####
au.Z <- MCMCglmm(Fisher_Z ~ 1, 
                 random = ~ animal + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = priors)
au.Z$DIC
#-245.6888

#Save the model for later
#save(au.Z, file = "au_Z.RDATA")

###### Run the model without weight ####
weight.Z <- MCMCglmm(Fisher_Z ~ 1, 
                 random = ~ animal + Authors, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = priors)
weight.Z$DIC
#171.2191

################ Random Effects Model subset by ####
###### Run the model subset by social rank uncontrolled removed ####
#Subset data by removing papers that did not control for social rank
rankcon <- metadatas[metadatas$Social_Rank_Controlled != "uncontrolled",]

#load the model
load("R Files/rsocial_Z.RDATA")

rsocial.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                     data=rankcon, pedigree = tree, 
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex)
rsocial.Z$DIC
#-212.0524

#Save the model
#save(rsocial.Z, file = "rsocial_Z.RDATA")

#### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rsocial.Z$VCV/apply(rsocial.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.008603109     0.018581734     0.969337643     0.003477514

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
# 1.858173   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.8603109

## total heterogeneity percent
(I2s*100)+(I2u*100)
#2.718484 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2805756

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rsocial.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.20316297 0.04774556 0.38257041  

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-rankcon$Weight
MR<-rankcon$Fisher_Z-pred_matrix[1:124]
zMR<-MR*Precision
rankcon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=rankcon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  1.25332    0.52901   2.369 0.019769 *  
#Precision   -0.27041    0.07167  -3.773 0.000275 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Trim_Fill <- meta::trimfill(MR, rankcon$SE_Z)
Trim_Fill

#Number of studies combined: k = 129 (with 5 added studies)

#                                         95%-CI     z  p-value
#Random effects model  -0.1044 [-0.1619; -0.0470] -3.56  0.0004

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.20316297)
FisherZInv(0.04774556)
FisherZInv(0.38257041)
FisherZInv(0.20316297-0.1044)
FisherZInv(0.04774556-0.1044)
FisherZInv(0.38257041-0.1044)

cols <- c(rep("black", 124), rep("blue", 5))
p <- c(rep(16, 124), rep(15,5))
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


###### Run the model subset by social rank uncontrolled removed and color class ####
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
#-207.2363

#Save the model
save(rsocolor.Z, file = "rsocolor_Z.RDATA")

#### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rsocolor.Z$VCV/apply(rsocolor.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.014828945     0.023710401     0.957448430     0.004012224 

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
# 2.37104    

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.482894

## total heterogeneity percent
(I2s*100)+(I2u*100)
#3.853935 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.3484935 

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rsocolor.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.21617793 -0.04608399  0.48594556  

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-rankcon$Weight
MR<-rankcon$Fisher_Z-pred_matrix[1:124]
zMR<-MR*Precision
rankcon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=rankcon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  -0.04622    0.29624  -0.156    0.876
#Precision   -0.02097    0.03309  -0.634    0.527
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Trim and Fill is only used for funnel plot creation since Egger's regression is non-sign

Trim_Fill <- meta::trimfill(MR, rankcon$SE_Z)
Trim_Fill

#Number of studies combined: k = 133 (with 9 added studies)

#                                         95%-CI     z  p-value
#Random effects model  0.0056 [-0.0438; 0.0551] 0.22  0.8230

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.21617793)
FisherZInv(-0.04608399)
FisherZInv(0.48594556)


mycol <- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
cols <- c(rep("black", 124), rep(mycol, 9))
p <- c(rep(16, 124), rep(15,9))
FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = p,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#Export 8x8


###### Run the model subset by condition uncontrolled removed ####
condcon <- metadatas[metadatas$Condition != "None measured",]

#load the model
load("R Files/rcond_Z.RDATA")

rcond.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                      data=condcon, pedigree = tree, 
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex)
rcond.Z$DIC
#-113.7745

#Save the model
#save(rcond.Z, file = "rcond_Z.RDATA")

#### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rcond.Z$VCV/apply(rcond.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.02934254      0.12013054      0.83805235      0.01247457

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
# 12.01305   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 2.934254

## total heterogeneity percent
(I2s*100)+(I2u*100)
#14.94731 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.1811853

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rcond.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.30746122 0.03218604 0.59945190  

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-condcon$Weight
MR<-condcon$Fisher_Z-pred_matrix[1:101]
zMR<-MR*Precision
condcon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=condcon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  1.25332    0.52901   2.369 0.019769 *  
#Precision   -0.27041    0.07167  -3.773 0.000275 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Trim_Fill <- meta::trimfill(MR, condcon$SE_Z)
Trim_Fill

#Number of studies combined: k = 111 (with 10 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.1437 [-0.2440; -0.0433] -2.81  0.0050

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.30746122)
FisherZInv(0.03218604)
FisherZInv(0.59945190)
FisherZInv(0.30746122-0.1437)
FisherZInv(0.03218604-0.1437)
FisherZInv(0.59945190-0.1437)

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


###### Run the model subset by condition uncontrolled removed and color class ####
#load the model
load("R Files/rcondcolor_Z.RDATA")

rcondcolor.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, random = ~animal + Authors + us(SE_Z):units,
                    data=condcon, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
rcondcolor.Z$DIC
#-111.1669

#Save the model
#save(rcondcolor.Z, file = "rcondcolor_Z.RDATA")

#### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rcondcolor.Z$VCV/apply(rcondcolor.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.04302625      0.12503543      0.81891425      0.01302407

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
# 12.50354    

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 4.302625

## total heterogeneity percent
(I2s*100)+(I2u*100)
#16.80617 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2376015

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rcondcolor.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.31436643 -0.07767286  0.73040592  

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-condcon$Weight
MR<-condcon$Fisher_Z-pred_matrix[1:101]
zMR<-MR*Precision
condcon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=condcon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.79896    0.52216   1.530   0.1292   
#Precision   -0.19519    0.07074  -2.759   0.0069 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Non-significant Egger's Regression = trim fill only for funnel plot
Trim_Fill <- meta::trimfill(MR, condcon$SE_Z)
Trim_Fill

#Number of studies combined: k = 108 (with 7 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.1215 [-0.2165; -0.0265] -2.51  0.0122

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.31436643)
FisherZInv(-0.07767286)
FisherZInv(0.73040592)

mycol <- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
cols <- c(rep("black", 101), rep(mycol, 7))
p <- c(rep(16, 101), rep(15, 7))
FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = p,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

#Export 8x8


###### Run the model subset by social rank and condition uncontrolled removed ####
randccon <- rankcon[rankcon$Condition != "None measured",]

#load the model
load("R Files/rsandc_Z.RDATA")

rsandc.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                      data=randccon, pedigree = tree, 
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex)
rsandc.Z$DIC
#-105.4071

#Save the model
#save(rsandc.Z, file = "rsandc_Z.RDATA")

#### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rsandc.Z$VCV/apply(rsandc.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.007797701     0.017569402     0.967664597     0.006968300

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
# 1.75694   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.7797701

## total heterogeneity percent
(I2s*100)+(I2u*100)
#2.53671  

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2411506

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rsandc.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.23269131 0.07336586 0.38822407 

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-randccon$Weight
MR<-randccon$Fisher_Z-pred_matrix[1:78]
zMR<-MR*Precision
randccon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=randccon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  1.44243    0.44360   3.252  0.00171 ** 
#Precision   -0.21982    0.05815  -3.780  0.00031 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Trim_Fill <- meta::trimfill(MR, randccon$SE_Z)
Trim_Fill

#Number of studies combined: k = 87 (with 9 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.0637 [-0.1396; 0.0122] -1.64  0.1001

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.23269131)
FisherZInv(0.07336586)
FisherZInv(0.38822407)
FisherZInv(0.23269131-0.0637)
FisherZInv(0.07336586-0.0637)
FisherZInv(0.38822407-0.0637)

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


###### Run the model subset by social rank and condition uncontrolled removed and color####
#load the model
load("R Files/rsandcolor_Z.RDATA")

rsandcolor.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, random = ~animal + Authors + us(SE_Z):units,
                     data=randccon, pedigree = tree, 
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)
rsandcolor.Z$DIC
#-102.8041

#Save the model
save(rsandcolor.Z, file = "rsandcolor_Z.RDATA")

#### calculate I^2 to quantify heterogeneity ####
# Proportion of variance explained by random factors
rand <- rsandcolor.Z$VCV/apply(rsandcolor.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors SE_Z:SE_Z.units           units 
#0.01091610      0.01795318      0.96417528      0.00695545

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
# 1.795318   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.09161 

## total heterogeneity percent
(I2s*100)+(I2u*100)
#2.886927  

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.3047084 

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rsandcolor.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.239133149 0.007324507 0.481847788 

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-randccon$Weight
MR<-randccon$Fisher_Z-pred_matrix[1:78]
zMR<-MR*Precision
randccon[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=randccon)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  1.26185    0.43412   2.907 0.004784 ** 
#Precision   -0.19502    0.05691  -3.427 0.000988 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Trim_Fill <- meta::trimfill(MR, randccon$SE_Z)
Trim_Fill

#Number of studies combined: k = 86 (with 8 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.0637 [-0.1396; 0.0122] -1.64  0.1001

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.239133149)
FisherZInv(0.007324507)
FisherZInv(0.481847788)
FisherZInv(0.239133149-0.0637)
FisherZInv(0.007324507-0.0637)
FisherZInv(0.481847788-0.0637)

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


###### Run the model subset by age uncontrolled removed ####
#Subset data by removing papers that did not control for social rank
agecon <- metadatas[metadatas$Age_Controlled != "uncontrolled",]

#load the model
load("R Files/rage_Z.RDATA")

rage.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                      data=agecon, pedigree = tree, 
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex)
rage.Z$DIC
#13.23489

#Save the model
#save(rage.Z, file = "rage_Z.RDATA")

#### calculate I^2 to quantify heterogeneity ####
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

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rage.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2598076 -0.4301027  1.0783987  

#### publication bias ####
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


###### Run the model subset by age uncontrolled removed and color class####

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

#### calculate I^2 to quantify heterogeneity ####
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

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(ragec.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2891451 -1.2252288  1.8440485  

#### publication bias ####
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


###### Run the model subset by age, social rank, and condition uncontrolled removed and color class####
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

#### calculate I^2 to quantify heterogeneity ####
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

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(rankagecond.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.1769460 -0.1807661  0.5232827  

#### publication bias ####
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


################ Mixed Effects Model and Tests for the Best Model ####
###### Run the model with Class only ####
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
#-253.01

#Save the model for later
save(class.Z, file = "class_Z.RDATA")

###### Run the model with Class No Unknowns coded as Eumelanic ####
metadat <- metadata
metadat$Eu_Pheomelanin <- ifelse(metadata$Eu_Pheomelanin =="unknown", yes = "eumelanin", no = metadata$Eu_Pheomelanin)
unique(metadat$Eu_Pheomelanin)

randnounkn <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                       data=metadat, pedigree = tree, 
                       nitt = 5000000, thin = 1000, burnin = 2500000, 
                       prior = prior.ex)
randnounkn$DIC  
#-258.9389

pred_matrix <- predict(randnounkn, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2502108 -0.7282278  1.2379403

cnounk.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadat, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
cnounk.Z$DIC
#-259.4144

#Save the model for later
#save(cnounk.Z, file = "cnounk_Z.RDATA")

# Proportion of variance explained by random factors
rand <- cnounk.Z$VCV/apply(cnounk.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.009253155     0.021019517     0.968085932     0.001641395

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
#  2.101952

## species level heterogeneity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.9253155

## total heterogeneity
(I2s*100) + (I2u*100)
#3.027267 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2899397 

# Proportion of variance explained by random factors
Sol <- cnounk.Z$Sol/apply(cnounk.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(cnounk.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2583312 -0.7554615  1.2765707

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(cnounk.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.26136008 -0.02045448  0.55866470  

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadat$Weight
MR<-metadat$Fisher_Z-pred_matrix[1:158]
zMR<-MR*Precision
metadat[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadat)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 0.52561    0.37104   1.417   0.1586  
#Precision   -0.08306    0.04399  -1.888   0.0608 .
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#### Since the intercept of Egger's regression is not significant 
#### there is no publication bias. Trim and fill is only used for 
#### funnel plot creation

Trim_Fill <- trimfill(MR, metadat$SE_Z)
Trim_Fill

#Number of studies combined: k = 158 (with 0 added studies)

#                                       95%-CI    z p-value
#Random effects model 0.0133 [-0.0611; 0.0877] 0.35  0.7259


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


#### Medians and 95% Credible Intervals ####
emmeans(cnounk.Z, ~ Eu_Pheomelanin, data=metadat, level = 0.95)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.142    -0.130     0.408
#eumelanin       0.345     0.101     0.624
#pheomelanin     0.222    -0.116     0.616
#structural      0.246    -0.161     0.695
#Point estimate displayed: median 
#HPD interval probability: 0.95
Zcar <- 0.142; lcar <- -0.130; ucar <- 0.408
Zeu <- 0.345;  leu <-   0.101; ueu <- 0.624
Zph <- 0.222;  lph <-  -0.116; uph <- 0.616
Zst <- 0.246;  lst <-  -0.161; ust <- 0.695


#Probability above 0
p_significance(cnounk.Z$Sol, threshold = 0)
#Parameter                 |   ps
#Eu_Pheomelanincarotenoid  | 0.89
#Eu_Pheomelanineumelanin   | 0.99
#Eu_Pheomelaninpheomelanin | 0.91
#Eu_Pheomelaninstructural  | 0.88



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
###### Run the model with Age Controlled only ####
#Run the model
acon.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
acon.Z$DIC
#-256.5689


###### Run the model with Plasticity only ####
#Run the model
plasticity.Z <- MCMCglmm(Fisher_Z ~ Plasticity - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
plasticity.Z$DIC
#-261.1454

#Save the model for later
#save(plasticity.Z, file = "plasticity_Z.RDATA")

###### Run the model with Time Lag Analysis ####
#Run the model with linear term
year.Z <- MCMCglmm(Fisher_Z ~ as.numeric(Publication.Year) - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadatas, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
year.Z$DIC
#-263.6191

pred_matrix <- predict(year.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#0.25081365 0.01858821 0.46312191 

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
#52.36772

#Save the model for later
save(year2.Z, file = "yearquad_Z.RDATA")

pred_matrix <- predict(year2.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#0.265441702 -0.002931191  0.590576417

###### Run the model with Sex only ####
#Run the model
sex.Z <- MCMCglmm(Fisher_Z ~ Sex - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
sex.Z$DIC
#-260.0345

#Save the model for later
save(sex.Z, file = "sex_Z.RDATA")

###### Run the model with Vert/Invert only ####
#Run the model
vert.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
vert.Z$DIC
#-260.7969

#Save the model for later
#save(vert.Z, file = "vert_Z.RDATA")

###### Run the model with Location only ####
#Run the model
loc.Z <- MCMCglmm(Fisher_Z ~ Location - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
loc.Z$DIC
#-260.2084

#Save the model for later
#save(loc.Z, file = "loc_Z.RDATA")

###### Run the model with Season only ####
#Run the model
season.Z <- MCMCglmm(Fisher_Z ~ Season - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
season.Z$DIC
#-253.636

###### Run the model with Life Stage only ####
age.Z <- MCMCglmm(Fisher_Z ~ Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
age.Z$DIC
#-259.9691

###### Run the model with Geography only ####
geo.Z <- MCMCglmm(Fisher_Z ~ Geographic - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
geo.Z$DIC
#-261.2155

#save model
#save(geo.Z, file = "geo_Z.RDATA")

###### Run the model with Social Rank only ####
social.Z <- MCMCglmm(Fisher_Z ~ Social_Rank_Controlled - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
social.Z$DIC
#-257.1873

#save model
#save(social.Z, file = "social_Z.RDATA")

###### Run the model with Obs vs Exp only ####
obsexp.Z <- MCMCglmm(Fisher_Z ~ Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
obsexp.Z$DIC
#-261.2969

#save model
#save(obsexp.Z, file = "obsexp_Z.RDATA")

###### Run the model with Condition only ####
condition.Z <- MCMCglmm(Fisher_Z ~ Condition - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
condition.Z$DIC
#-259.8766

#save model
#save(condition.Z, file = "condition_Z.RDATA")

###### Run the model with Aggression Measures only ####
aggr.Z <- MCMCglmm(Fisher_Z ~ Aggression - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
aggr.Z$DIC
#-262.4415

#Save the model
#save(aggr.Z, file = "agg_Z.RDATA")

###### Run the model with Class and Vert/Invert ####
#Run the model
cv.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree, 
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)
cv.Z$DIC
#-258.4293

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
#-257.9002

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
#-251.133

###### Run the model with Class and Life Stage ####
#Run the model
ca.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Age - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
ca.Z$DIC
#-257.8608

###### Run the model with Plasticity and Life Stage ####
#Run the model
pa.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
pa.Z$DIC
#-259.9565

###### Run the model with Sex and Life Stage ####
#Run the model
sa.Z <- MCMCglmm(Fisher_Z ~ Sex + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
sa.Z$DIC
#-259.0631

###### Run the model with Vert/Invert and Life Stage ####
#Run the model
verta.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
verta.Z$DIC
#-259.2963

###### Run the model with Location and Life Stage ####
#Run the model
la.Z <- MCMCglmm(Fisher_Z ~ Location + Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
la.Z$DIC
#-259.6543

###### Run the model with Season and Life Stage ####
#Run the model
seaa.Z <- MCMCglmm(Fisher_Z ~ Season + Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
seaa.Z$DIC
#-254.3463

###### Run the model with Class and Plasticity ####
#Run the model
pc.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
pc.Z$DIC
#-257.7965

###### Run the model with Class and Sex ####
#Run the model
cs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
cs.Z$DIC
#-257.1659

###### Run the model with Plasticity and Sex ####
#Run the model
ps.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
ps.Z$DIC
#-259.4866

###### Run the model with Plasticity and Vert/Invert ####
#Run the model
pv.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
pv.Z$DIC
#-260.3315

###### Run the model with Plasticity and Location ####
#Run the model
ploc.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Location - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
ploc.Z$DIC
#-261.0393

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
#-252.1122

###### Run the model with Sex and Vert/Invert ####
#Run the model
svert.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
svert.Z$DIC
#-259.3059

###### Run the model with Sex and Location ####
#Run the model
seloc.Z <- MCMCglmm(Fisher_Z ~ Sex + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
seloc.Z$DIC
#-259.0473

###### Run the model with Sex and Season ####
#Run the model
ss.Z <- MCMCglmm(Fisher_Z ~ Sex + Season - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
ss.Z$DIC
#-252.897

###### Run the model with Vert/Invert and Location ####
#Run the model
vl.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
vl.Z$DIC
#-259.4229

###### Run the model with Vert/Invert and Season ####
#Run the model
vsea.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Season - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
vsea.Z$DIC
#-253.587

###### Run the model with Location and Season ####
#Run the model
sealoc.Z <- MCMCglmm(Fisher_Z ~ Season + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
sealoc.Z$DIC
#-254.2403


###### Run the model with Color and Geography ####
cg.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Geographic - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
cg.Z$DIC
#-257.8996

###### Run the model with Color and Social Rank ####
csoc.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Social_Rank_Controlled - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
csoc.Z$DIC
#-256.1298

###### Run the model with Color and Agg Measure ####
cagg.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Aggression - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
cagg.Z$DIC
#-259.6632

#Save the model for later
#save(cagg.Z, file = "cagg_Z.RDATA")

###### Run the model with Color and Condition ####
ccond.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Condition - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
ccond.Z$DIC
#-258.124
###### Run the model with Color and Obs vs Exp ####
cobs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Obs_vs_Exp - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
cobs.Z$DIC
#-258.1487

###### Run the model with Agg Measures and Plasticity ####
ap.Z <- MCMCglmm(Fisher_Z ~ Aggression + Plasticity - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

ap.Z$DIC
#-262.1645

#Save the model
save(ap.Z, file = "ap_Z.RDATA")

###### Run the model with Agg Measures and Sex ####
as.Z <- MCMCglmm(Fisher_Z ~ Aggression + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

as.Z$DIC
#-260.9065

###### Run the model with Agg Measures and Vert/Invert ####
av.Z <- MCMCglmm(Fisher_Z ~ Aggression + Vert_Invert - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

av.Z$DIC
#-261.3037

###### Run the model with Agg Measures and Location ####
al.Z <- MCMCglmm(Fisher_Z ~ Aggression + Location - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

al.Z$DIC
#-261.519

###### Run the model with Agg Measures and Season ####
asea.Z <- MCMCglmm(Fisher_Z ~ Aggression + Season - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

asea.Z$DIC
#-254.7479

###### Run the model with Agg Measures and Life Stage ####
aage.Z <- MCMCglmm(Fisher_Z ~ Aggression + Age - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

aage.Z$DIC
#-261.1512

###### Run the model with Agg Measures and Geography ####
aggeo.Z <- MCMCglmm(Fisher_Z ~ Aggression + Geographic - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

aggeo.Z$DIC
#-262.6289

#Save the model for later
save(aggeo.Z, file = "aggeo_Z.RDATA")

###### Run the model with Agg Measures and Social Rank ####
agsoc.Z <- MCMCglmm(Fisher_Z ~ Aggression + Social_Rank_Controlled - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

agsoc.Z$DIC
#-259.3838

###### Run the model with Agg Measures and Obs vs Exp ####
agob.Z <- MCMCglmm(Fisher_Z ~ Aggression + Obs_vs_Exp - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

agob.Z$DIC
#-262.0606

#Save the model for later
save(agob.Z, file = "agob_Z.RDATA")
###### Run the model with Agg Measures and Condition ####
agcond.Z <- MCMCglmm(Fisher_Z ~ Aggression + Condition - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

agcond.Z$DIC
#-259.8311

###### Run the model with Plasticity and Geography ####
pg.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Geographic - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree,
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)

pg.Z$DIC
#-261.0748

###### Run the model with Plasticity and Social Rank ####
psoc.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Social_Rank_Controlled - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

psoc.Z$DIC
#-257.6511

###### Run the model with Plasticity and Obs vs Exp ####
pob.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Obs_vs_Exp - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

pob.Z$DIC
#-261.1292

###### Run the model with Plasticity and Condition ####
pcond.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Condition - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree,
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)

pcond.Z$DIC
#-260.0619

###### Run the model with Sex and Geography ####
sg.Z <- MCMCglmm(Fisher_Z ~ Sex + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

sg.Z$DIC
#-260.1671

###### Run the model with Sex and Social Rank ####
ssoc.Z <- MCMCglmm(Fisher_Z ~ Sex + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

ssoc.Z$DIC
#-257.6619

###### Run the model with Sex and Obs vs Exp ####
sob.Z <- MCMCglmm(Fisher_Z ~ Sex + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

sob.Z$DIC
#-259.3934

###### Run the model with Sex and Condition ####
scond.Z <- MCMCglmm(Fisher_Z ~ Sex + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

scond.Z$DIC
#-258.7696

###### Run the model with Vert_Invert and Geography ####
vg.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

vg.Z$DIC
#-260.7141

###### Run the model with Vert_Invert and Social Rank ####
vsoc.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

vsoc.Z$DIC
#-256.7141

###### Run the model with Vert_Invert and Obs vs Exp ####
vob.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

vob.Z$DIC
#-260.4961

###### Run the model with Vert_Invert and Condition ####
vcond.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

vcond.Z$DIC
#-259.4847

###### Run the model with Location and Geography ####
lg.Z <- MCMCglmm(Fisher_Z ~ Location + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

lg.Z$DIC
#-260.1976

###### Run the model with Location and Social Rank ####
lsoc.Z <- MCMCglmm(Fisher_Z ~ Location + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

lsoc.Z$DIC
#-256.6105

###### Run the model with Location and Obs vs Exp ####
lob.Z <- MCMCglmm(Fisher_Z ~ Location + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

lob.Z$DIC
#-261.1541

###### Run the model with Location and Condition ####
lcond.Z <- MCMCglmm(Fisher_Z ~ Location + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

lcond.Z$DIC
#-258.8596

###### Run the model with Season and Geography ####
seag.Z <- MCMCglmm(Fisher_Z ~ Season + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

seag.Z$DIC
#-253.9413

###### Run the model with Season and Social Rank ####
seasoc.Z <- MCMCglmm(Fisher_Z ~ Season + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

seasoc.Z$DIC
#-251.773

###### Run the model with Season and Obs vs Exp ####
seaob.Z <- MCMCglmm(Fisher_Z ~ Season + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

seaob.Z$DIC
#-253.9706

###### Run the model with Season and Condition ####
seacond.Z <- MCMCglmm(Fisher_Z ~ Season + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

seacond.Z$DIC
#-253.2819

###### Run the model with Life Stage and Geography ####
ag.Z <- MCMCglmm(Fisher_Z ~ Age + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

ag.Z$DIC
#-259.8103

###### Run the model with Life Stage and Social Rank ####
asoc.Z <- MCMCglmm(Fisher_Z ~ Age + Social_Rank_Controlled - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree,
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)

asoc.Z$DIC
#-256.9167

###### Run the model with Life Stage and Obs vs Exp ####
aob.Z <- MCMCglmm(Fisher_Z ~ Age + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

aob.Z$DIC
#-260.1696

###### Run the model with Life Stage and Condition ####
acond.Z <- MCMCglmm(Fisher_Z ~ Age + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

acond.Z$DIC
#-259.29

###### Run the model with Geography and Social Rank ####
gsoc.Z <- MCMCglmm(Fisher_Z ~ Social_Rank_Controlled + Geographic - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)

gsoc.Z$DIC
#-258.682

###### Run the model with Geography and Obs vs Exp ####
gob.Z <- MCMCglmm(Fisher_Z ~ Geographic + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

gob.Z$DIC
#-258.4777

###### Run the model with Geography and Condition ####
gcond.Z <- MCMCglmm(Fisher_Z ~ Geographic + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

gcond.Z$DIC
#-260.5629

###### Run the model with Social Rank and Obs vs Exp ####
socob.Z <- MCMCglmm(Fisher_Z ~ Social_Rank_Controlled + Obs_vs_Exp - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)

socob.Z$DIC
#-257.7749

###### Run the model with Social Rank and Condition ####
soccond.Z <- MCMCglmm(Fisher_Z ~ Social_Rank_Controlled + Condition - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)

soccond.Z$DIC
#-256.129

###### Run the model with Obs_vs_Exp and Condition ####
obscond.Z <- MCMCglmm(Fisher_Z ~ Obs_vs_Exp + Condition - 1, 
                      random = ~animal + Authors + us(SE_Z):units, 
                      data=metadata, pedigree = tree,
                      nitt = 5000000, thin = 1000, burnin = 2500000, 
                      prior = prior.ex2)

obscond.Z$DIC
#-259.735

###### Run the model with Age controlled and Color ####
agc.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Eu_Pheomelanin - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
agc.Z$DIC
#-252.7446

###### Run the model with Age controlled and Agg Measure ####
agagg.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Aggression - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
agagg.Z$DIC
#-256.9074

###### Run the model with Age controlled and Plasticity ####
agp.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Plasticity - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
agp.Z$DIC
#-255.3949

###### Run the model with Age controlled and Sex ####
ags.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Sex - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
ags.Z$DIC
#-254.9366

###### Run the model with Age controlled and Vert_Invert ####
agv.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Vert_Invert - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
agv.Z$DIC
#-256.3538

###### Run the model with Age controlled and Location ####
agl.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Location - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
agl.Z$DIC
#-255.7213

###### Run the model with Age controlled and Season ####
agsea.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Season - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree,
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
agsea.Z$DIC
#-248.6002

###### Run the model with Age controlled and Life Stage ####
agage.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
agage.Z$DIC
#-255.4617

###### Run the model with Age controlled and Geography ####
aggeo.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Geographic - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
aggeo.Z$DIC
#-255.1007

###### Run the model with Age controlled and Social Rank ####
agsoc.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Social_Rank_Controlled - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
agsoc.Z$DIC
#-252.4422

###### Run the model with Age controlled and Obs vs. Exp ####
agobs.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Obs_vs_Exp - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
agobs.Z$DIC
#-255.7944

###### Run the model with Age controlled and Condition ####
agcon.Z <- MCMCglmm(Fisher_Z ~ Age_Controlled + Obs_vs_Exp - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
agcon.Z$DIC
#-256.1426

###### Mixed Effects Model Plasticity, Sex, and Plasticity x Sex ####
#Run the model
pxs.Z <- MCMCglmm(Fisher_Z ~ Sex + Plasticity + Sex*Plasticity - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pxs.Z$DIC
#-253.2273

###### Mixed Effects Model Class, Plasticity, and Class x Plasticity ####
#Run the model
ebyp.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity 
                   + Eu_Pheomelanin*Plasticity - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 2000000, thin = 1000, burnin = 1000000, 
                   prior = prior.ex2)
ebyp.Z$DIC
#-245.0806

###### Mixed Effects Model Class, Sex, and Class x Sex ####
#Run the model
ebys.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Sex*Eu_Pheomelanin - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 2000000, thin = 1000, burnin = 1000000, 
                   prior = prior.ex2)
ebys.Z$DIC
#-246.4917


###### COMBOS OF THREES ####
###### Run the model Class, Plasticity, and Sex ####
#Run the model
cps.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
cps.Z$DIC
#-257.1607

###### Run the model Class, Plasticity, and Vert/Invert ####
#Run the model
cpv.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity + Vert_Invert -1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 2000000, thin = 1000, burnin = 1000000, 
                   prior = prior.ex2)
cpv.Z$DIC
#-257.7878

###### Run the model Class, Sex, and Season ####
#Run the model
cssea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Season - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cssea.Z$DIC
#-250.2955

###### Run the model Class, Plasticity, and Location ####
#Run the model
cpl.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity + Location -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cpl.Z$DIC
#-258.4716

###### Run the model Class, Plasticity, and Season ####
#Run the model
cpsea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cpsea.Z$DIC
#-250.0485

###### Run the model Class, Sex, and Vert/Invert ####
#Run the model
csv.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Vert_Invert -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
csv.Z$DIC
#-256.2235

###### Run the model Class, Sex, and Location ####
#Run the model
csl.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Location -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
csl.Z$DIC
#-257.2403

###### Run the model Class, Vert/Invert, and Location ####
#Run the model
cvl.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Location -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cvl.Z$DIC
#-257.8201

###### Run the model Class, Vert/Invert, and Season ####
#Run the model
cvsea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cvsea.Z$DIC
#-250.1862

###### Run the model Class, Location, and Season ####
#Run the model
clsea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Location + Season -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
clsea.Z$DIC
#-251.5832

###### Run the model Class, Plasticity, and Age ####
#Run the model
cpa.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity + Age -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
cpa.Z$DIC
#-257.6487

###### Run the model Class, Sex, and Age ####
#Run the model
csa.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
csa.Z$DIC
#-255.9653

###### Run the model Class, Vert/Invert, and Age ####
#Run the model
cva.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cva.Z$DIC
#-257.9217

###### Run the model Class, Location, and Age ####
#Run the model
cla.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Location + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cla.Z$DIC
#-257.721

###### Run the model Class, Season, and Age ####
#Run the model
cseaa.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Season + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cseaa.Z$DIC
#-251.9882

###### Run the model Plasticity, Sex, and Age ####
#Run the model
psa.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex + Age -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
psa.Z$DIC
#-256.7983

###### Run the model Plasticity, Vert/Invert, and Age ####
#Run the model
pva.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pva.Z$DIC
#-258.5672

###### Run the model Plasticity, Location, and Age ####
#Run the model
pla.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Location + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
pla.Z$DIC
#-259.2949

###### Run the model Plasticity, Season, and Age ####
#Run the model
pseaa.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Season + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pseaa.Z$DIC
#-252.482

###### Run the model Plasticity, Sex, and Vert/Invert ####
#Run the model
psv.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex + Vert_Invert -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
psv.Z$DIC
#-259.1313

###### Run the model Plasticity, Sex, and Location ####
#Run the model
psl.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex + Location -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
psl.Z$DIC
#-260.092

###### Run the model Plasticity, Sex, and Season ####
#Run the model
pssea.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pssea.Z$DIC
#-252.2781

###### Run the model Plasticity, Vert/Invert, and Location ####
#Run the model
pvl.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert + Location -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
pvl.Z$DIC
#-260.258

###### Run the model Plasticity, Vert/Invert, and Season ####
#Run the model
pvsea.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pvsea.Z$DIC
#-252.3633

###### Run the model Plasticity, Location, and Season ####
#Run the model
plsea.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Location + Season -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
plsea.Z$DIC
#-253.9275

###### Run the model Sex, Vert/Invert, and Location ####
#Run the model
svl.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert + Location -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
svl.Z$DIC
#-258.5529

###### Run the model Sex, Vert/Invert, and Season ####
#Run the model
svsea.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
svsea.Z$DIC
#-252.0772

###### Run the model Sex, Location, and Season ####
#Run the model
slsea.Z <- MCMCglmm(Fisher_Z ~ Sex + Location + Season -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
slsea.Z$DIC
#-252.9407

###### Run the model Sex, Vert/Invert, and Age ####
#Run the model
sva.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert + Age -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
sva.Z$DIC
#-256.9187

###### Run the model Sex, Location, and Age ####
#Run the model
sla.Z <- MCMCglmm(Fisher_Z ~ Sex + Location + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
sla.Z$DIC
#-258.4774

###### Run the model Sex, Season, and Age ####
#Run the model
sseaa.Z <- MCMCglmm(Fisher_Z ~ Sex + Season + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
sseaa.Z$DIC
#-252.3519

###### Run the model Vert/Invert, Location, and Age ####
#Run the model
vla.Z <- MCMCglmm(Fisher_Z ~ Location + Vert_Invert + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
vla.Z$DIC
#-258.7068

###### Run the model Season, Vert/Invert, and Age ####
#Run the model
seava.Z <- MCMCglmm(Fisher_Z ~ Season + Vert_Invert + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
seava.Z$DIC
#-254.1434

###### Run the model Location, Season, and Age ####
#Run the model
lseaa.Z <- MCMCglmm(Fisher_Z ~ Location + Season + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
lseaa.Z$DIC
#-254.265

###### Run the model Vert/Invert, Location, and Season ####
#Run the model
vlsea.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Location + Season -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
vlsea.Z$DIC
#-253.476

###### Run the model Class, Plasticity, Sex, and Vert/Invert ####
#Run the model
vert.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity + Vert_Invert -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
vert.Z$DIC
#-255.7414

###### Run the model Class, Plasticity, Sex and Location ####
#Run the model
loc.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity + Location -1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 2000000, thin = 1000, burnin = 1000000, 
                   prior = prior.ex2)
loc.Z$DIC
#-257.3939

###### Mixed Effects Model Class, Plasticity, Sex, and Plasticity x Sex ####
#Run the model
pxs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity 
                  + Sex*Plasticity - 1, 
                random = ~animal + Authors + us(SE_Z):units, 
                data=metadata, pedigree = tree, 
                nitt = 2000000, thin = 1000, burnin = 1000000, 
                prior = prior.ex2)
pxs.Z$DIC
#-256.1967

###### Mixed Effects Model Class, Plasticity, Sex, and Class x Plasticity ####
#Run the model
ebyp.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity 
                   + Eu_Pheomelanin*Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 2000000, thin = 1000, burnin = 1000000, 
                 prior = prior.ex2)
ebyp.Z$DIC
#-248.9188

###### Mixed Effects Model Class, Plasticity, Sex, and Class x Sex ####
#Run the model
ebys.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity
                 + Sex*Eu_Pheomelanin - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 2000000, thin = 1000, burnin = 1000000, 
                 prior = prior.ex2)
ebys.Z$DIC
#-252.0912


################## "BEST" MODEL (class) ##############
load("R Files/class_Z.RDATA")

color.Z <- class.Z

summary(color.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -258.4055 

#Location effects: Fisher_Z ~ Eu_Pheomelanin - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#Eu_Pheomelanincarotenoid    0.14136 -0.09293  0.43880     2258 0.2152  
#Eu_Pheomelanineumelanin     0.32579  0.05894  0.59783     2651 0.0176 *
#Eu_Pheomelaninpheomelanin   0.21991 -0.13056  0.61024     2333 0.2040  
#Eu_Pheomelaninstructural    0.23796 -0.21834  0.67230     2500 0.2840  
#Eu_Pheomelaninunknown       0.43468  0.07568  0.83967     2500 0.0392 *
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

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.008658133     0.020632088     0.968998857     0.001710922

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
# 2.063209

## species level heterogeneity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.8658133 

## total heterogeneity
(I2s*100) + (I2u*100)
#2.929022 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2792843 

# Proportion of variance explained by random factors
Sol <- color.Z$Sol/apply(color.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(color.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2551818 -0.7526698  1.2739840

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(color.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.25697879 -0.05568362  0.57507478  

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$Weight
MR<-metadata$Fisher_Z-pred_matrix[1:158]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 0.42765    0.36746   1.164    0.246
#Precision   -0.06134    0.04395  -1.396    0.165
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#### Since the intercept of Egger's regression is not significant 
#### there is no publication bias. Trim and fill is only used for 
#### funnel plot creation

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 158 (with 0 added studies)

#                                       95%-CI    z p-value
#Random effects model 0.0185 [-0.0538; 0.0907] 0.50  0.6158


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


#### Medians and 95% Credible Intervals ####
emmeans(color.Z, ~ Eu_Pheomelanin, data=metadata, level = 0.95)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.1206   -0.1645     0.388
#eumelanin       0.3282    0.0559     0.608
#pheomelanin     0.2166   -0.1038     0.624
#pteridine      -0.0944   -0.8887     0.782
#structural      0.2365   -0.2244     0.660
#unknown         0.4328    0.0613     0.827
#Point estimate displayed: median 
#HPD interval probability: 0.95
Zcar <- 0.137; lcar <- -0.0929; ucar <- 0.439
Zeu <- 0.315;  leu <-   0.0589; ueu <- 0.598
Zph <- 0.210;  lph <-  -0.1306; uph <- 0.610
Zst <- 0.236;  lst <-  -0.2183; ust <- 0.672
Zun <- 0.441;  lun <-   0.0757; uun <- 0.840

#Probability above 0
p_significance(color.Z$Sol, threshold = 0)
#Parameter                 |   ps
#Eu_Pheomelanincarotenoid  | 0.89
#Eu_Pheomelanineumelanin   | 0.99
#Eu_Pheomelaninpheomelanin | 0.90
#Eu_Pheomelaninstructural  | 0.86
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
################## "BEST" MODEL (aggression measures) ##############
load("R Files/agg_Z.RDATA")

summary(aggr.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC:  -262.4415
#Location effects: Fisher_Z ~ Aggression - 1 
#                   post.mean l-95% CI u-95% CI eff.samp  pMCMC   
#AggressionDirect     0.31278  0.12481  0.54175     2500 0.0088 **
#AggressionIndirect   0.18700 -0.01168  0.41621     2500 0.0568 . 
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

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.005254161     0.023307691     0.969755897     0.001682251

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
#  2.330769

## species level heterogeneity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.5254161

## total heterogeneity
(I2s*100) + (I2u*100)
#2.856185 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.1737251 

# Proportion of variance explained by random factors
Sol <- aggr.Z$Sol/apply(aggr.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(aggr.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2494707 -0.7141408  1.2149455

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(aggr.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.25228065 0.05915659 0.48136754  

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$Weight
MR<-metadata$Fisher_Z-pred_matrix[1:158]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 0.46103    0.36651   1.258     0.21
#Precision   -0.06981    0.04345  -1.607     0.11
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#### Since the intercept of Egger's regression is not significant 
#### there is no publication bias. Trim and fill is only used for 
#### funnel plot creation

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 158 (with 0 added studies)

#                                       95%-CI    z p-value
#Random effects model 0.0156 [-0.0570; 0.0882] 0.42  0.6737

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
                  
                  
#### Medians and 95% Credible Intervals ####
emmeans(aggr.Z, ~ Aggression, data=metadata, level = 0.95)
#Aggression emmean lower.HPD upper.HPD
#Direct      0.309    0.1248     0.542
#Indirect    0.181   -0.0117     0.416                  

Zdir <- 0.309; ldir <- 0.1248; udir <- 0.542
Zindir <- 0.181; lindir <- -0.0117; uindir <- 0.416

#Probability above 0
p_significance(aggr.Z$Sol, threshold = 0)
#Parameter          |   ps
#AggressionDirect   | 1.00
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
################## "BEST" MODEL (plasticity) ##############
load("R Files/plasticity_Z.RDATA")

plas.Z <- plasticity.Z

summary(plas.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -238.9878 

#Location effects: Fisher_Z ~ Plasticity - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#PlasticityNo        0.29045  0.03976  0.53710     2326 0.0216 *
#PlasticityPlastic   0.20638 -0.06833  0.49038     2500 0.1272  
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(plas.Z$Sol)
autocorr(plas.Z$Sol)
autocorr(plas.Z$VCV)

xsim <- simulate(plas.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- plas.Z$VCV/apply(plas.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

####calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.008941009     0.023053009     0.966301737     0.001704245

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

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
#  2.305301

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.8941009 

## total heterogeneity
(I2s*100) + (I2u*100)
#3.199402

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2653255

# Proportion of variance explained by random factors
Sol <- plas.Z$Sol/apply(plas.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(plas.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2650784 -0.7688504  1.3029157

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(plas.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#       fit        lwr        upr 
#0.26224001 0.00348823 0.52142155

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$Weight
MR<-metadata$Fisher_Z-pred_matrix[1:149]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Call:
#glm(formula = zMR ~ Precision, family = "gaussian", data = metadata)

#Deviance Residuals: 
#  Min      1Q    Median    3Q      Max  
#-6.2276  -1.4924  -0.5019   1.1344  13.4256   

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  1.00845    0.38510   2.619 0.009752 **
#Precision   -0.15645    0.04585  -3.412 0.000833 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 6.250279)

#Null deviance: 991.56  on 148  degrees of freedom
#Residual deviance: 918.79  on 147  degrees of freedom
#AIC: 699.89

#Number of Fisher Scoring iterations: 2

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 156 (with 7 added studies)

#                                         95%-CI     z  p-value
#Random effects model   -0.0545 [-0.1171; 0.0081] -1.71  0.0878

#Quantifying heterogeneity:
#tau^2 = 0.1214 [0.2164; 0.3776]; tau = 0.3484 [0.4651; 0.6145]
#I^2 = 89.0% [87.6%; 90.3%]; H = 3.02 [2.84; 3.21]

#Test of heterogeneity:
#     Q  d.f. p-value
#1414.02  155 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry

#### Adjusted mean and 95% confidence interval #####
#              fit               lwr              upr 
#0.26224001-0.0545 0.00348823-0.0545 0.52142155-0.0545
#          0.20774       -0.05101177         0.4669216
0.26224001-0.0545
0.00348823-0.0545
0.52142155-0.0545

#Means and 95%CI in rho rather than Fisher Z
FisherZInv(0.26224001)
FisherZInv(0.00348823)
FisherZInv(0.52142155)
FisherZInv(0.26224001-0.0545)
FisherZInv(0.00348823-0.0545)
FisherZInv(0.52142155-0.0545)



cols <- rep("black",157)
shapes <- rep(0,157)
for (i in 1:157){
  if (i >= 150){
    shapes[i] <- 8
  } else if (metadata$Plasticity[i] == "Plastic"){
    shapes[i] <- 15
  } else if (metadata$Plasticity[i] == "No"){
    shapes[i] <- 16
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Plastic", "Non-Plastic", "Trim and Fill Points")
points <- c(15,16,8)
ys <- c(20.1,19.3,18.5)
for(i in 1:3){
  points(x=-2.9, y=ys[i], pch=points[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Medians and 95% Credible Intervals ####
#Probability above 0
p_significance(plas.Z$Sol, threshold = 0)
#Parameter         |   ps
#PlasticityNo      | 0.99
#PlasticityPlastic | 0.94


emmeans(plas.Z, ~ Plasticity, data=metadata, level = 0.95)
#Plasticity emmean lower.HPD upper.HPD
#No          0.282    0.0398     0.537
#Plastic     0.209   -0.0683     0.490
Znpl <- 0.282; lnpl <- 0.0398; unpl <- 0.537
Zpl <- 0.209;  lpl <- -0.0683;  upl <- 0.490

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.3,0.7),axes=F,ann=F)
axis(1)
#Fisher Z
#Plastic
segments(lpl,0.6,upl,0.6);
points(Zpl,0.6,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl,0.4,unpl,0.4);
points(Znpl,0.4,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.3,0.7),axes=F,ann=F)
axis(1)
#Fisher Z
#Plastic
segments(lpl-0.0545,0.6,upl-0.0545,0.6);
points(Zpl-0.0545,0.6,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl-0.0545,0.4,unpl-0.0545,0.4);
points(Znpl-0.0545,0.4,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))

#Export 6x6
######## Rho Figure (Figure 2)#######
plot(NA,xlim=c(-1.1,1.1),ylim=c(-0.1,0.8),axes=F,ann=F)
axis(1)

#### Random Effects Model
#Mean No correction
segments(FisherZInv(0.04646076),0.7,FisherZInv(0.48347851),0.7);
points(FisherZInv(0.26002580),0.7,pch=16,col = "black",xpd=NA)

#### Random Effects Model
#Mean with correction
segments(FisherZInv(0.04646076-0.0255),0.6,FisherZInv(0.48347851-0.0255),0.6);
points(FisherZInv(0.26002580-0.0255),0.6,pch=15,col = "black",xpd=NA)

#### Mixed: Plasticity
#Plastic no correction
segments(FisherZInv(lpl),0.5,FisherZInv(upl),0.5);
points(FisherZInv(Zpl),0.5,pch=16,col = "black",xpd=NA)
#Non-Plastic no correction
segments(FisherZInv(lnpl),0.3,FisherZInv(unpl),0.3);
points(FisherZInv(Znpl),0.3,pch=16,col = "black",xpd=NA)
#Mean Plasticity Model no correction
segments(0.003,0.1,0.479,0.1);
points(0.256,0.1,pch=16,col = "black",xpd=NA)

#### Mixed: Plasticity
#Plastic with correction
segments(FisherZInv(lpl-0.0545),0.4,FisherZInv(upl-0.0545),0.4);
points(FisherZInv(Zpl-0.0545),0.4,pch=15,col = "black",xpd=NA)
#Non-Plastic with correction
segments(FisherZInv(lnpl-0.0545),0.2,FisherZInv(unpl-0.0545),0.2);
points(FisherZInv(Znpl-0.0545),0.2,pch=15,col = "black",xpd=NA)
#Mean Plasticity Model with correction
segments(-0.051,0,0.436,0);
points(0.205,0,pch=15,col = "black",xpd=NA)

#Add line at 0 and separate models
abline(h = 0.55, lty = 1)
abline(v = 0, lty = 2)

#Add axis labels
title(xlab = "Correlation Coefficient")
text(-1,0.65,"Random Effects Model", cex = 0.9, adj = c(0,0))
text(FisherZInv(0.26002580),0.7,"*", cex = 1, adj = c(0.5,0))
text(FisherZInv(0.26002580-0.0255),0.6,"*", cex = 1, adj = c(0.5,0))
text(-1,0.45,"Plastic", cex = 0.9, adj = c(0,0))
text(-1,0.25,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(FisherZInv(Znpl),0.3,"*", cex = 1, adj = c(0.4,0))
text(-1,0.05,"Overall Model", cex = 0.9, adj = c(0,0))
text(0.256,0.1,"*", cex = 1, adj = c(0.4,0))

#legend
words <- c("Before Correction", "After Correction")
Cols <- c("black", "black")
points <- c(16,15)
ys <- c(0.8, 0.75)
for(i in 1:2){
  points(x=0.7, y=ys[i], pch=points[i], col=Cols[i])
  text(x=0.7,y=ys[i], labels=words[i], pos=4, cex=.75, font = 2)
}

#Export 7x7

################## Fisher Z with Pub Bias correction subset combined ####
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
################## Fisher Z with Publication Bias correction all combined ####
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

################## Fisher Z for each Study #####
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

plot(NA,xlim=c(-4,4),ylim=c(0,330),axes=F,ann=F)
axis(1)
polygon(x = c(-4,-4,4,4), y = c(8,14,14,8), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-4,-4,4,4), y = c(16,40,40,16), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-4,-4,4,4), y = c(162,214,214,162), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-4,-4,4,4), y = c(238,320,320,238), col = alpha("gray", 0.8), density = NA)

abline(v=0)
for (i in 1:length(order$Study)){
  segments(order$lci[i],i*2,order$uci[i],i*2);
  points(order$Fisher_Z[i],i*2,pch=shapesa[i], col = colsa[i],xpd=NA)
  text(-4, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
}

#Random effects model before and then after
points(x=0.25088699, y = 330,pch = 17)
segments(0.01803682, 330, 0.45812805, 330)
points(x=0.25088699+0.0233, y= 325,pch = 18, cex=1.2)
segments(0.01803682+0.0233, 325, 0.45812805+0.0233, 325)

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
text(-4.075,328,cex=0.5, pos=4,"Random Effects Model", font=2)


#Export 15x10

################## Fisher Z by Phylogeny ####
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

