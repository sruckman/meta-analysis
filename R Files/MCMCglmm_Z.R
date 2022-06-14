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

#### Phylogeny ####
#First, we need to read in the tree and root it to convert it from an 
#ultrametric tree to a rooted tree for the analysis
#Our outgroup is the sea anemone, Phymactis clematis 
tree <- read.tree("Excel Sheets/list.nwk")
root(tree, outgroup = "Phymactis_clematis")
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
metadata <- read.csv("Excel Sheets/meta_complete_data.csv", header=T)

#Change Eu_Phae to represent all pigments
metadata$Eu_Pheomelanin <- ifelse(metadata$Eu_Pheomelanin=="N/A", metadata$Classification, metadata$Eu_Pheomelanin)
unique(metadata$Eu_Pheomelanin)

#we must change species to animal for the analysis to run
names(metadata)[5] <- "animal"
#Check to see if the labels on our tree match the animal column
tree$tip.label<-gsub("_"," ",tree$tip.label)
tips <- TipLabels(tree)
tips %in% metadata$animal

#store for later figures
metadatas <- metadata
#remove pterinidine point (only one)
metadata <- metadata[c(1:125, 127:150),]

################ Random Effects Model with strict prior ####
load("R Files/allRnd_Z.RDATA")

#priors:
prior.ex <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model:
allRnd.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + us(SE_Z):units,
                        data=metadata, pedigree = tree, 
                        nitt = 800000, thin = 100, burnin = 600000, 
                        prior = prior.ex)

#Save the model for later
#save(allRnd.Z, file = "allRnd_Z.RDATA")

#Summary of Results
summary(allRnd.Z)
#Iterations = 600001:799901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -238.297 

#G-structure:  ~animal
#       post.mean  l-95% CI u-95% CI eff.samp
#animal   0.02486 3.514e-09   0.1046     1439

#~Authors
#        post.mean  l-95% CI u-95% CI eff.samp
#Authors   0.08147  0.02262   0.1505     2000

#~us(SE_Z):units
#              post.mean l-95% CI u-95% CI eff.samp
#SE_Z:SE_Z.units   3.643    2.428    5.016     2000

#R-structure:  ~units
#          post.mean l-95% CI u-95% CI eff.samp
#units      0.006539   0.0014   0.0143     2000

#Location effects: Fisher_Z ~ 1  
#             post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)    0.26257  0.05691  0.48449     2000 0.029 *


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
#0.006815841     0.023035157     0.968337959     0.001811042

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
# 2.303516   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.6815841

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2152685 

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(allRnd.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.26256669 0.05691203 0.48449461  

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
#  Min        1Q    Median        3Q       Max  
#-6.5967  -1.4894  -0.3931   1.0565  13.5862

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.86894    0.38531   2.255  0.02560 * 
#Precision   -0.13223    0.04588  -2.882  0.00454 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 6.257106)

#Null deviance: 971.77  on 148  degrees of freedom
#Residual deviance: 919.79  on 147  degrees of freedom
#AIC: 700.05

#Number of Fisher Scoring iterations: 2

Trim_Fill <- meta::trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 153 (with 4 added studies)

#                                         95%-CI     z  p-value
#Random effects model   -0.0280 [-0.0862; 0.0302] -0.94  0.3458

#Quantifying heterogeneity:
#tau^2 = 0.0993 [0.1824; 0.3213]; tau = 0.3152 [0.4271; 0.5668]
#I^2 = 87.1% [85.3%; 88.7%]; H = 2.78 [2.61; 2.97]

#Test of heterogeneity:
#      Q  d.f. p-value
#1178.21  152 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry
cols <- rep(0,154)
for (i in 1:154){
  if (i > 149){
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

shapes <- rep(0,154)
for (i in 1:154){
  if (i >= 150){
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
for(i in 1:9){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 6x6


###### Run the model without species ####
priors <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))

sp.Z <- MCMCglmm(Fisher_Z ~ 1, 
                    random = ~ Authors + us(SE_Z):units, 
                    data=metadata, 
                    nitt = 600000, thin = 100, burnin = 400000, 
                    prior = priors)
sp.Z$DIC
#-235.9753

###### Run the model without Authors ####
au.Z <- MCMCglmm(Fisher_Z ~ 1, 
                 random = ~ animal + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = priors)
au.Z$DIC
#-221.2108

###### Run the model without Authors ####
weight.Z <- MCMCglmm(Fisher_Z ~ 1, 
                 random = ~ animal + Authors, 
                 data=metadata, pedigree = tree,
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = priors)
weight.Z$DIC
#-173.1499

################ Mixed Effects Model and Tests for the Best Model ####
###### Run the model with Class only ####
prior.ex2 <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
                           G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))

class.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree,
                     nitt = 600000, thin = 100, burnin = 400000, 
                     prior = prior.ex2)
class.Z$DIC
#-235.6355

#Compare to Random Model:
teststat <- (as.numeric(allRnd.Z$DIC)-as.numeric(class.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.2486678, Not significantly better

###### Run the model with Plasticity only ####
#Run the model
plasticity.Z <- MCMCglmm(Fisher_Z ~ Plasticity - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 600000, thin = 100, burnin = 400000, 
                         prior = prior.ex2)
plasticity.Z$DIC
#-236.9118

#Compare to Random Model:
teststat <- (as.numeric(allRnd.Z$DIC)-as.numeric(plasticity.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.4052804, Not significantly better

###### Run the model with Sex only ####
#Run the model
sex.Z <- MCMCglmm(Fisher_Z ~ Sex - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 600000, thin = 100, burnin = 400000, 
                         prior = prior.ex2)
sex.Z$DIC
#-237.1848

#Compare to Random Model:
teststat <- (as.numeric(allRnd.Z$DIC)-as.numeric(sex.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.4558261, Not significantly better

###### Run the model with Class and Plasticity ####
#Run the model
plastic.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2)
plastic.Z$DIC
#-234.8729

#Compare to Random Model:
teststat <- (as.numeric(allRnd.Z$DIC)-as.numeric(plastic.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.1907189, Not significantly better

#Compare to Class Only Model:
teststat <- (as.numeric(class.Z$DIC)-as.numeric(plastic.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.5369173, Not significantly better

#Compare to Plasticity Only Model:
teststat <- (as.numeric(plasticity.Z$DIC)-as.numeric(plastic.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.3126489, Not significantly better

###### Run the model with Plasticity and Sex ####
#Run the model
ps.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex - 1, 
                       random = ~animal + Authors + us(SE_Z):units, 
                       data=metadata, pedigree = tree, 
                       nitt = 600000, thin = 100, burnin = 400000, 
                       prior = prior.ex2)
ps.Z$DIC
#-238.3838

#Compare to Random Model:
teststat <- (as.numeric(allRnd.Z$DIC)-as.numeric(ps.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 1, Not significantly better

#Compare to Class Only Model:
teststat <- (as.numeric(sex.Z$DIC)-as.numeric(ps.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 1, Not significantly better

#Compare to Plasticity Only Model:
teststat <- (as.numeric(plasticity.Z$DIC)-as.numeric(ps.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 1, Not significantly better

###### Run the model with Class and Sex ####
#Run the model
cs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = prior.ex2)
cs.Z$DIC
#-234.2826

#Compare to Random Model:
teststat <- (as.numeric(allRnd.Z$DIC)-as.numeric(cs.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.1565517, Not significantly better

#Compare to Class Only Model:
teststat <- (as.numeric(class.Z$DIC)-as.numeric(cs.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.4108119, Not significantly better

#Compare to Sex Only Model:
teststat <- (as.numeric(sex.Z$DIC)-as.numeric(cs.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.2283518, Not significantly better

###### Run the model with Class, Plasticity, and Sex ####
#Run the model
cps.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = prior.ex2)
cps.Z$DIC
#-234.0214

#Compare to Random Model:
teststat <- (as.numeric(allRnd.Z$DIC)-as.numeric(cps.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.1437089, Not significantly better

#Compare to Class and Plastic Model:
teststat <- (as.numeric(plastic.Z$DIC)-as.numeric(cps.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.5140902, Not significantly better

#Compare to Plasticity and Sex Model:
teststat <- (as.numeric(ps.Z$DIC)-as.numeric(cps.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.1397066, Not significantly better

#Compare to Class and Sex Model:
teststat <- (as.numeric(cs.Z$DIC)-as.numeric(cps.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.7178457, Not significantly better

###### Run the model add Vert/Invert ####
#Run the model
vert.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity + Vert_Invert -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 600000, thin = 100, burnin = 400000, 
                  prior = prior.ex2)
vert.Z$DIC
#-233.4535

#Compare to Class, Plasticity, and Sex Model:
teststat <- (as.numeric(cps.Z$DIC)-as.numeric(vert.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.5941223, Not significantly better

###### Run the model add Location - Vert/Invert ####
#Run the model
loc.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity + Location -1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2)
loc.Z$DIC
#-234.7651

#Compare to Class, Plasticity, and Sex Model:
teststat <- (as.numeric(cps.Z$DIC)-as.numeric(loc.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 1, Not significantly better

###### Run the model with Season not Plasticity ####
#Run the model
sea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Season - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 600000, thin = 100, burnin = 400000, 
                  prior = prior.ex2)
sea.Z$DIC
#-225.1154

#Compare to Class and Sex Model:
teststat <- (as.numeric(cs.Z$DIC)-as.numeric(sea.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.03227908, Significantly worse

###### Mixed Effects Model with Plasticity x Sex ####
#Run the model
pxs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity 
                  + Sex*Plasticity - 1, 
                random = ~animal + Authors + us(SE_Z):units, 
                data=metadata, pedigree = tree, 
                nitt = 600000, thin = 100, burnin = 400000, 
                prior = prior.ex2)
pxs.Z$DIC
#-232.2638

#Compare to Class, Plasticity, and Sex Model:
teststat <- (as.numeric(cps.Z$DIC)-as.numeric(pxs.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.3485277, Not significantly better

###### Mixed Effects Model with Eu_Pheomelanin x Plasticity ####
#Run the model
ebyp.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity 
                   + Eu_Pheomelanin*Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = prior.ex2)
ebyp.Z$DIC
#-225.5206

#Compare to Class, Plasticity, and Sex Model:
teststat <- (as.numeric(cps.Z$DIC)-as.numeric(ebyp.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.03924045, Significantly worse than cps model

###### Mixed Effects Model with Eu_Pheomelanin x Sex ####
#Run the model
ebys.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity
                 + Sex*Eu_Pheomelanin - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = prior.ex2)
ebys.Z$DIC
#-230.5733

#Compare to Class, Plasticity, and Sex Model:
teststat <- (as.numeric(cps.Z$DIC)-as.numeric(ebys.Z$DIC))/-2
p.val <- pchisq(teststat, df = 1, lower.tail = TRUE)
1-p.val
#p-value = 0.1891704, Not significantly better

################## "BEST" MODEL (class and plasticity) ##############
load("R Files/mixed_Z.RDATA")

#prior with expanded parameters (Cuachy Distribution close to a Fisher Z)
prior.ex2<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model
mixed.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 800000, thin = 100, burnin = 600000, 
                   prior = prior.ex2)
#Save the model for later
#save(mixed.Z, file = "mixed_Z.RDATA")

summary(mixed.Z)
#Iterations = 600001:799901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -234.0923 

#Location effects: Fisher_Z ~ Eu_Pheomelanin + Plasticity - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#Eu_Pheomelanincarotenoid    0.15742 -0.18712  0.53343     2000 0.323  
#Eu_Pheomelanineumelanin     0.34961  0.03481  0.71041     2154 0.043 *
#Eu_Pheomelaninpheomelanin   0.22786 -0.17709  0.66700     2000 0.238  
#Eu_Pheomelaninstructural    0.34005 -0.08905  0.81139     1766 0.127  
#Eu_Pheomelaninunknown       0.41805 -0.03605  0.84654     2000 0.074 .
#PlasticityPlastic          -0.01255 -0.26089  0.22384     2000 0.940  
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(mixed.Z$Sol)
autocorr(mixed.Z$Sol)
autocorr(mixed.Z$VCV)

xsim <- simulate(mixed.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- mixed.Z$VCV/apply(mixed.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.015918122     0.023260701     0.958985695     0.001835482

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
#  2.32607

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.591812 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.3881115

# Proportion of variance explained by random factors
Sol <- mixed.Z$Sol/apply(mixed.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(mixed.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2681985 -0.8466088  1.3802318

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(mixed.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2662392 -0.1043692  0.6625217 

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
#  Min        1Q    Median      3Q       Max  
#-5.8672  -1.3691  -0.2501   1.1773  13.1018    

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.55101    0.37024   1.488    0.139  
#Precision   -0.07960    0.04408  -1.806    0.073 .
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 5.671266)

#Null deviance: 868.11  on 148  degrees of freedom
#Residual deviance: 849.27  on 147  degrees of freedom
#AIC: 688.17

#Number of Fisher Scoring iterations: 2

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 149 (with 0 added studies)

#                                         95%-CI     z  p-value
#Random effects model     0.0100 [-0.0413; 0.0613] 0.38  0.7018

#Quantifying heterogeneity:
#tau^2 = 0.0694 [0.1142; 0.2065]; tau = 0.2634 [0.3380; 0.4544]
#I^2 = 82.8% [80.2%; 85.1%]; H = 2.41 [2.25; 2.59]

#Test of heterogeneity:
#     Q  d.f. p-value
#862.07  148 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry

cols <- rep(0,150)
for (i in 1:150){
  if (i > 150){
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
  if (i >= 151){
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
#Export 6x6


#### Medians and 95% Credible Intervals ####
emmeans(mixed.Z, ~ Eu_Pheomelanin, data=metadata, level = 0.95)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.151  -0.16446     0.512
#eumelanin       0.336  -0.00644     0.701
#pheomelanin     0.216  -0.19159     0.684
#structural      0.330  -0.12731     0.799
#unknown         0.417  -0.07124     0.820
#Point estimate displayed: median 
#HPD interval probability: 0.95
Zcar <- 0.151; lcar <- -0.1645; ucar <- 0.512
Zeu <- 0.336;  leu <-  -0.0064; ueu <- 0.701
Zph <- 0.216;  lph <-  -0.1916; uph <- 0.684
Zst <- 0.330;  lst <-  -0.1273; ust <- 0.799
Zun <- 0.417;  lun <-  -0.0712; uun <- 0.820

#Probability above 0
p_significance(mixed.Z$Sol, threshold = 0)
#Parameter                 |   ps
#Eu_Pheomelanincarotenoid  | 0.84
#Eu_Pheomelanineumelanin   | 0.98
#Eu_Pheomelaninpheomelanin | 0.88
#Eu_Pheomelaninstructural  | 0.94
#Eu_Pheomelaninunknown     | 0.96
#PlasticityPlastic         | 0.53



emmeans(mixed.Z, ~ Plasticity, data=metadata, level = 0.95)
#Plasticity emmean lower.HPD upper.HPD
#No          0.291   -0.0320     0.604
#Plastic     0.284   -0.0631     0.674
Znpl <- 0.291; lnpl <- -0.0320; unpl <- 0.604
Zpl <- 0.284;  lpl <- -0.0631;  upl <- 0.674

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.3,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
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
text(-2.1,1.6,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Eumelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1.2,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Structural", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Unknown", cex = 0.9, adj = c(0,0))
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.3,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
#Carotenoid
segments(lcar+0.01,1.6,ucar+0.01,1.6);
points(Zcar+0.01,1.6,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu+0.01,1.4,ueu+0.01,1.4);
points(Zeu+0.01,1.4,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph+0.01,1.2,uph+0.01,1.2);
points(Zph+0.01,1.2,pch=16,col = "black",xpd=NA)
#Structural 
segments(lst+0.01,1,ust+0.01,1);
points(Zst+0.01,1,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun+0.01,0.8,uun+0.01,0.8);
points(Zun+0.01,0.8,pch=16,col = "black",xpd=NA)
#Plastic
segments(lpl+0.01,0.6,upl+0.01,0.6);
points(Zpl+0.01,0.6,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl+0.01,0.4,unpl+0.01,0.4);
points(Znpl+0.01,0.4,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Eumelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1.2,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Structural", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Unknown", cex = 0.9, adj = c(0,0))
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))

#Export 6x6
################## "BEST" MODEL (plasticity and sex) ##############
load("R Files/mixed1_Z.RDATA")

#prior with expanded parameters (Cuachy Distribution close to a Fisher Z)
prior.ex2<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model
mixed1.Z <- MCMCglmm(Fisher_Z ~ Sex + Plasticity - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 800000, thin = 100, burnin = 600000, 
                    prior = prior.ex2)
#Save the model for later
#save(mixed1.Z, file = "mixed1_Z.RDATA")

summary(mixed1.Z)
#Iterations = 600001:799901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -236.0665 

#Location effects: Fisher_Z ~ Sex + Plasticity 

#                  post.mean  l-95% CI  u-95% CI eff.samp pMCMC  
#Sexboth            0.337465 -0.008471  0.692250     2209 0.053 .
#Sexfemales         0.275037 -0.060811  0.641119     2100 0.089 .
#Sexmales           0.281392 -0.012476  0.599365     2180 0.061 .
#PlasticityPlastic -0.091031 -0.304225  0.122307     1876 0.424  

#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(mixed1.Z$Sol)
autocorr(mixed1.Z$Sol)
autocorr(mixed1.Z$VCV)

xsim <- simulate(mixed1.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- mixed1.Z$VCV/apply(mixed1.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.013716100     0.025128789     0.959374242     0.001780869

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
#  2.512879

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.37161 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.3376208

# Proportion of variance explained by random factors
Sol <- mixed1.Z$Sol/apply(mixed1.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(mixed1.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2682486 -0.8267097  1.3604641

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(mixed1.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.26399349 -0.06253622  0.59979836 

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
#  Min        1Q    Median       3Q      Max  
#-6.1078  -1.3492  -0.2740   0.9869  13.5119

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   1.21735    0.37983   3.205  0.00166 **  
#Precision    -0.19242    0.04522  -4.255 3.71e-05 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 5.671266)

#Null deviance: 1003.9  on 148  degrees of freedom
#Residual deviance: 893.8  on 147  degrees of freedom
#AIC: 695.78

#Number of Fisher Scoring iterations: 2

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 161 (with 12 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.0857 [-0.1489; -0.0224] -2.66  0.0079

#Quantifying heterogeneity:
#tau^2 = 0.1291 [0.2261; 0.3910]; tau = 0.3593 [0.4755; 0.6253]
#I^2 = 89.4% [88.1%; 90.6%]; H = 3.08 [2.90; 3.26]

#Test of heterogeneity:
#      Q  d.f. p-value
#1513.07  160 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry

cols <- rep(0,162)
for (i in 1:162){
  if (i > 150){
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

shapes <- rep(0,162)
for (i in 1:162){
  if (i >= 151){
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
for(i in 1:9){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 6x7


#### Medians and 95% Credible Intervals ####
emmeans(mixed1.Z, ~ Sex, data=metadata, level = 0.95)
#Sex     emmean lower.HPD upper.HPD
#both     0.293   -0.0256     0.640
#females  0.227   -0.1135     0.571
#males    0.237   -0.0555     0.556
Zb <-  0.293;  lb <- -0.0256; ub <- 0.640
Zf <-  0.227;  lf <- -0.1135; uf <- 0.571
Zm <-  0.237;  lm <- -0.0555; um <- 0.556


emmeans(mixed1.Z, ~ Plasticity, data=metadata, level = 0.95)
#Plasticity emmean lower.HPD upper.HPD
#No          0.291   -0.0130     0.584
#Plastic     0.211   -0.0993     0.528
Znpl <- 0.291; lnpl <-  -0.0130; unpl <- 0.584
Zpl <- 0.211;  lpl <- -0.0993;  upl <- 0.528

p_significance(mixed1.Z$Sol, threshold = 0)
#Parameter         |   ps
#Sexboth           | 0.97
#Sexfemales        | 0.96
#Sexmales          | 0.97
#PlasticityPlastic | 0.79


#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.7,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
#Female
segments(lf,1.6,uf,1.6);
points(Zf,1.6,pch=16,col = "black",xpd=NA)
#Male
segments(lm,1.4,um,1.4);
points(Zm,1.4,pch=16,col = "black",xpd=NA)
#Both Sexes
segments(lb,1.2,ub,1.2);
points(Zb,1.2,pch=16,col = "black",xpd=NA)
#Plastic
segments(lpl,1,upl,1);
points(Zpl,1,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl,0.8,unpl,0.8);
points(Znpl,0.8,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Males", cex = 0.9, adj = c(0,0))
text(-2.1,1.2,"Both Sexes", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Non-Plastic", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.7,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
#Female
segments(lf-0.0857,1.6,uf-0.0857,1.6);
points(Zf-0.0857,1.6,pch=16,col = "black",xpd=NA)
#Male
segments(lm-0.0857,1.4,um-0.0857,1.4);
points(Zm-0.0857,1.4,pch=16,col = "black",xpd=NA)
#Both
segments(lb-0.0857,1.2,ub-0.0857,1.2);
points(Zb-0.0857,1.2,pch=16,col = "black",xpd=NA)
#Plastic
segments(lpl-0.0857,1,upl-0.0857,1);
points(Zpl-0.0857,1,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl-0.0857,0.8,unpl-0.0857,0.8);
points(Znpl-0.0857,0.8,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Males", cex = 0.9, adj = c(0,0))
text(-2.1,1.2,"Both Sexes", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Non-Plastic", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z for each Study #####
order <- read.csv("Excel Sheets/species.csv")
order <- cbind(order,rep(0,length(order$Study)),rep(0,length(order$Study)))
names(order)[41] <- "lci"
names(order)[42] <- "uci"

for (i in 1:length(order$Study)){
  order$uci[i] <- order$Fisher_Z[i] + 1.96*order$SE_Z[i]
  order$lci[i] <- order$Fisher_Z[i] - 1.96*order$SE_Z[i]
}

#Fisher Z Plot
colsa <- rep(0,length(order$Classification))
for (i in 1:length(order$Classification)){
  if (order$Classification[i] == "carotenoid"){
    colsa[i] <- "orange"
  } else if (order$Eu_Pheomelanin[i] == "eumelanin"){
    colsa[i] <- "black"
  } else if (order$Eu_Pheomelanin[i] == "pheomelanin"){
    colsa[i] <- "orangered3"
  } else if (order$Classification[i] == "unknown"){
    colsa[i] <- "darkorchid4"
  }  else if (order$Classification[i] == "pteridine"){
    colsa[i] <- "deeppink"
  } else if (order$Classification[i] == "structural"){
    colsa[i] <- "cornflowerblue"
  } 
}

shapesa <- rep(0,length(order$Plasticity))
for (i in 1:length(order$Vert_Invert)){
  if (order$Plasticity[i] == "Plastic"){
    shapesa[i] <- 20
  } else if (order$Plasticity[i] == "No"){
    shapesa[i] <- 18
  } 
}

plot(NA,xlim=c(-4,4),ylim=c(0,310),axes=F,ann=F)
axis(1)
abline(v=0)
for (i in 1:length(order$Study)){
  segments(order$lci[i],i*2,order$uci[i],i*2);
  points(order$Fisher_Z[i],i*2,pch=shapesa[i], col = colsa[i],xpd=NA)
  text(-4.1, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
}

#Export 15x10

#### Fisher Z by Phylogeny ####
plot(NA,xlim=c(-5,5),ylim=c(0,300),axes=F,ann=F)
axis(1)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(2,8,8,2), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(12,14,14,12), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(16,40,40,16), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(150,196,196,150), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(220,302,302,220), col = alpha("gray", 0.8), density = NA)

abline(v=0)
for (i in 1:length(order$Study)){
  segments(order$lci[i],i*2,order$uci[i],i*2);
  points(order$Fisher_Z[i],i*2,pch=shapesa[i],col = colsa[i],xpd=NA)
  text(-5.1, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
  text(3.5, i*2, order$Class[i], cex= 0.5, adj = c(0,0), font = 2)
}

#Export 15x10

