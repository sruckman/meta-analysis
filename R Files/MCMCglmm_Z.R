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
                        nitt = 5000000, thin = 1000, burnin = 2500000, 
                        prior = prior.ex)


#Save the model for later
#save(allRnd.Z, file = "allRnd_Z.RDATA")

#Summary of Results
summary(allRnd.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -238.5263 

#G-structure:  ~animal
#       post.mean  l-95% CI u-95% CI eff.samp
#animal   0.0251 9.512e-09  0.09833     2463

#~Authors
#        post.mean  l-95% CI u-95% CI eff.samp
#Authors    0.0824  0.01902   0.1487     2500

#~us(SE_Z):units
#              post.mean l-95% CI u-95% CI eff.samp
#SE_Z:SE_Z.units   3.671    2.474    4.995     2500

#R-structure:  ~units
#          post.mean l-95% CI u-95% CI eff.samp
#units      0.006336 0.001487  0.01327     2500

#Location effects: Fisher_Z ~ 1  
#             post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)    0.26003  0.04646  0.48348     2500 0.032 *


plot(allRnd.Z$Sol)

#proportion of samples above 0
p_significance(allRnd.Z$Sol,threshold = 0)
#Practical Significance (threshold: 0.00)
#Parameter   |   ps
#(Intercept) | 0.98


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
#0.006732780     0.023056499     0.968476479     0.001734243

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
# 2.30565   

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.673278

## total heterogeneity percent
(I2s*100)+(I2u*100)
#2.978928

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2135796 

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(allRnd.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.26002580 0.04646076 0.48347851  

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
#Precision   -0.12969    0.04588  -2.827  0.00536 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 6.257106)

#Null deviance: 969.79  on 148  degrees of freedom
#Residual deviance: 919.79  on 147  degrees of freedom
#AIC: 700.05

#Number of Fisher Scoring iterations: 2

Trim_Fill <- meta::trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 153 (with 4 added studies)

#                                         95%-CI     z  p-value
#Random effects model   -0.0255 [-0.0837; 0.0327] -0.86  0.3913

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

#### Adjusted mean and 95% confidence interval #####
#              fit               lwr              upr 
#0.26002580-0.0255 0.04646076-0.0255 0.48347851-0.0255
#        0.2345258        0.02096076         0.4579785
0.26002580-0.0255
0.04646076-0.0255
0.48347851-0.0255

#Mean and 95% CI as rho rather than Fisher Z
FisherZInv(0.26002580)
FisherZInv(0.04646076)
FisherZInv(0.48347851)
FisherZInv(0.26002580-0.0255)
FisherZInv(0.04646076-0.0255)
FisherZInv(0.48347851-0.0255)

cols <- rep("black",154)

shapes <- c(rep(16,149),rep(8,4))

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = shapes,
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


###### Run the model without species ####
priors <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))

species.Z <- MCMCglmm(Fisher_Z ~ 1, 
                    random = ~ Authors + us(SE_Z):units, 
                    data=metadata, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = priors)
species.Z$DIC
#-235.9753

###### Run the model without Authors ####
au.Z <- MCMCglmm(Fisher_Z ~ 1, 
                 random = ~ animal + us(SE_Z):units, 
                 data=metadata, pedigree = tree,
                 nitt = 2000000, thin = 1000, burnin = 1000000, 
                 prior = priors)
au.Z$DIC
#-223.0189

###### Run the model without weight ####
weight.Z <- MCMCglmm(Fisher_Z ~ 1, 
                 random = ~ animal + Authors, 
                 data=metadata, pedigree = tree,
                 nitt = 2000000, thin = 1000, burnin = 1000000, 
                 prior = priors)
weight.Z$DIC
#173.1812

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
#-235.5343

#Save the model for later
#save(class.Z, file = "class_Z.RDATA")

###### Run the model with Plasticity only ####
#Run the model
plasticity.Z <- MCMCglmm(Fisher_Z ~ Plasticity - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
plasticity.Z$DIC
#-238.9878

#Save the model for later
#save(plasticity.Z, file = "plasticity_Z.RDATA")

###### Run the model with Sex only ####
#Run the model
sex.Z <- MCMCglmm(Fisher_Z ~ Sex - 1, 
                         random = ~animal + Authors + us(SE_Z):units, 
                         data=metadata, pedigree = tree, 
                         nitt = 5000000, thin = 1000, burnin = 2500000, 
                         prior = prior.ex2)
sex.Z$DIC
#-237.0872

#Save the model for later
#save(sex.Z, file = "sex_Z.RDATA")

###### Run the model with Vert/Invert only ####
#Run the model
vert.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
vert.Z$DIC
#-238.1355

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
#-237.8867

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
#-230.2013

###### Run the model with Age only ####
age.Z <- MCMCglmm(Fisher_Z ~ Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
age.Z$DIC
#-236.7384

###### Run the model with Class and Vert/Invert ####
#Run the model
cv.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert - 1, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree, 
                     nitt = 5000000, thin = 1000, burnin = 2500000, 
                     prior = prior.ex2)
cv.Z$DIC
#-235.4292

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
#-235.1119

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
#-226.9562

###### Run the model with Class and Age ####
#Run the model
ca.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Age - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
ca.Z$DIC
#-234.4946

###### Run the model with Plasticity and Age ####
#Run the model
pa.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
pa.Z$DIC
#-236.5053

###### Run the model with Sex and Age ####
#Run the model
sa.Z <- MCMCglmm(Fisher_Z ~ Sex + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
sa.Z$DIC
#-235.4912

###### Run the model with Vert/Invert and Age ####
#Run the model
verta.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Age - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
verta.Z$DIC
#-236.2052

###### Run the model with Location and Age ####
#Run the model
la.Z <- MCMCglmm(Fisher_Z ~ Location + Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
la.Z$DIC
#-236.9072

###### Run the model with Season and Age ####
#Run the model
seaa.Z <- MCMCglmm(Fisher_Z ~ Season + Age - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
seaa.Z$DIC
#-230.7736

###### Run the model with Class and Plasticity ####
#Run the model
pc.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
pc.Z$DIC
#-234.8311

###### Run the model with Class and Sex ####
#Run the model
cs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
cs.Z$DIC
#-233.7354

###### Run the model with Plasticity and Sex ####
#Run the model
ps.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
ps.Z$DIC
#-237.1174

#Save the model for later
#save(ps.Z, file = "ps_Z.RDATA")

###### Run the model with Plasticity and Vert/Invert ####
#Run the model
pv.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
pv.Z$DIC
#-238.3022

#Save the model for later
#save(pv.Z, file = "pv_Z.RDATA")

###### Run the model with Plasticity and Location ####
#Run the model
ploc.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Location - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
ploc.Z$DIC
#-238.1164

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
#-229.9228

###### Run the model with Sex and Vert/Invert ####
#Run the model
svert.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert - 1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 5000000, thin = 1000, burnin = 2500000, 
                   prior = prior.ex2)
svert.Z$DIC
#-236.1474

###### Run the model with Sex and Location ####
#Run the model
seloc.Z <- MCMCglmm(Fisher_Z ~ Sex + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
seloc.Z$DIC
#-236.6283

###### Run the model with Sex and Season ####
#Run the model
ss.Z <- MCMCglmm(Fisher_Z ~ Sex + Season - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
ss.Z$DIC
#-230.3209

###### Run the model with Vert/Invert and Location ####
#Run the model
vl.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
vl.Z$DIC
#-236.6824

###### Run the model with Vert/Invert and Season ####
#Run the model
vsea.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Season - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
vsea.Z$DIC
#-230.2726

###### Run the model with Location and Season ####
#Run the model
sealoc.Z <- MCMCglmm(Fisher_Z ~ Season + Location - 1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
sealoc.Z$DIC
#-231.5948


###### Run the model with Class, Plasticity, and Sex ####
#Run the model
cps.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 5000000, thin = 1000, burnin = 2500000, 
                 prior = prior.ex2)
cps.Z$DIC
#-233.927

###### Run the model Class, Plasticity, and Vert/Invert ####
#Run the model
cpv.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity + Vert_Invert -1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 2000000, thin = 1000, burnin = 1000000, 
                   prior = prior.ex2)
cpv.Z$DIC
#-235.3293

###### Run the model Class, Sex, and Season ####
#Run the model
cssea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Season - 1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cssea.Z$DIC
#-226.4387

###### Run the model Class, Plasticity, and Location ####
#Run the model
cpl.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity + Location -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cpl.Z$DIC
#-235.0754

###### Run the model Class, Plasticity, and Season ####
#Run the model
cpsea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cpsea.Z$DIC
#-225.9835

###### Run the model Class, Sex, and Vert/Invert ####
#Run the model
csv.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Vert_Invert -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
csv.Z$DIC
#-232.5056

###### Run the model Class, Sex, and Location ####
#Run the model
csl.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Location -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
csl.Z$DIC
#-233.454

###### Run the model Class, Vert/Invert, and Location ####
#Run the model
cvl.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Location -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cvl.Z$DIC
#-234.8134

###### Run the model Class, Vert/Invert, and Season ####
#Run the model
cvsea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cvsea.Z$DIC
#-226.9591

###### Run the model Class, Location, and Season ####
#Run the model
clsea.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Location + Season -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
clsea.Z$DIC
#-229.3553

###### Run the model Class, Plasticity, and Age ####
#Run the model
cpa.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Plasticity + Age -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
cpa.Z$DIC
#-234.7057

###### Run the model Class, Sex, and Age ####
#Run the model
csa.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
csa.Z$DIC
#-233.6389

###### Run the model Class, Vert/Invert, and Age ####
#Run the model
cva.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cva.Z$DIC
#-232.0138

###### Run the model Class, Location, and Age ####
#Run the model
cla.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Location + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cla.Z$DIC
#-234.6989

###### Run the model Class, Season, and Age ####
#Run the model
cseaa.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Season + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
cseaa.Z$DIC
#-227.7061

###### Run the model Plasticity, Sex, and Age ####
#Run the model
psa.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex + Age -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
psa.Z$DIC
#-234.6685

###### Run the model Plasticity, Vert/Invert, and Age ####
#Run the model
pva.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pva.Z$DIC
#-235.7323

###### Run the model Plasticity, Location, and Age ####
#Run the model
pla.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Location + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
pla.Z$DIC
#-236.8264

###### Run the model Plasticity, Season, and Age ####
#Run the model
pseaa.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Season + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pseaa.Z$DIC
#-229.492

###### Run the model Plasticity, Sex, and Vert/Invert ####
#Run the model
psv.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex + Vert_Invert -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
psv.Z$DIC
#-236.4338

###### Run the model Plasticity, Sex, and Location ####
#Run the model
psl.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex + Location -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
psl.Z$DIC
#-236.4045

###### Run the model Plasticity, Sex, and Season ####
#Run the model
pssea.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Sex + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pssea.Z$DIC
#-229.5225

###### Run the model Plasticity, Vert/Invert, and Location ####
#Run the model
pvl.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert + Location -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
pvl.Z$DIC
#-237.9993

###### Run the model Plasticity, Vert/Invert, and Season ####
#Run the model
pvsea.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Vert_Invert + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
pvsea.Z$DIC
#-228.8145

###### Run the model Plasticity, Location, and Season ####
#Run the model
plsea.Z <- MCMCglmm(Fisher_Z ~ Plasticity + Location + Season -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
plsea.Z$DIC
#-230.6133

###### Run the model Sex, Vert/Invert, and Location ####
#Run the model
svl.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert + Location -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 5000000, thin = 1000, burnin = 2500000, 
                    prior = prior.ex2)
svl.Z$DIC
#-235.6853

###### Run the model Sex, Vert/Invert, and Season ####
#Run the model
svsea.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert + Season -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
svsea.Z$DIC
#-230.3167

###### Run the model Sex, Location, and Season ####
#Run the model
slsea.Z <- MCMCglmm(Fisher_Z ~ Sex + Location + Season -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
slsea.Z$DIC
#-230.5738

###### Run the model Sex, Vert/Invert, and Age ####
#Run the model
sva.Z <- MCMCglmm(Fisher_Z ~ Sex + Vert_Invert + Age -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
sva.Z$DIC
#-233.9941

###### Run the model Sex, Location, and Age ####
#Run the model
sla.Z <- MCMCglmm(Fisher_Z ~ Sex + Location + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
sla.Z$DIC
#-235.7001

###### Run the model Sex, Season, and Age ####
#Run the model
sseaa.Z <- MCMCglmm(Fisher_Z ~ Sex + Season + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
sseaa.Z$DIC
#-231.1525

###### Run the model Vert/Invert, Location, and Age ####
#Run the model
vla.Z <- MCMCglmm(Fisher_Z ~ Location + Vert_Invert + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 5000000, thin = 1000, burnin = 2500000, 
                  prior = prior.ex2)
vla.Z$DIC
#-236.0379

###### Run the model Season, Vert/Invert, and Age ####
#Run the model
seava.Z <- MCMCglmm(Fisher_Z ~ Season + Vert_Invert + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
seava.Z$DIC
#-230.5236

###### Run the model Location, Season, and Age ####
#Run the model
lseaa.Z <- MCMCglmm(Fisher_Z ~ Location + Season + Age -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
lseaa.Z$DIC
#-231.9506

###### Run the model Vert/Invert, Location, and Season ####
#Run the model
vlsea.Z <- MCMCglmm(Fisher_Z ~ Vert_Invert + Location + Season -1, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 2000000, thin = 1000, burnin = 1000000, 
                    prior = prior.ex2)
vlsea.Z$DIC
#-231.7949

###### Run the model Class, Plasticity, Sex, and Vert/Invert ####
#Run the model
vert.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity + Vert_Invert -1, 
                  random = ~animal + Authors + us(SE_Z):units, 
                  data=metadata, pedigree = tree, 
                  nitt = 2000000, thin = 1000, burnin = 1000000, 
                  prior = prior.ex2)
vert.Z$DIC
#-232.5562

###### Run the model Class, Plasticity, Sex and Location ####
#Run the model
loc.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity + Location -1, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 2000000, thin = 1000, burnin = 1000000, 
                   prior = prior.ex2)
loc.Z$DIC
#-233.2415

###### Mixed Effects Model Class, Plasticity, Sex, and Plasticity x Sex ####
#Run the model
pxs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity 
                  + Sex*Plasticity - 1, 
                random = ~animal + Authors + us(SE_Z):units, 
                data=metadata, pedigree = tree, 
                nitt = 2000000, thin = 1000, burnin = 1000000, 
                prior = prior.ex2)
pxs.Z$DIC
#-232.9638

###### Mixed Effects Model Class, Plasticity, Sex, and Class x Plasticity ####
#Run the model
ebyp.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity 
                   + Eu_Pheomelanin*Plasticity - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 2000000, thin = 1000, burnin = 1000000, 
                 prior = prior.ex2)
ebyp.Z$DIC
#-226.7052

###### Mixed Effects Model Class, Plasticity, Sex, and Class x Sex ####
#Run the model
ebys.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Sex + Plasticity
                 + Sex*Eu_Pheomelanin - 1, 
                 random = ~animal + Authors + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 2000000, thin = 1000, burnin = 1000000, 
                 prior = prior.ex2)
ebys.Z$DIC
#-229.9616


################## "BEST" MODEL (class) ##############
load("R Files/class_Z.RDATA")

color.Z <- class.Z

summary(color.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -235.5343 

#Location effects: Fisher_Z ~ Eu_Pheomelanin - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#Eu_Pheomelanincarotenoid    0.14759 -0.16082  0.39257     2500 0.2064  
#Eu_Pheomelanineumelanin     0.34029  0.07157  0.62119     2500 0.0168 *
#Eu_Pheomelaninpheomelanin   0.21560 -0.13966  0.57279     2500 0.1952  
#Eu_Pheomelaninstructural    0.33618 -0.05140  0.77266     2500 0.0952 .
#Eu_Pheomelaninunknown       0.43495  0.04792  0.82925     2500 0.0384 *
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
#0.008709969     0.020457038     0.969109368     0.001723625

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
#  2.045704

## species level heterogeneity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.8709969

## total heterogeneity
(I2s*100) + (I2u*100)
#2.916701

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2819615

# Proportion of variance explained by random factors
Sol <- color.Z$Sol/apply(color.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(color.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2612002 -0.7902991  1.3135728

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(color.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.27294299 -0.08404796  0.62026270 

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
#-5.9201  -1.4008  -0.2601   1.2131  13.1621    

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.44368    0.37079   1.197    0.233  
#Precision   -0.05648    0.04415  -1.279    0.203
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 5.794479)

#Null deviance: 861.27 on 148  degrees of freedom
#Residual deviance: 851.79  on 147  degrees of freedom
#AIC: 688.61

#Number of Fisher Scoring iterations: 2

#### Since the intercept of Egger's regression is not significant 
#### there is no publication bias. Trim and fill is only used for 
#### funnel plot creation

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 149 (with 0 added studies)

#                                       95%-CI    z p-value
#Random effects model 0.0160 [-0.0352; 0.0673] 0.61  0.5395

#Quantifying heterogeneity:
#tau^2 = 0.0692 [0.1131; 0.2046]; tau = 0.2630 [0.3363; 0.4523]
#I^2 = 82.8% [80.2%; 85.1%]; H = 2.41 [2.25; 2.59]

#Test of heterogeneity:
#     Q  d.f.  p-value
#860.08  148 < 0.0001

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
#Export 8x8


#### Medians and 95% Credible Intervals ####
emmeans(color.Z, ~ Eu_Pheomelanin, data=metadata, level = 0.95)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.143   -0.1608     0.393
#eumelanin       0.330    0.0716     0.621
#pheomelanin     0.203   -0.1397     0.573
#structural      0.334   -0.0514     0.773
#unknown         0.439    0.0479     0.829
#Point estimate displayed: median 
#HPD interval probability: 0.95
Zcar <- 0.143; lcar <- -0.1608; ucar <- 0.393
Zeu <- 0.330;  leu <-   0.0716; ueu <- 0.621
Zph <- 0.203;  lph <-  -0.1397; uph <- 0.573
Zst <- 0.334;  lst <-  -0.0514; ust <- 0.773
Zun <- 0.439;  lun <-   0.0479; uun <- 0.829

#Probability above 0
p_significance(color.Z$Sol, threshold = 0)
#Parameter                 |   ps
#Eu_Pheomelanincarotenoid  | 0.90
#Eu_Pheomelanineumelanin   | 0.99
#Eu_Pheomelaninpheomelanin | 0.90
#Eu_Pheomelaninstructural  | 0.95
#Eu_Pheomelaninunknown     | 0.98


#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.7,1.9),axes=F,ann=F)
axis(1)
#Fisher Z
#Overall Model
segments(-0.08404796,1.8,0.62026270,1.8)
points(0.27294299,1.8,pch=16,col="black",xpd=NA)
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
text(-2.1,0.8,"Unknown", cex = 0.9, adj = c(0,0))

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

#### calculate I^2 to quantify heterogeneity ####
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
################## "BEST" MODEL (plasticity and vert/invert) ##############
load("R Files/pv_Z.RDATA")

summary(pv.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -238.3022 

#Location effects: Fisher_Z ~ Plasticity - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#PlasticityNo            0.30746 -0.12834  0.75337     2500 0.168
#PlasticityPlastic       0.22155 -0.22211  0.73506     2500 0.350
#Vert_Invertvertebrate  -0.02325 -0.56563  0.55058     2500 0.929
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(pv.Z$Sol)
autocorr(pv.Z$Sol)
autocorr(pv.Z$VCV)

xsim <- simulate(pv.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- pv.Z$VCV/apply(pv.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.013461734     0.023363354     0.961446888     0.001728025

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
#  2.336335

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.346173 

## total heterogeneity
(I2s*100) + (I2u*100)
#3.682509

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.3491737

# Proportion of variance explained by random factors
Sol <- pv.Z$Sol/apply(pv.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(pv.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2647184 -0.8420542  1.3734664

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(pv.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#       fit        lwr        upr 
# 0.2564761 -0.1410408  0.6338541

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
#-6.2131  -1.4895  -0.5078   1.1383  13.4286    

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  1.00562    0.38521   2.611  0.00997 **
#Precision   -0.15023    0.04587  -3.275  0.00132 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 6.253749)

#Null deviance: 986.39  on 148  degrees of freedom
#Residual deviance: 919.30  on 147  degrees of freedom
#AIC: 699.97

#Number of Fisher Scoring iterations: 2

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 156 (with 7 added studies)

#                                         95%-CI     z  p-value
#Random effects model   -0.0488 [-0.1114; 0.0138] -1.53  0.1266

#Quantifying heterogeneity:
#tau^2 = 0.1214 [0.2165; 0.3778]; tau = 0.3485 [0.4653; 0.6147]
#I^2 = 89.0% [87.6%; 90.3%]; H = 3.02 [2.84; 3.21]

#Test of heterogeneity:
#     Q  d.f. p-value
#1414.78  155 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry

#### Adjusted mean and 95% confidence interval #####
#              fit               lwr              upr 
#0.2564761-0.0488   -0.1410408-0.0488 0.6338541-0.0488
#       0.2076761          -0.1898408        0.5850541
0.2564761-0.0488
-0.1410408-0.0488
0.6338541-0.0488

cols <- rep(0,156)
for (i in 1:156){
  if (i > 149){
    cols[i] <- "black"
  } else if (metadata$Plasticity[i] == "Plastic"){
    cols[i] <- "darkturquoise"
  } else if (metadata$Plasticity[i] == "No"){
    cols[i] <- "darkorchid4"
  } 
}

shapes <- rep(0,156)
for (i in 1:156){
  if (i > 149){
    shapes[i] <- 8
  } else if (metadata$Vert_Invert[i] == "vertebrate"){
    shapes[i] <- 15
  } else if (metadata$Vert_Invert[i] == "invertebrate"){
    shapes[i] <- 16
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Vertebrate", "Invertebrate",  "Plastic", "Non-Plastic", 
           "Trim and Fill Points")
Cols <- c("black", "black", "darkturquoise", "darkorchid4", "black")
points <- c(15,16,15,15,8)
ys <- c(20.1,19.3,18.5,17.7,16.9)
for(i in 1:5){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Medians and 95% Credible Intervals ####
emmeans(pv.Z, ~ Plasticity, data=metadata, level = 0.95)
#Plasticity emmean lower.HPD upper.HPD
#No          0.295  -0.00335     0.582
#Plastic     0.216  -0.13340     0.523
Znpla <- 0.282; lnpla <- 0.0398; unpla <- 0.537
Zpla <- 0.209;  lpla <- -0.0683;  upla <- 0.490

emmeans(pv.Z, ~ Vert_Invert, data = metadata, level = 0.95)
#Vert_Invert  emmean lower.HPD upper.HPD
#invertebrate  0.270    -0.171     0.710
#vertebrate    0.245    -0.139     0.618
Zin <- 0.270; lin <- -0.171; uin <- 0.710
Zve <- 0.245; lve <- -0.139; uve <- 0.618


#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.1,0.7),axes=F,ann=F)
axis(1)
#Fisher Z
#Plastic
segments(lpla,0.6,upla,0.6);
points(Zpla,0.6,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpla,0.4,unpla,0.4);
points(Znpla,0.4,pch=16,col = "black",xpd=NA)
#Vertebrate
segments(lve,0.2,uve,0.2);
points(Zve,0.2,pch=16,col = "black",xpd=NA)
#Invert
segments(lin,0,uin,0);
points(Zin,0,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.2,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,0,"Invertebrate", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.1,0.7),axes=F,ann=F)
axis(1)
#Fisher Z
#Plastic
segments(lpla-0.0488,0.6,upla-0.0488,0.6);
points(Zpla-0.0488,0.6,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpla-0.0488,0.4,unpla-0.0488,0.4);
points(Znpla-0.0488,0.4,pch=16,col = "black",xpd=NA)
#Vert
segments(lve-0.0488,0.2,uve-0.0488,0.2);
points(Zve-0.0488,0.2,pch=16,col = "black",xpd=NA)
#Invert
segments(lin-0.0488,0,uin-0.0488,0);
points(Zin-0.0488,0,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.6,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Invertebrate", cex = 0.9, adj = c(0,0))

#Export 6x6
################### THE FOLLOWING MODELS WERE NOT DESCRIBED IN PAPER ####
################## "BEST" MODEL (vert/invert) ##############
load("R Files/vert_Z.RDATA")

summary(vert.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -238.1355

#Location effects: Fisher_Z ~ Eu_Pheomelanin - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#Vert_Invertinvertebrate   0.29724 -0.13230  0.72998     2500 0.170  
#Vert_Invertvertebrate     0.25936 -0.04434  0.60560     2500 0.088 .
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(vert.Z$Sol)
autocorr(vert.Z$Sol)
autocorr(vert.Z$VCV)

xsim <- simulate(vert.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- vert.Z$VCV/apply(vert.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.009661749     0.023963123     0.964621051     0.001754078

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
#  2.396312

## species level heterogeneity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.9661749

## total heterogeneity
(I2s*100) + (I2u*100)
#3.362487 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2730931

# Proportion of variance explained by random factors
Sol <- vert.Z$Sol/apply(vert.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(vert.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2607961 -0.7997447  1.3241824

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(vert.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.26114108 -0.04846914  0.61144020 

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
#-6.5861  -1.4792  -0.3826   1.0668  13.5966    

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.85936    0.38540   2.230  0.02728 * 
#Precision   -0.12918    0.04589  -2.815  0.00555 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 6.259928)

#Null deviance: 969.82  on 148  degrees of freedom
#Residual deviance: 920.21  on 147  degrees of freedom
#AIC: 700.12

#Number of Fisher Scoring iterations: 2

#### Since the intercept of Egger's regression is not significant 
#### there is no publication bias. Trim and fill is only used for 
#### funnel plot creation

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 152 (with 3 added studies)

#                                       95%-CI    z p-value
#Random effects model -0.0165 [-0.0736; 0.0407] -0.56  0.5725

#Quantifying heterogeneity:
#tau^2 = 0.0943 [0.1697; 0.3000]; tau = 0.3070 [0.4119; 0.5477]
#I^2 = 86.6% [84.7%; 88.2%]; H = 2.73 [2.56; 2.91]

#Test of heterogeneity:
#     Q  d.f.  p-value
#1123.97  151 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry

cols <- rep(0,154)
for (i in 1:154){
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

shapes <- rep(0,154)
for (i in 1:154){
  if (i >= 151){
    shapes[i] <- 8
  } else if (metadatas$Vert_Invert[i] == "vertebrate"){
    shapes[i] <- 20
  } else if (metadatas$Vert_Invert[i] == "invertebrate"){
    shapes[i] <- 18
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown",
           "Structural", "Pteridine",  "Vertebrate", "Invertebrate", 
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


#### Adjusted mean and 95% confidence interval #####
#              fit               lwr              upr 
#0.26114108-0.0165 -0.04846914-0.0165 0.61144020-0.0165
#        0.2446411        -0.06496914         0.5949402
0.26114108-0.0165
-0.04846914-0.0165
0.61144020-0.0165

#### Medians and 95% Credible Intervals ####
emmeans(vert.Z, ~ Vert_Invert, data=metadata, level = 0.95)
#Vert_Invert emmean lower.HPD upper.HPD
#invertebrate  0.296   -0.1323     0.730
#vertebrate    0.256   -0.0443     0.606
#Point estimate displayed: median 
#HPD interval probability: 0.95
Zi <- 0.296; li <- -0.1323; ui <- 0.730
Zv <- 0.256; lv <- -0.0443; uv <- 0.606

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(1.3,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
#Vertebrate
segments(lv,1.6,uv,1.6);
points(Zv,1.6,pch=16,col = "black",xpd=NA)
#Invertebrate
segments(li,1.4,ui,1.4);
points(Zi,1.4,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Invertebrate", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.3,0.7),axes=F,ann=F)
axis(1)
#Fisher Z
#Vert
segments(lv-0.0165,0.6,uv-0.0165,0.6);
points(Zv-0.0165,0.6,pch=16,col = "black",xpd=NA)
#Invert
segments(li-0.0165,0.4,ui-0.0165,0.4);
points(Zi-0.0165,0.4,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Invertebrate", cex = 0.9, adj = c(0,0))

#Export 6x6
################## "BEST" MODEL (sex) ##############
load("R Files/sex_Z.RDATA")

summary(sex.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -237.0872 

#Location effects: Fisher_Z ~ Sex - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#Sexboth      0.30714  0.03720  0.59197     2256 0.0336 *
#Sexfemales   0.23147 -0.04780  0.51257     2333 0.0880 .
#Sexmales     0.25289  0.03239  0.52864     2786 0.0472 *  
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(sex.Z$Sol)
autocorr(sex.Z$Sol)
autocorr(sex.Z$VCV)

xsim <- simulate(sex.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- sex.Z$VCV/apply(sex.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.006366161     0.023249182     0.968664284     0.001720373

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
#  2.324918

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.6366161 

## total heterogeneity
(I2s*100) + (I2u*100)
#2.961534

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.2031599 

# Proportion of variance explained by random factors
Sol <- sex.Z$Sol/apply(sex.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(sex.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2636701 -0.7566854  1.2864617

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(sex.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.26284439 0.01908368 0.54187907 

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
#-6.5264  -1.4665  -0.2528   0.9647  13.7547     

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  1.06859    0.37884   2.821 0.005455 ** 
#Precision   -0.16523    0.04511  -3.663 0.000347 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 6.048785)

#Null deviance: 970.34  on 148  degrees of freedom
#Residual deviance: 889.17  on 147  degrees of freedom
#AIC: 695.01

#Number of Fisher Scoring iterations: 2

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 159 (with 10 added studies)

#                                         95%-CI     z  p-value
#Random effects model     -0.0695 [-0.1318; -0.0072] -2.19  0.0289

#Quantifying heterogeneity:
#tau^2 = 0.1228 [0.2148; 0.3733]; tau = 0.3504 [0.4635; 0.6109]
#I^2 = 89.0% [87.6%; 90.3%]; H = 3.02 [2.84; 3.21]

#Test of heterogeneity:
#     Q  d.f. p-value
#1439.20  158 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry

#### Adjusted mean and 95% confidence interval #####
#              fit               lwr              upr 
#0.26284439-0.0695 0.01908368-0.0695 0.54187907-0.0695
#        0.1933444       -0.05041632         0.4723791
0.26284439-0.0695
0.01908368-0.0695
0.54187907-0.0695

cols <- rep(0,159)
for (i in 1:159){
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

shapes <- rep(0,159)
for (i in 1:159){
  if (i >= 151){
    shapes[i] <- 8
  } else if (metadatas$Sex[i] == "both"){
    shapes[i] <- 16
  } else if (metadatas$Sex[i] == "females"){
    shapes[i] <- 17
  } else if (metadatas$Sex[i] == "males"){
    shapes[i] <- 18
    }
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
                         col = cols, pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown",
           "Structural", "Pteridine",  "Both Sexes", "Female", 
           "Male", "Trim and Fill Points")
Cols <- c("orange","black", "orangered3", "darkorchid4", 
          "cornflowerblue", "deeppink", "grey50", "grey50", "grey50",
          "forestgreen")
points <- c(15,15,15,15,15,15,16,17,18,8)
ys <- c(20.1,19.3,18.5,17.7,16.9,16.1,15.3,14.5,13.7,12.9)
for(i in 1:10){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Medians and 95% Credible Intervals ####
emmeans(sex.Z, ~ Sex, data=metadata, level = 0.95)
#Sex emmean lower.HPD upper.HPD
#both     0.306    0.0372     0.592
#females  0.227   -0.0478     0.513
#males    0.251    0.0324     0.529
#Point estimate displayed: median 
#HPD interval probability: 0.95
Zb <- 0.306;  lb <-  0.0372; ub <- 0.592
Zf <- 0.227;  lf <- -0.0478; uf <- 0.513
Zm <- 0.251;  lm <-  0.0324; um <- 0.529

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(1.1,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
#Both
segments(lb,1.6,ub,1.6);
points(Zb,1.6,pch=16,col = "black",xpd=NA)
#Females
segments(lf,1.4,uf,1.4);
points(Zf,1.4,pch=16,col = "black",xpd=NA)
#Males
segments(lm,1.2,um,1.2);
points(Zm,1.2,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Both Sexes", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,1.2,"Males", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(1.1,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
#Both
segments(lb-0.0695,1.6,ub-0.0695,1.6);
points(Zb-0.0695,1.6,pch=16,col = "black",xpd=NA)
#Females
segments(lf-0.0695,1.4,uf-0.0695,1.4);
points(Zf-0.0695,1.4,pch=16,col = "black",xpd=NA)
#Males
segments(lm-0.0695,1.2,um-0.0695,1.2);
points(Zm-0.0695,1.2,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Both Sexes", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,1.2,"Males", cex = 0.9, adj = c(0,0))

#Export 6x6
################## "BEST" MODEL (plasticity and sex) ##############
load("R Files/ps_Z.RDATA")

mixed1.Z <- ps.Z

summary(mixed1.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500

#DIC: -237.1174 

#Location effects: Fisher_Z ~ Sex + Plasticity - 1

#                  post.mean  l-95% CI  u-95% CI eff.samp pMCMC  
#PlasticityNo       0.350598  0.007813  0.645716     2500 0.032 *
#PlasticityPlastic  0.255839 -0.031520  0.616604     2500 0.110  
#Sexfemales        -0.077449 -0.372538  0.213645     2500 0.590  
#Sexmales          -0.070627 -0.301090  0.210541     2500 0.570   

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
#0.010340460     0.023339057     0.964584781     0.001735702

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
#  2.333906

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.034046 

## total heterogeneity
(I2s*100) + (I2u*100)
#3.367952 

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.2919779

# Proportion of variance explained by random factors
Sol <- mixed1.Z$Sol/apply(mixed1.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(mixed1.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2684798 -0.7873152  1.3294955

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(mixed1.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.26494888 -0.03240881  0.56967204 

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
#-6.0670  -1.4016  -0.3685   0.9692  13.5250

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   1.27409    0.37905   3.361 0.000989 ***
#Precision    -0.20269    0.04513  -4.491 1.42e-05 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 6.055397)

#Null deviance: 1012.28  on 148  degrees of freedom
#Residual deviance: 890.14  on 147  degrees of freedom
#AIC: 695.17

#Number of Fisher Scoring iterations: 2

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 162 (with 13 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.0927 [-0.1561; -0.0293] -2.86  0.0042

#Quantifying heterogeneity:
#tau^2 = 0.1310 [0.2281; 0.3937]; tau = 0.3620 [0.4776; 0.6275]
#I^2 = 89.5% [88.2%; 90.7%]; H = 3.09 [2.91; 3.28]

#Test of heterogeneity:
#      Q  d.f. p-value
#1537.61  161 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry

#### Adjusted mean and 95% confidence interval #####
#              fit               lwr              upr 
#0.26494888-0.0927 -0.03240881-0.0927 0.56967204-0.0927
#        0.1722489         -0.1251088         0.476972
0.26494888-0.0927
-0.03240881-0.0927
0.56967204-0.0927

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
#Export 8x8


#### Medians and 95% Credible Intervals ####
emmeans(mixed1.Z, ~ Sex, data=metadata, level = 0.95)
#Sex     emmean lower.HPD upper.HPD
#both     0.308    0.0115     0.616
#females  0.225   -0.0969     0.509
#males    0.232   -0.0582     0.492
Zbs <-  0.308;  lbs <-  0.0115; ubs <- 0.616
Zfs <-  0.225;  lfs <- -0.0969; ufs <- 0.509
Zms <-  0.232;  lms <- -0.0582; ums <- 0.492


emmeans(mixed1.Z, ~ Plasticity, data=metadata, level = 0.95)
#Plasticity emmean lower.HPD upper.HPD
#No          0.297    0.0548     0.592
#Plastic     0.210   -0.0834     0.489
Znp <- 0.297; lnp <-  0.0548; unp <- 0.592
Zp <-  0.210;  lp <- -0.0834;  up <- 0.489

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.7,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
#Female
segments(lfs,1.6,ufs,1.6);
points(Zfs,1.6,pch=16,col = "black",xpd=NA)
#Male
segments(lms,1.4,ums,1.4);
points(Zms,1.4,pch=16,col = "black",xpd=NA)
#Both Sexes
segments(lbs,1.2,ubs,1.2);
points(Zbs,1.2,pch=16,col = "black",xpd=NA)
#Plastic
segments(lp,1,up,1);
points(Zp,1,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnp,0.8,unp,0.8);
points(Znp,0.8,pch=16,col = "black",xpd=NA)

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
segments(lfs-0.0927,1.6,ufs-0.0927,1.6);
points(Zfs-0.0927,1.6,pch=16,col = "black",xpd=NA)
#Male
segments(lms-0.0927,1.4,ums-0.0927,1.4);
points(Zms-0.0927,1.4,pch=16,col = "black",xpd=NA)
#Both
segments(lbs-0.0927,1.2,ubs-0.0927,1.2);
points(Zbs-0.0927,1.2,pch=16,col = "black",xpd=NA)
#Plastic
segments(lp-0.0927,1,up-0.0927,1);
points(Zp-0.0927,1,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnp-0.0927,0.8,unp-0.0927,0.8);
points(Znp-0.0927,0.8,pch=16,col = "black",xpd=NA)

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
################## "BEST" MODEL (plasticity and vert/invert) ##############
load("R Files/pv_Z.RDATA")

summary(pv.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -238.3022 

#Location effects: Fisher_Z ~ Plasticity - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#PlasticityNo            0.30746 -0.12834  0.75337     2500 0.168
#PlasticityPlastic       0.22155 -0.22211  0.73506     2500 0.350
#Vert_Invertvertebrate  -0.02325 -0.56563  0.55058     2500 0.929
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(pv.Z$Sol)
autocorr(pv.Z$Sol)
autocorr(pv.Z$VCV)

xsim <- simulate(pv.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- pv.Z$VCV/apply(pv.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.013461734     0.023363354     0.961446888     0.001728025

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
#  2.336335

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.346173 

## total heterogeneity
(I2s*100) + (I2u*100)
#3.682509

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.3491737

# Proportion of variance explained by random factors
Sol <- pv.Z$Sol/apply(pv.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(pv.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2647184 -0.8420542  1.3734664

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(pv.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#       fit        lwr        upr 
# 0.2564761 -0.1410408  0.6338541

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
#-6.2131  -1.4895  -0.5078   1.1383  13.4286    

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  1.00562    0.38521   2.611  0.00997 **
#Precision   -0.15023    0.04587  -3.275  0.00132 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 6.253749)

#Null deviance: 986.39  on 148  degrees of freedom
#Residual deviance: 919.30  on 147  degrees of freedom
#AIC: 699.97

#Number of Fisher Scoring iterations: 2

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 156 (with 7 added studies)

#                                         95%-CI     z  p-value
#Random effects model   -0.0488 [-0.1114; 0.0138] -1.53  0.1266

#Quantifying heterogeneity:
#tau^2 = 0.1214 [0.2165; 0.3778]; tau = 0.3485 [0.4653; 0.6147]
#I^2 = 89.0% [87.6%; 90.3%]; H = 3.02 [2.84; 3.21]

#Test of heterogeneity:
#     Q  d.f. p-value
#1414.78  155 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry

#### Adjusted mean and 95% confidence interval #####
#              fit               lwr              upr 
#0.2564761-0.0488   -0.1410408-0.0488 0.6338541-0.0488
#       0.2076761          -0.1898408        0.5850541
0.2564761-0.0488
-0.1410408-0.0488
0.6338541-0.0488

cols <- rep(0,158)
for (i in 1:158){
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

shapes <- rep(0,157)
for (i in 1:157){
  if (i > 150){
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
emmeans(pv.Z, ~ Plasticity, data=metadata, level = 0.95)
#Plasticity emmean lower.HPD upper.HPD
#No          0.295  -0.00335     0.582
#Plastic     0.216  -0.13340     0.523
Znpla <- 0.282; lnpla <- 0.0398; unpla <- 0.537
Zpla <- 0.209;  lpla <- -0.0683;  upla <- 0.490

emmeans(pv.Z, ~ Vert_Invert, data = metadata, level = 0.95)
#Vert_Invert  emmean lower.HPD upper.HPD
#invertebrate  0.270    -0.171     0.710
#vertebrate    0.245    -0.139     0.618
Zin <- 0.270; lin <- -0.171; uin <- 0.710
Zve <- 0.245; lve <- -0.139; uve <- 0.618


#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.1,0.7),axes=F,ann=F)
axis(1)
#Fisher Z
#Plastic
segments(lpla,0.6,upla,0.6);
points(Zpla,0.6,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpla,0.4,unpla,0.4);
points(Znpla,0.4,pch=16,col = "black",xpd=NA)
#Vertebrate
segments(lve,0.2,uve,0.2);
points(Zve,0.2,pch=16,col = "black",xpd=NA)
#Invert
segments(lin,0,uin,0);
points(Zin,0,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.2,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,0,"Invertebrate", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.1,0.7),axes=F,ann=F)
axis(1)
#Fisher Z
#Plastic
segments(lpla-0.0488,0.6,upla-0.0488,0.6);
points(Zpla-0.0488,0.6,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpla-0.0488,0.4,unpla-0.0488,0.4);
points(Znpla-0.0488,0.4,pch=16,col = "black",xpd=NA)
#Vert
segments(lve-0.0488,0.2,uve-0.0488,0.2);
points(Zve-0.0488,0.2,pch=16,col = "black",xpd=NA)
#Invert
segments(lin-0.0488,0,uin-0.0488,0);
points(Zin-0.0488,0,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.6,"Vertebrate", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Invertebrate", cex = 0.9, adj = c(0,0))

#Export 6x6
################## "BEST" MODEL (plasticity and location) ##############
load("R Files/ploc_Z.RDATA")

summary(ploc.Z)
#Iterations = 2500001:4999001
#Thinning interval  = 1000
#Sample size  = 2500 

#DIC: -238.1164 

#Location effects: Fisher_Z ~ Plasticity - 1 
#                           post.mean l-95% CI u-95% CI eff.samp pMCMC  
#PlasticityNo        0.38021  0.13032  0.63781     2500 0.0088 **
#PlasticityPlastic   0.21526 -0.05218  0.45287     2500 0.0960 . 
#Locationdomestic   -0.34289 -0.74245  0.04987     2500 0.0840 . 
#Locationfield      -0.15402 -0.39670  0.08811     2500 0.2080   
#Locationlab         0.10155 -0.16918  0.42983     2500 0.4992
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(ploc.Z$Sol)
autocorr(ploc.Z$Sol)
autocorr(ploc.Z$VCV)

xsim <- simulate(ploc.Z)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- ploc.Z$VCV/apply(ploc.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calculate I^2 to quantify heterogeneity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal         Authors         SE_Z:SE_Z.units   units 
#0.005001111     0.028965095     0.964190292     0.001843502

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
#  2.89651 

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.5001111  

## total heterogeneity
(I2s*100) + (I2u*100)
#3.396621

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
#0.139658

# Proportion of variance explained by random factors
Sol <- ploc.Z$Sol/apply(ploc.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(ploc.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2374498 -0.7880431  1.2602636

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(ploc.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#       fit        lwr        upr 
# 0.23722721 -0.03642905  0.51857244

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
#-6.1704  -1.5852  -0.4085   1.3894  12.2355    

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  -0.08943    0.41406  -0.216    0.829
#Precision    0.05054    0.04930   1.025    0.307
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 7.225614)

#Null deviance: 1069.8  on 148  degrees of freedom
#Residual deviance: 1062.2  on 147  degrees of freedom
#AIC: 721.5

#Number of Fisher Scoring iterations: 2

#### Since the intercept of Egger's regression is not significant 
#### there is no publication bias. Trim and fill is only used for 
#### funnel plot creation

Trim_Fill <- trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 162 (with 13 added studies)

#                                         95%-CI     z  p-value
#Random effects model   0.1000 [0.0426; 0.1574] 3.41  0.0006

#Quantifying heterogeneity:
#tau^2 = 0.1028 [0.1490; 0.2604]; tau = 0.3207 [0.3860; 0.5103]
#I^2 = 87.2% [85.5%; 88.7%]; H = 2.79 [2.62; 2.97]

#Test of heterogeneity:
#     Q  d.f. p-value
#1254.28  161 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry


cols <- rep(0,150)
for (i in 1:150){
  if (i > 150){
    cols[i] <- "forestgreen"
  } else if (metadatas$Location[i] == "field"){
    cols[i] <- "orange"
  } else if (metadatas$Location[i] == "captivity"){
    cols[i] <- "black"
  } else if (metadatas$Location[i] == "lab"){
    cols[i] <- "dodgerblue"
  } else if (metadatas$Location[i] == "domestic"){
    cols[i] <- "darkorchid4"
  } 
}

shapes <- rep(0,150)
for (i in 1:150){
  if (i > 150){
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
words <- c("Field", "Captivity", "Lab", "Domestic",
           "Plastic", "Non-Plastic")
Cols <- c("orange","black", "dodgerblue", "darkorchid4", 
          "grey50", "grey50")
points <- c(15,15,15,15,20,18)
ys <- c(20.1,19.3,18.5,17.7,16.9,16.1)
for(i in 1:8){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 8x8


#### Medians and 95% Credible Intervals ####
emmeans(ploc.Z, ~ Plasticity, data=metadata, level = 0.95)
#Plasticity emmean lower.HPD upper.HPD
#No          0.282    0.0701     0.507
#Plastic     0.113   -0.1222     0.399
Znplas <- 0.282; lnplas <- 0.0701; unplas <- 0.507
Zplas <- 0.113;  lplas <- -0.1222;  uplas <- 0.399

emmeans(ploc.Z, ~ Location, data = metadata, level = 0.95)
#Vert_Invert  emmean lower.HPD upper.HPD
#Location  emmean lower.HPD upper.HPD
#captivity  0.298    0.0771     0.533
#domestic  -0.046   -0.4365     0.356
#field      0.138   -0.0889     0.395
#lab        0.401    0.0710     0.678
Zcap <- 0.298; lcap <- 0.0771; ucap <- 0.533
Zdom <- -0.045; ldom <- -0.4365; udom <- 0.356
Zfi <- 0.138; lfi <- -0.0889; ufi <- 0.395
Zlab <- 0.401; llab <- 0.0710; ulab <- 0.678

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.5,0.7),axes=F,ann=F)
axis(1)
#Fisher Z
#Plastic
segments(lplas,0.6,uplas,0.6);
points(Zplas,0.6,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnplas,0.4,unplas,0.4);
points(Znplas,0.4,pch=16,col = "black",xpd=NA)
#Captivity
segments(lcap,0.2,ucap,0.2);
points(Zcap,0.2,pch=16,col = "black",xpd=NA)
#Domestic
segments(ldom,0,udom,0);
points(Zdom,0,pch=16,col = "black",xpd=NA)
#Lab
segments(llab,-0.2,ulab,-0.2);
points(Zlab,-0.2,pch=16,col = "black",xpd=NA)
#Field
segments(lfi,-0.4,ufi,-0.4);
points(Zfi,-0.4,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,0.6,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.4,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,0.2,"Captivity", cex = 0.9, adj = c(0,0))
text(-2.1,0,"Domestic", cex = 0.9, adj = c(0,0))
text(-2.1,-0.2,"Lab", cex = 0.9, adj = c(0,0))
text(-2.1,-0.4,"Field", cex = 0.9, adj = c(0,0))

#Export 6x6
 
######## Rho Figure #######
plot(NA,xlim=c(-1.1,1.1),ylim=c(-0.1,0.8),axes=F,ann=F)
axis(1)

#### Random Effects Model
#Mean No correction
segments(FisherZInv(0.04646076),0.7,FisherZInv(0.48347851),0.7);
points(FisherZInv(0.26002580),0.7,pch=16,col = "black",xpd=NA)

#### Random Effects Model
#Mean with correction
segments(FisherZInv(0.04646076-0.0255),0.6,FisherZInv(0.48347851-0.0255),0.6);
points(FisherZInv(0.26002580-0.0255),0.6,pch=15,col = "grey75",xpd=NA)

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
points(FisherZInv(Zpl-0.0545),0.4,pch=15,col = "grey75",xpd=NA)
#Non-Plastic with correction
segments(FisherZInv(lnpl-0.0545),0.2,FisherZInv(unpl-0.0545),0.2);
points(FisherZInv(Znpl-0.0545),0.2,pch=15,col = "grey75",xpd=NA)
#Mean Plasticity Model with correction
segments(-0.051,0,0.436,0);
points(0.205,0,pch=15,col = "grey75",xpd=NA)

#Add line at 0 and separate models
abline(h = 0.55, lty = 1)
abline(v = 0, lty = 2)

#Add axis labels
title(xlab = "Correlation Coefficient")
text(-1,0.65,"Random Effects Model*", cex = 0.9, adj = c(0,0), font = 2)
text(-1,0.45,"Plastic", cex = 0.9, adj = c(0,0))
text(-1,0.25,"Non-Plastic*", cex = 0.9, adj = c(0,0))
text(-1,0.05,"Overall Model*", cex = 0.9, adj = c(0,0))

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

plot(NA,xlim=c(-4,4),ylim=c(0,310),axes=F,ann=F)
axis(1)
abline(v=0)
for (i in 1:length(order$Study)){
  segments(order$lci[i],i*2,order$uci[i],i*2);
  points(order$Fisher_Z[i],i*2,pch=shapesa[i], col = colsa[i],xpd=NA)
  text(-4.1, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
}

#Export 15x10

################## Fisher Z by Phylogeny ####
order <- read.csv("Excel Sheets/species.csv")
order <- cbind(order,rep(0,length(order$Study)),rep(0,length(order$Study)))
names(order)[41] <- "lci"
names(order)[42] <- "uci"

for (i in 1:length(order$Study)){
  order$uci[i] <- order$rho[i] + 1.96*order$SE_r[i]
  order$lci[i] <- order$rho[i] - 1.96*order$SE_r[i]
}

plot(NA,xlim=c(-2,2),ylim=c(1.9,300),axes=F,ann=F)
axis(1)
polygon(x = c(-2,-2,2,2), y = c(2,8,8,2), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-2,-2,2,2), y = c(12,14,14,12), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-2,-2,2,2), y = c(16,40,40,16), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-2,-2,2,2), y = c(150,196,196,150), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-2,-2,2,2), y = c(220,302,302,220), col = alpha("gray", 0.8), density = NA)

abline(v=0)
for (i in 1:length(order$Study)){
  segments((order$lci[i]),i*2,(order$uci[i]),i*2);
  points((order$rho[i]),i*2,pch=shapesa[i],col = colsa[i],xpd=NA)
  text(-2, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
  text(1.75, i*2, order$Class[i], cex= 0.5, adj = c(0,0), font = 2)
}

#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Structural", 
           "Pteridine", "Unknown", "Plastic", "Non-Plastic")
Cols <- c("darkorange","black", "tan4", "violet", 
          "turquoise4", "slateblue4", "black", "black")
points <- c(17,17,17,17,17,17,15,16)
ys <- c(140,135,130,125,120,115,110,105)
for(i in 1:8){
  points(x=1.75, y=ys[i], pch=points[i], col=Cols[i])
  text(x=1.75,y=ys[i], labels=words[i], pos=4, cex=.5, font = 2)
}
title(xlab = "Correlation Coefficient")

#Export 12x15

