#MCMCglmm
setwd("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis")

library(ape)
library(MCMCglmm)
library(phytools)
library(TreeTools)
library(ggplot2)
library(emmeans)

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


#### Data ####
#read in the dataset with the calculated effect sizes
metadata <- read.csv("Excel Sheets/metafull.csv", header=T)

#Change Eu_Phae to represent all pigments
metadata$Eu_Pheomelanin <- ifelse(metadata$Eu_Pheomelanin=="N/A", metadata$Classification, metadata$Eu_Pheomelanin)
unique(metadata$Eu_Pheomelanin)

#we must change species to animal for the analysis to run
names(metadata)[5] <- "animal"
#Check to see if the labels on our tree match the animal column
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
##### Random Effects Model with strict prior ####
#priors:
prior.ex <- list(G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                         G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                         G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
                         G4 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                         G5 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                         G6 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                         G7 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5), 
                         G8 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
                         G9 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
                         G10 = list(V = 1, nu = 0.02)), 
                 R = list(V=1, nu=0.02))

#Run the model:
allRnd.mcmc <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + Pattern + 
                          Age + Sex + Location + Season + Plasticity 
                        + us(Weight):units + Eu_Pheomelanin,
                        data=metadata, pedigree = tree, 
                        nitt = 600000, thin = 100, burnin = 400000, 
                        prior = prior.ex, trunc = T)

#Save the model for later
save(allRnd.mcmc, file = "allRnd.RDATA")
load("allRnd.RDATA")

#Summary of Results
summary(allRnd.mcmc)
#Iterations = 400001:599901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: 192.6642 

#G-structure:  ~animal
#post.mean  l-95% CI u-95% CI eff.samp
#animal   0.02282 7.686e-10  0.09499     1745

#~Authors
#post.mean  l-95% CI u-95% CI eff.samp
#Authors   0.01556 3.635e-08  0.05394     2443

#~Pattern
#post.mean  l-95% CI u-95% CI eff.samp
#Pattern  0.009694 6.755e-10  0.03818     2000

#~Age
#post.mean  l-95% CI u-95% CI eff.samp
#Age    1.086 6.403e-12    0.373     2000

#~Sex
#post.mean  l-95% CI u-95% CI eff.samp
#Sex   0.06059 1.271e-09   0.1034     2000

#~Location
#post.mean  l-95% CI u-95% CI eff.samp
#Location   0.03652 3.907e-10   0.1173     2000

#~Season
#post.mean  l-95% CI u-95% CI eff.samp
#Season    1.175 1.801e-06     2.12     2000

#~Plasticity
#post.mean  l-95% CI u-95% CI eff.samp
#Plasticity     20.5 2.527e-09    1.389     2000

#~us(Weight):units
#post.mean l-95% CI u-95% CI eff.samp
#Weight:Weight.units 4.308e-05 2.102e-11 0.000166     2000

#~Eu_Pheomelanin
#post.mean l-95% CI u-95% CI eff.samp
#Eu_Pheomelanin   0.109 0.001323   0.2908     2000

#R-structure:  ~units
#post.mean l-95% CI u-95% CI eff.samp
#units     0.1815   0.1364   0.2313     2000

#Location effects: Fisher_Z ~ 1 
#post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)    0.4238  -0.7812   1.8894     2133 0.329


plot(allRnd.mcmc$Sol)

# Proportion of variance explained by random factors
rand <- allRnd.mcmc$VCV/apply(allRnd.mcmc$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calcualte I^2 to quantify heterogenity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's

# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#animal             Authors             Pattern 
#2.909328e-02        2.195610e-02        1.305600e-02 
#Age                 Sex            Location 
#6.888222e-02        2.685282e-02        3.589917e-02 
#Season          Plasticity Weight:Weight.units 
#3.864346e-01        9.602642e-02        5.735686e-05 
#Eu_Pheomelanin               units 
#8.577614e-02        2.359659e-01  

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[11]

#varM = samlng error effect
varM = Rand_Var[9]

#Other variance terms for total
var3 <- Rand_Var[3] # Pattern
var4 <- Rand_Var[4] #Age
var5 <- Rand_Var[5] #Sex
var6 <- Rand_Var[6] #Location
var7 <- Rand_Var[7] #Season
var8 <- Rand_Var[8] #Plasticity
var10 <- Rand_Var[10] #Eu_Pheomelanin

##Total variance equation
varT = varA + varS + varE + varM + var3 + var4 + var5 + var6 + var7 + var8 + var10


## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
#  2.19561  

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 2.909328

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.1013649

##### Testing Random Effects ####
prior.ex3 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5),
                          G4 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#### Run the model without phylogeny ####
animal <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Season, 
                   random = ~animal + Authors + Location + us(SE):units, 
                   data=metadata,  
                   nitt = 300000, thin = 100, burnin = 200000, 
                   prior = prior.ex3, trunc = T)
animal$DIC

## Heterogenetiy
# Proportion of variance explained by random factors
rand <- animal$VCV/apply(animal$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal     Authors    Location SE:SE.units       units 
#0.001836771 0.001697928 0.005380878 0.987884999 0.003199423

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[5]

#varM = samlng error effect
varM = Rand_Var[4]

#Other variance terms for total
var3 <- Rand_Var[3] # Location

##Total variance equation
varT = varA + varS + varE + varM + var3 

## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
#  0.1697928  

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 2.909328

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2


#### Run the model without Authors ####
prior.ex4<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

authors <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Season, 
                   random = ~animal + Location + us(SE):units, 
                   data=metadata, pedigree = tree,
                   nitt = 300000, thin = 100, burnin = 200000, 
                   prior = prior.ex4, trunc = T)
authors$DIC

#### Run the model without Location ####
location <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Season, 
                    random = ~animal + Authors + us(SE):units, 
                    data=metadata, pedigree = tree,
                    nitt = 300000, thin = 100, burnin = 200000, 
                    prior = prior.ex4, trunc = T)
location$DIC
#### Run the model unweighted ####
weight <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Season, 
                     random = ~animal + Authors + Location, 
                     data=metadata, pedigree = tree,
                     nitt = 300000, thin = 100, burnin = 200000, 
                     prior = prior.ex4, trunc = T)
weight$DIC
#### Run the model using Season rather than Plasticity ####
#prior with expanded parameters (Cuachy Distribution close to a Fisher Z)
prior.ex2<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5),
                          G4 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model
mixeds <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Season, 
                   random = ~animal + Authors + Location + us(SE):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2, trunc = T)
mixeds$DIC
##### Mixed Effects Model and strong prior using Plasticity ####
#prior with expanded parameters (Cuachy Distribution close to a Fisher Z)
prior.ex2<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5),
                          G4 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model
mixed1 <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                   random = ~animal + Authors + Location + us(SE):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2, trunc = T)
#Save the model for later
#save(mixed1, file = "mixed1.RDATA")
load("R Files/mixed1.RDATA")

summary(mixed1)
#Iterations = 400001:599901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -179.7628 

#G-structure:  ~animal
#post.mean  l-95% CI u-95% CI eff.samp
#animal   0.02989 6.629e-08   0.1062     1414

#~Authors
#post.mean  l-95% CI u-95% CI eff.samp
#Authors   0.01292 1.002e-09  0.03757     1707

#~Location
#post.mean  l-95% CI u-95% CI eff.samp
#Location    0.0287 5.758e-09   0.1024     2000

#~us(SE):units
#post.mean l-95% CI u-95% CI eff.samp
#SE:SE.units       5.1    3.606    6.696     2000

#R-structure:  ~units
#post.mean l-95% CI u-95% CI eff.samp
#units   0.01113 0.001534  0.02656     1778

#Location effects: Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity 
#                          post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)                 0.27409 -0.16427  0.72086     2164 0.211
#Eu_Pheomelanineumelanin     0.11803 -0.04916  0.29166     2000 0.179
#Eu_Pheomelaninpheomelanin   0.05907 -0.20167  0.29972     2000 0.643
#Eu_Pheomelaninunknown       0.42004  0.11640  0.70004     2000 0.008 **
#Vert_Invertvertebrate      -0.11121 -0.57538  0.34663     2000 0.579
#Sexfemales                 -0.16919 -0.41832  0.06963     2000 0.183
#Sexmales                   -0.10077 -0.32084  0.09929     2166 0.330
#PlasticityPlastic           0.06858 -0.10308  0.22533     2000 0.413
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Model diagnostics
plot(mixed1$Sol)
autocorr(mixed1$Sol)
autocorr(mixed1$VCV)

xsim <- simulate(mixed1)

model_test <- data.frame(metadata$Fisher_Z, metadata$Weight, xsim)
colnames(model_test) <- c("Fisher_Z", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=Fisher_Z, y=invSE)) +
  geom_point(data=model_test, aes(x=sim, y=invSE), color="red")

# Proportion of variance explained by random factors
rand <- mixed1$VCV/apply(mixed1$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### calcualte I^2 to quantify heterogenity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#     animal     Authors    Location SE:SE.units       units 
#0.005925896 0.002647712 0.005332847 0.983787733 0.002305812

#varA = phylogenetic level variance
varA <- Rand_Var[1] 

#varS = study level variance
varS = Rand_Var[2]

#varE = effect-size-specific effect
varE = Rand_Var[5]

#varM = samlng error effect
varM = Rand_Var[4]

#Other variance terms
var3 <- Rand_Var[3] # Location

##Total variance equation
varT = varA + varS + varE + varM + var3 


## study level heterogenity I^2s = varS/varT
I2s <- varS/varT
I2s*100
#  Authors 
#  0.2647712

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.5925896

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.5446886

# Proportion of variance explained by random factors
Sol <- mixed1$Sol/apply(mixed1$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(mixed1, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2045801 -0.8814606  1.2939119 

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(mixed1, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2040512 -0.1575069  0.5787719

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
#-10.0723   -0.9282    0.1388    1.2390    9.3065  

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.77841    0.37017   2.103   0.0372 *
#  Precision   -0.09287    0.04310  -2.155   0.0328 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 5.728923)

#Null deviance: 868.76  on 148  degrees of freedom
#Residual deviance: 842.15  on 147  degrees of freedom
#AIC: 686.91

#Number of Fisher Scoring iterations: 2

Trim_Fill <- meta::trimfill(MR, metadata$SE)
Trim_Fill

#Number of studies combined: k = 174 (with 25 added studies)

#95%-CI     z p-value
#Random effects model -0.0648 [-0.1193; -0.0102] -2.33  0.0200

#Quantifying heterogeneity:
#  tau^2 = 0.0957 [0.1436; 0.2486]; tau = 0.3093 [0.3789; 0.4986]
#I^2 = 86.1% [84.2%; 87.7%]; H = 2.68 [2.52; 2.85]

#Test of heterogeneity:
#  Q d.f.  p-value
#1242.67  173 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry
cols <- rep(0,length(metadata$Classification))
for (i in 1:length(metadata$Classification)){
  if (metadata$Classification[i] == "carotenoid"){
    cols[i] <- "orange"
  } else if (metadata$Eu_Pheomelanin[i] == "eumelanin"){
    cols[i] <- "black"
  } else if (metadata$Eu_Pheomelanin[i] == "pheomelanin"){
    cols[i] <- "orangered3"
  } else if (metadata$Classification[i] == "unknown"){
    cols[i] <- "darkorchid4"
  } 
}

shapes <- rep(0,length(metadata$Plasticity))
for (i in 1:length(metadata$Vert_Invert)){
  if (metadata$Plasticity[i] == "Plastic"){
    shapes[i] <- 20
  } else if (metadata$Plasticity[i] == "No"){
    shapes[i] <- 18
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Weight", 
                         col = alpha(cols, 0.8), pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown", 
           "Plastic", "Non-Plastic")
Cols <- c("orange","black", "orangered3", "darkorchid4", 
          "dodgerblue", "dodgerblue")
points <- c(15,15,15,15,20,18)
ys <- c(20.1,19.3,18.5,17.7,16.9,16.1)
for(i in 1:6){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 6x6

#### Medians and 95% Credible Intervals ####
emmeans(mixed1, ~ Eu_Pheomelanin, data=metadata)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.165    -0.119     0.471
#eumelanin       0.277    -0.018     0.599
#pheomelanin     0.217    -0.131     0.574
#unknown         0.582     0.218     0.988
Zcar <- 0.165; lcar <- -0.119; ucar <- 0.471
Zeu <- 0.277; leu <- -0.018; ueu <- 0.599
Zph <- 0.217; lph <- -0.131; uph <- 0.574
Zun <- 0.582; lun <- 0.218; uun <- 0.988

emmeans(mixed1, ~ Vert_Invert, data=metadata)
#Vert_Invert  emmean lower.HPD upper.HPD
#invertebrate  0.367   -0.0528     0.769
#vertebrate    0.253   -0.1049     0.559
Zin <- 0.367; lin <- -0.0528; uin <- 0.769
Zvert <- 0.253; lvert <- -0.1049; uvert <- 0.559

emmeans(mixed1, ~ Sex, data=metadata)
#Sex     emmean lower.HPD upper.HPD
#both     0.400   0.11828     0.730
#females  0.233  -0.10962     0.597
#males    0.293   0.00889     0.607
Zboth <- 0.400; lboth <- 0.11828; uboth <- 0.730
Zf <- 0.233; lf <- -0.10962; uf <- 0.597
Zm <- 0.293; lm <- 0.00889; um <- 0.607

emmeans(mixed1, ~ Plasticity, data=metadata)
#Plasticity emmean lower.HPD upper.HPD
#No          0.275  -0.00413     0.573
#Plastic     0.347   0.03384     0.670
Znpl <- 0.275; lnpl <- -0.00413; unpl <- 0.573
Zpl <- 0.347; lpl <- 0.03384; upl <- 0.670

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.9,1.3),axes=F,ann=F)
axis(1)
#Fisher Z
#Carotenoid
segments(lcar,1.2,ucar,1.2);
points(Zcar,1.2,pch=16,col = "orange",xpd=NA)
#Eumelanin
segments(leu,1,ueu,1);
points(Zeu,1,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph,0.8,uph,0.8);
points(Zph,0.8,pch=16,col = "orangered3",xpd=NA)
#Unknown
segments(lun,0.6,uun,0.6);
points(Zun,0.6,pch=16,col = "darkorchid4",xpd=NA)
#Vertebrates
segments(lvert,0.4,uvert,0.4);
points(Zvert,0.4,pch=16,col = "dodgerblue4",xpd=NA)
#Invertebrates
segments(lin,0.2,uin,0.2);
points(Zin,0.2,pch=16,col = "dodgerblue2",xpd=NA)
#Plastic
segments(lpl,0,upl,0);
points(Zpl,0,pch=16,col = "darkgreen",xpd=NA)
#Non-Plastic
segments(lnpl,-0.2,unpl,-0.2);
points(Znpl,-0.2,pch=16,col = "mediumseagreen",xpd=NA)
#Sex Both 
segments(lboth,-0.4,uboth,-0.4);
points(Zboth,-0.4,pch=16,col = "darkorchid2",xpd=NA)
#Females
segments(lf,-0.6,uf,-0.6);
points(Zf,-0.6,pch=16,col = "red",xpd=NA)
#Males
segments(lm,-0.8,um,-0.8);
points(Zm,-0.8,pch=16,col = "dodgerblue",xpd=NA)

#Add point for significant values


#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.2,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Eumelanin", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,0.6,"Unknown*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,0.4,"Vertebrates", cex = 0.9, adj = c(0,0))
text(-2.1,0.2,"Invertebrates", cex = 0.9, adj = c(0,0))
text(-2.1,0,"Plastic*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.2,"Non-Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,-0.4,"Both Sexes*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.6,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,-0.8,"Males*", cex = 0.9, adj = c(0,0), font = 2)

#Export 6x6
#### Fisher Z for each Study #####
order <- read.csv("Excel Sheets/species order.csv")
order <- cbind(order,rep(0,length(order$Study)),rep(0,length(order$Study)))
names(order)[40] <- "lci"
names(order)[41] <- "uci"

for (i in 1:length(order$Study)){
  order$uci[i] <- order$Fisher_Z[i] + 1.96*order$SE[i]
  order$lci[i] <- order$Fisher_Z[i] - 1.96*order$SE[i]
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

plot(NA,xlim=c(-4,4),ylim=c(0,160),axes=F,ann=F)
axis(1)
for (i in 1:length(order$Study)){
  segments(order$lci[i],i+3,order$uci[i],i+3);
  points(order$Fisher_Z[i],i+3,pch=shapesa[i],col = colsa[i],xpd=NA)
  text(-4.1, i+3, order$Study[i], cex = 0.5, adj = c(0,0))
}
abline(v=0)

#### Fisher Z by Phylogeny ####