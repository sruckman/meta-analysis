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

#DIC: -182.2326 

#G-structure:  ~animal
#       post.mean  l-95% CI u-95% CI eff.samp
#animal   0.03103 1.35e-07   0.1157     1544

#~Authors
#       post.mean  l-95% CI u-95% CI eff.samp
#Authors   0.0131 4.504e-10  0.03838     1469

#~Location
#          post.mean  l-95% CI u-95% CI eff.samp
#Location    0.03135 1.164e-10   0.1183     1060

#~us(SE):units
#              post.mean l-95% CI u-95% CI eff.samp
#SE:SE.units       5.058    3.548    6.596     2000

#R-structure:  ~units
#      post.mean l-95% CI u-95% CI eff.samp
#units   0.01112 0.001884  0.02751     1875

#Location effects: Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity 
#                          post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)                 0.27314 -0.15456  0.71858     2000 0.226
#Eu_Pheomelanineumelanin     0.11840 -0.05627  0.28893     2000 0.171
#Eu_Pheomelaninpheomelanin   0.05717 -0.18334  0.30560     2000 0.642
#Eu_Pheomelaninunknown       0.41694  0.13685  0.74891     2000 0.018 *
#Vert_Invertvertebrate      -0.10697 -0.56453  0.31942     2000 0.603
#Sexfemales                 -0.16732 -0.42735  0.06138     2000 0.196
#Sexmales                   -0.10210 -0.30223  0.10108     1535 0.332
#PlasticityPlastic           0.06515 -0.10563  0.21659     2135 0.438
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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
#0.006145046 0.002693605 0.005588416 0.983245529 0.002327404 

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
#  0.2693605

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.6145046

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.5503327

# Proportion of variance explained by random factors
Sol <- mixed1$Sol/apply(mixed1$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(mixed1, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.1975966 -0.8985029  1.2915561 

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
#-10.1320   -0.9268    0.0719    1.1480    9.2613  

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.83680    0.37016   2.261   0.0253 *
# Precision  -0.09595    0.04310  -2.226   0.0275 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 5.728923)

#Null deviance: 870.54  on 148  degrees of freedom
#Residual deviance: 842.14  on 147  degrees of freedom
#AIC: 686.91

#Number of Fisher Scoring iterations: 2

Trim_Fill <- meta::trimfill(MR, metadata$SE)
Trim_Fill

#Number of studies combined: k = 176 (with 27 added studies)

#                                         95%-CI     z p-value
#Random effects model -0.0652 [-0.1198; -0.0107] -2.34  0.0190

#Quantifying heterogeneity:
#  tau^2 = 0.0971 [0.1457; 0.2514]; tau = 0.3116 [0.3817; 0.5014]
#I^2 = 86.2% [84.4%; 87.8%]; H = 2.69 [2.53; 2.86]

#Test of heterogeneity:
#      Q d.f.  p-value
#1266.03  175 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry
cols <- rep(0,174)
for (i in 1:174){
  if (i >= 150){
    cols[i] <- "dodgerblue"
  } else if (metadata$Classification[i] == "carotenoid"){
    cols[i] <- "orange"
  } else if (metadata$Eu_Pheomelanin[i] == "eumelanin"){
    cols[i] <- "black"
  } else if (metadata$Eu_Pheomelanin[i] == "pheomelanin"){
    cols[i] <- "orangered3"
  } else if (metadata$Classification[i] == "unknown"){
    cols[i] <- "darkorchid4"
  } 
}

shapes <- rep(0,174)
for (i in 1:174){
  if (i >= 150){
    shapes[i] <- 8
  } else if (metadata$Plasticity[i] == "Plastic"){
    shapes[i] <- 20
  } else if (metadata$Plasticity[i] == "No"){
    shapes[i] <- 18
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Weight", 
                         col = cols, pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown", 
           "Plastic", "Non-Plastic", "Trim and Fill Points")
Cols <- c("orange","black", "orangered3", "darkorchid4", 
          "grey50", "grey50", "dodgerblue")
points <- c(15,15,15,15,20,18,8)
ys <- c(20.1,19.3,18.5,17.7,16.9,16.1,15.3)
for(i in 1:7){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 6x6

#### Medians and 95% Credible Intervals ####
emmeans(mixed1, ~ Eu_Pheomelanin, data=metadata)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.162    -0.1143     0.489
#eumelanin       0.283   -0.0533     0.570
#pheomelanin     0.212   -0.1404     0.581
#unknown         0.580    0.1668     0.976
Zcar <- 0.162; lcar <- -0.1143; ucar <- 0.489
Zeu <- 0.283; leu <- -0.0533; ueu <- 0.570
Zph <- 0.212; lph <- -0.1404; uph <- 0.581
Zun <- 0.580; lun <- 0.1668; uun <- 0.976

emmeans(mixed1, ~ Vert_Invert, data=metadata)
#Vert_Invert  emmean lower.HPD upper.HPD
#invertebrate  0.361   -0.0417     0.756
#vertebrate    0.255   -0.0653     0.606
Zin <- 0.361; lin <- -0.0417; uin <- 0.756
Zvert <- 0.255; lvert <- -0.0653; uvert <- 0.606

emmeans(mixed1, ~ Sex, data=metadata)
#Sex     emmean lower.HPD upper.HPD
#both     0.399    0.0959     0.702
#females  0.234   -0.1100     0.592
#males    0.298   -0.0258     0.618
Zboth <- 0.399; lboth <- 0.0959; uboth <- 0.702
Zf <- 0.234; lf <- -0.1100; uf <- 0.592
Zm <- 0.298; lm <- -0.0258; um <- 0.618

emmeans(mixed1, ~ Plasticity, data=metadata)
#Plasticity emmean lower.HPD upper.HPD
#No          0.279  -0.0476     0.575
#Plastic     0.339   0.0252     0.670
Znpl <- 0.279; lnpl <- -0.0476; unpl <- 0.575
Zpl <- 0.339; lpl <- 0.0252; upl <- 0.670

#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.9,1.3),axes=F,ann=F)
axis(1)
#Fisher Z
#Carotenoid
segments(lcar,1.2,ucar,1.2);
points(Zcar,1.2,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu,1,ueu,1);
points(Zeu,1,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph,0.8,uph,0.8);
points(Zph,0.8,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun,0.6,uun,0.6);
points(Zun,0.6,pch=16,col = "black",xpd=NA)
#Vertebrates
segments(lvert,0.4,uvert,0.4);
points(Zvert,0.4,pch=16,col = "black",xpd=NA)
#Invertebrates
segments(lin,0.2,uin,0.2);
points(Zin,0.2,pch=16,col = "black",xpd=NA)
#Plastic
segments(lpl,0,upl,0);
points(Zpl,0,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl,-0.2,unpl,-0.2);
points(Znpl,-0.2,pch=16,col = "black",xpd=NA)
#Sex Both 
segments(lboth,-0.4,uboth,-0.4);
points(Zboth,-0.4,pch=16,col = "black",xpd=NA)
#Females
segments(lf,-0.6,uf,-0.6);
points(Zf,-0.6,pch=16,col = "black",xpd=NA)
#Males
segments(lm,-0.8,um,-0.8);
points(Zm,-0.8,pch=16,col = "black",xpd=NA)

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
text(-2.1,-0.8,"Males", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.9,1.3),axes=F,ann=F)
axis(1)
#Fisher Z
#Carotenoid
segments(lcar-0.0652,1.2,ucar-0.0652,1.2);
points(Zcar-0.0652,1.2,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu-0.0652,1,ueu-0.0652,1);
points(Zeu-0.0652,1,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph-0.0652,0.8,uph-0.0652,0.8);
points(Zph-0.0652,0.8,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun-0.0652,0.6,uun-0.0652,0.6);
points(Zun-0.0652,0.6,pch=16,col = "black",xpd=NA)
#Vertebrates
segments(lvert-0.0652,0.4,uvert-0.0652,0.4);
points(Zvert-0.0652,0.4,pch=16,col = "black",xpd=NA)
#Invertebrates
segments(lin-0.0652,0.2,uin-0.0652,0.2);
points(Zin-0.0652,0.2,pch=16,col = "black",xpd=NA)
#Plastic
segments(lpl-0.0652,0,upl-0.0652,0);
points(Zpl-0.0652,0,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl-0.0652,-0.2,unpl-0.0652,-0.2);
points(Znpl-0.0652,-0.2,pch=16,col = "black",xpd=NA)
#Sex Both 
segments(lboth-0.0652,-0.4,uboth-0.0652,-0.4);
points(Zboth-0.0652,-0.4,pch=16,col = "black",xpd=NA)
#Females
segments(lf-0.0652,-0.6,uf-0.0652,-0.6);
points(Zf-0.0652,-0.6,pch=16,col = "black",xpd=NA)
#Males
segments(lm-0.0652,-0.8,um-0.0652,-0.8);
points(Zm-0.0652,-0.8,pch=16,col = "black",xpd=NA)

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
text(-2.1,-0.8,"Males", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z for each Study #####
order <- read.csv("Excel Sheets/species order.csv")
order <- cbind(order,rep(0,length(order$Study)),rep(0,length(order$Study)))
names(order)[41] <- "lci"
names(order)[42] <- "uci"

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

plot(NA,xlim=c(-4,4),ylim=c(0,310),axes=F,ann=F)
axis(1)
abline(v=0)
for (i in 1:length(order$Study)){
  segments(order$lci[i],i*2,order$uci[i],i*2);
  points(order$Fisher_Z[i],i*2,pch=shapesa[i],col = colsa[i],xpd=NA)
  text(-4.1, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
}

#Export 15x10

#### Fisher Z by Phylogeny ####
plot(NA,xlim=c(-5,5),ylim=c(0,300),axes=F,ann=F)
axis(1)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(2,8,8,2), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(10,16,16,10), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(18,46,46,18), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(150,196,196,150), col = alpha("gray", 0.8), density = NA)
polygon(x = c(-5.1,-5.1,4.5,4.5), y = c(220,300,300,220), col = alpha("gray", 0.8), density = NA)

abline(v=0)
for (i in 1:length(order$Study)){
  segments(order$lci[i],i*2,order$uci[i],i*2);
  points(order$Fisher_Z[i],i*2,pch=shapesa[i],col = colsa[i],xpd=NA)
  text(-5.1, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
  text(3.5, i*2, order$Class[i], cex= 0.5, adj = c(0,0), font = 2)
}

#Export 15x10