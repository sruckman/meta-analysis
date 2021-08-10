#MCMCglmm using Fisher Z not rho

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


################ Random Effects Model with strict prior ####
load("R Files/allRnd_Z.RDATA")

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
allRnd.Z <- MCMCglmm(Fisher_Z ~ 1, random = ~animal + Authors + Pattern + 
                          Age + Sex + Location + Season + Plasticity +
                          us(SE_Z):units + Eu_Pheomelanin,
                        data=metadata, pedigree = tree, 
                        nitt = 600000, thin = 100, burnin = 400000, 
                        prior = prior.ex, trunc = T)

#Save the model for later
#save(allRnd.Z, file = "allRnd_Z.RDATA")

#Summary of Results
summary(allRnd.Z)
#Iterations = 400001:599901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -231.5795 

#G-structure:  ~animal
#post.mean  l-95% CI u-95% CI eff.samp
#animal   0.01511 7.47e-10     0.06     1752

#~Authors
#post.mean  l-95% CI u-95% CI eff.samp
#Authors   0.03254 1.173e-07  0.09368     1574

#~Pattern
#post.mean  l-95% CI u-95% CI eff.samp
#Pattern  0.008495 1.361e-12  0.03236     1848

#~Age
#post.mean  l-95% CI u-95% CI eff.samp
#Age     0.1553 4.184e-11   0.4447     2000

#~Sex
#post.mean  l-95% CI u-95% CI eff.samp
#Sex   0.117 4.755e-08   0.2814     2000

#~Location
#post.mean  l-95% CI u-95% CI eff.samp
#Location   0.03207 3.153e-10   0.1243     1354

#~Season
#post.mean  l-95% CI u-95% CI eff.samp
#Season    0.05818 5.194e-10   0.2078     2000

#~Plasticity
#post.mean  l-95% CI u-95% CI eff.samp
#Plasticity     5.458 9.008e-11    1.122     2000

#~us(SE_Z):units
#post.mean l-95% CI u-95% CI eff.samp
#SE_Z:SE_Z.units   4.817    3.306     6.63     1519

#~Eu_Pheomelanin
#post.mean l-95% CI u-95% CI eff.samp
#Eu_Pheomelanin   0.07839 0.001814   0.2552     1781

#R-structure:  ~units
#post.mean l-95% CI u-95% CI eff.samp
#units      0.006606 0.001662  0.01393     2000

#Location effects: Fisher_Z ~ 1  
#post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)    0.3079  -0.9005   1.1110     2000 0.372


plot(allRnd.Z$Sol)

# Proportion of variance explained by random factors
rand <- allRnd.Z$VCV/apply(allRnd.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

### calcualte I^2 to quantify heterogenity ###
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#animal         Authors         Pattern         Age 
#0.002818506    0.006639465     0.001608806     0.019756405 
#Sex             Location        Season         Plasticity 
#0.013909832     0.005746989     0.009723418    0.036494379 
#SE_Z:SE_Z.units  Eu_Pheomelanin  units 
#0.888004816      0.014044925     0.001252459 

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
#  0.6639465  

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.2818506

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.2631553 

################ Mixed Effects Model and Tests for the Best Model ####
prior.ex3 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5),
                           G4 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))
###### Run the model without phylogeny ####
load("R Files/nophylogeny_Z.RDATA")

animal.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                   random = ~animal + Authors + Location + us(SE_Z):units, 
                   data=metadata,  
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex3, trunc = T)
animal.Z$DIC
#-233.1065

#Save the model for later
#save(animal.Z, file = "nophylogeny_Z.RDATA")

## Heterogenetiy
# Proportion of variance explained by random factors
rand <- animal.Z$VCV/apply(animal.Z$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
Rand_Var <- apply(rand,2,mean)
Rand_Var
#animal         Authors        Location SE_Z:SE_Z.units 
#0.002363271     0.011138863     0.011293295     0.973732329 
#units 
#0.001472242 

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
#  1.113886  

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.2363271

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
#animal 
#0.157821

###### Run the model without Authors ####
prior.ex4<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

authors.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                    random = ~animal + Location + us(SE_Z):units, 
                    data=metadata, pedigree = tree,
                    nitt = 600000, thin = 100, burnin = 400000, 
                    prior = prior.ex4, trunc = T)
authors.Z$DIC
#-225.2394

###### Run the model without Location ####
location.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                     random = ~animal + Authors + us(SE_Z):units, 
                     data=metadata, pedigree = tree,
                     nitt = 600000, thin = 100, burnin = 400000, 
                     prior = prior.ex4, trunc = T)
location.Z$DIC
#-234.5162

###### Run the model unweighted ####
weight.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                   random = ~animal + Authors + Location, 
                   data=metadata, pedigree = tree,
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex4, trunc = T)
weight.Z$DIC
#199.8237

###### Run the model using Season rather than Plasticity ####
#Run the model
season.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Season, 
                   random = ~animal + Authors + Location + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex3, trunc = T)
season.Z$DIC
#-230.4465

###### Run the model using Plasticity ####
#Run the model
plasticity.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                       random = ~animal + Authors + Location + us(SE_Z):units, 
                       data=metadata, pedigree = tree, 
                       nitt = 600000, thin = 100, burnin = 400000, 
                       prior = prior.ex3, trunc = T)
plasticity.Z$DIC
#-231.0236

###### Mixed Effects Model with Plasticity x Sex ####
#Run the model
pxs.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex 
                  + Plasticity + Sex*Plasticity, 
                random = ~animal + Authors + Location + us(SE_Z):units, 
                data=metadata, pedigree = tree, 
                nitt = 600000, thin = 100, burnin = 400000, 
                prior = prior.ex3, trunc = T)
pxs.Z$DIC
#-228.2787

###### Mixed Effects Model with Eu_Pheomelanin x Plasticity ####
#Run the model
ebyp.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex 
                   + Plasticity + Eu_Pheomelanin*Plasticity, 
                 random = ~animal + Authors + Location + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = prior.ex3, trunc = T)
ebyp.Z$DIC
#-227.6096

###### Mixed Effects Model with Eu_Pheomelanin x Sex ####
#Run the model
ebys.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity
                 + Sex*Eu_Pheomelanin, 
                 random = ~animal + Authors + Location + us(SE_Z):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = prior.ex3, trunc = T)
ebys.Z$DIC
#-229.7641

################## BEST MODEL ##############
load("R Files/mixed_Z.RDATA")

#prior with expanded parameters (Cuachy Distribution close to a Fisher Z)
prior.ex2<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model
mixed.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                   random = ~animal + Authors + us(SE_Z):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2, trunc = T)
#Save the model for later
#save(mixed.Z, file = "mixed_Z.RDATA")

summary(mixed.Z)
#Iterations = 400001:599901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -235.9867 

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

#### calcualte I^2 to quantify heterogenity ####
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#animal         Authors         SE_Z:SE_Z.units   units 
#0.006035244    0.012476836     0.980042588       0.001445332 

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
#  1.247684 

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 0.6035244

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.3024061

# Proportion of variance explained by random factors
Sol <- mixed.Z$Sol/apply(mixed.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(mixed.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2295228 -0.8336396  1.2975919 

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(mixed.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2366212 -0.1155468  0.5919508

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
#-15.7954   -1.1780   -0.0061    1.3155    6.6361 

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.87449    0.37307   2.344  0.02042 * 
#Precision   -0.11386    0.04343  -2.622  0.00968 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 5.819291)

#Null deviance: 895.43  on 148  degrees of freedom
#Residual deviance: 855.44  on 147  degrees of freedom
#AIC: 689.25

#Number of Fisher Scoring iterations: 2

Trim_Fill <- meta::trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 167 (with 18 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.0614 [-0.1164; -0.0063] -2.18  0.0289

#Quantifying heterogeneity:
#tau^2 = 0.0938 [0.1691; 0.2934]; tau = 0.3062 [0.4112; 0.5417]
#I^2 = 86.1% [84.3%; 87.8%]; H = 2.68 [2.52; 2.86]

#Test of heterogeneity:
#      Q  d.f. p-value
#1196.12  166 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry
cols <- rep(0,167)
for (i in 1:167){
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

shapes <- rep(0,167)
for (i in 1:167){
  if (i >= 150){
    shapes[i] <- 8
  } else if (metadata$Plasticity[i] == "Plastic"){
    shapes[i] <- 20
  } else if (metadata$Plasticity[i] == "No"){
    shapes[i] <- 18
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "Fisher Z", ylab = "Inverse SE", 
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
emmeans(mixed.Z, ~ Eu_Pheomelanin, data=metadata)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.159   -0.0998     0.434
#eumelanin       0.258   -0.0187     0.565
#pheomelanin     0.252   -0.0967     0.589
#unknown         0.672    0.3022     1.061
Zcar <- 0.159; lcar <- -0.0998; ucar <- 0.434
Zeu <- 0.258;  leu <-  -0.0187; ueu <- 0.565
Zph <- 0.252;  lph <-  -0.0967; uph <- 0.589
Zun <- 0.672;  lun <-   0.3022; uun <- 1.061

emmeans(mixed.Z, ~ Vert_Invert, data=metadata)
#Vert_Invert  emmean lower.HPD upper.HPD
#invertebrate  0.410    0.0480     0.865
#vertebrate    0.259   -0.0227     0.561
Zin <- 0.410;     lin <-  0.0480;   uin <- 0.865
Zvert <- 0.259; lvert <- -0.0227; uvert <- 0.561

emmeans(mixed.Z, ~ Sex, data=metadata)
#Sex     emmean lower.HPD upper.HPD
#both     0.414    0.1268     0.707
#females  0.216   -0.0979     0.539
#males    0.381    0.1084     0.654
Zboth <- 0.414; lboth <- 0.1268; uboth <- 0.707
Zf <- 0.216;    lf <- -0.0979;    uf <- 0.539
Zm <- 0.381;    lm <- 0.1084;    um <- 0.654

emmeans(mixed.Z, ~ Plasticity, data=metadata)
#Plasticity emmean lower.HPD upper.HPD
#No          0.337    0.0899     0.607
#Plastic     0.335    0.0229     0.630
Znpl <- 0.337; lnpl <- 0.0899; unpl <- 0.607
Zpl <- 0.335;  lpl <- 0.0229;  upl <- 0.630

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
text(-2.1,0.2,"Invertebrates", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,0,"Plastic", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.2,"Non-Plastic*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.4,"Both Sexes*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.6,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,-0.8,"Males*", cex = 0.9, adj = c(0,0), font = 2)

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.9,1.3),axes=F,ann=F)
axis(1)
#Fisher Z
#Carotenoid
segments(lcar-0.0614,1.2,ucar-0.0614,1.2);
points(Zcar-0.0614,1.2,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu-0.0614,1,ueu-0.0614,1);
points(Zeu-0.0614,1,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph-0.0614,0.8,uph-0.0614,0.8);
points(Zph-0.0614,0.8,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun-0.0614,0.6,uun-0.0614,0.6);
points(Zun-0.0614,0.6,pch=16,col = "black",xpd=NA)
#Vertebrates
segments(lvert-0.0614,0.4,uvert-0.0614,0.4);
points(Zvert-0.0614,0.4,pch=16,col = "black",xpd=NA)
#Invertebrates
segments(lin-0.0614,0.2,uin-0.0614,0.2);
points(Zin-0.0614,0.2,pch=16,col = "black",xpd=NA)
#Plastic
segments(lpl-0.0614,0,upl-0.0614,0);
points(Zpl-0.0614,0,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl-0.0614,-0.2,unpl-0.0614,-0.2);
points(Znpl-0.0614,-0.2,pch=16,col = "black",xpd=NA)
#Sex Both 
segments(lboth-0.0614,-0.4,uboth-0.0614,-0.4);
points(Zboth-0.0614,-0.4,pch=16,col = "black",xpd=NA)
#Females
segments(lf-0.0614,-0.6,uf-0.0614,-0.6);
points(Zf-0.0614,-0.6,pch=16,col = "black",xpd=NA)
#Males
segments(lm-0.0614,-0.8,um-0.0614,-0.8);
points(Zm-0.0614,-0.8,pch=16,col = "black",xpd=NA)

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
text(-2.1,0,"Plastic", cex = 0.9, adj = c(0,0))
text(-2.1,-0.2,"Non-Plastic*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.4,"Both Sexes*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.6,"Females", cex = 0.9, adj = c(0,0))
text(-2.1,-0.8,"Males*", cex = 0.9, adj = c(0,0), font = 2)

#Export 6x6
#### Fisher Z for each Study #####
order <- read.csv("Excel Sheets/species_order.csv")
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
