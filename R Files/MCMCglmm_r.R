#MCMCglmm using r's only not Fisher Z

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
load("R Files/allRnd.RDATA")

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
allRnd.mcmc <- MCMCglmm(rho ~ 1, random = ~animal + Authors + Pattern + 
                          Age + Sex + Location + Season + Plasticity +
                          us(SE_r):units + Eu_Pheomelanin,
                        data=metadata, pedigree = tree, 
                        nitt = 600000, thin = 100, burnin = 400000, 
                        prior = prior.ex, trunc = T)

#Save the model for later
#save(allRnd.mcmc, file = "allRnd.RDATA")

#Summary of Results
summary(allRnd.mcmc)
#Iterations = 400001:599901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -9.884181 

#G-structure:  ~animal
#post.mean  l-95% CI u-95% CI eff.samp
#animal   0.0113 1.784e-09  0.04585     2193

#~Authors
#post.mean  l-95% CI u-95% CI eff.samp
#Authors   0.07029  0.03397   0.1094     1842

#~Pattern
#post.mean  l-95% CI u-95% CI eff.samp
#Pattern  0.006171 6.084e-09  0.02293     2000

#~Age
#post.mean  l-95% CI u-95% CI eff.samp
#Age    0.07098 3.192e-09   0.2498     2000

#~Sex
#post.mean  l-95% CI u-95% CI eff.samp
#Sex   2.41 2.912e-10   0.1032     2000

#~Location
#post.mean  l-95% CI u-95% CI eff.samp
#Location   0.02415 1.81e-08  0.08719     2000

#~Season
#post.mean  l-95% CI u-95% CI eff.samp
#Season    0.06076 1.61e-08   0.2303     2000

#~Plasticity
#post.mean  l-95% CI u-95% CI eff.samp
#Plasticity     2.853 1.283e-08     1.06     2000

#~us(SE_r):units
#post.mean l-95% CI u-95% CI eff.samp
#SE_r:SE_r.units   0.04064 4.206e-08   0.1794     1298

#~Eu_Pheomelanin
#post.mean l-95% CI u-95% CI eff.samp
#Eu_Pheomelanin   0.0773 0.002681   0.2749     2000

#R-structure:  ~units
#post.mean l-95% CI u-95% CI eff.samp
#units     0.03656  0.02533  0.05104     2000

#Location effects: rho ~ 1 
#post.mean l-95% CI u-95% CI eff.samp pMCMC
#(Intercept)    0.3457  -0.5180   1.0974     2000 0.246


plot(allRnd.mcmc$Sol)

# Proportion of variance explained by random factors
rand <- allRnd.mcmc$VCV/apply(allRnd.mcmc$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

### calcualte I^2 to quantify heterogenity ###
## As per Nakagawa & Santos 2012, calculate the different random effect-level I^2's
# Get the mean value of the variance for all the random effects
Rand_Var <- apply(rand,2,mean)
Rand_Var
#animal         Authors         Pattern             Age             Sex 
#0.02943910      0.20213982      0.01701554      0.08898916      0.05102322 
#Location          Season      Plasticity SE_r:SE_r.units  Eu_Pheomelanin 
#0.05188869      0.08697531      0.12488384      0.08836948      0.15286509 
#units 
#0.10641076 

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
#  20.21398  

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 2.94391

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.08710059

################ Mixed Effects Model and Tests for the Best Model ####
prior.ex3 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                           G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5),
                           G4 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                  R = list(V=1, nu=0.02))
###### Run the model without phylogeny ####
load("R Files/nophylogeny.RDATA")

animal <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                   random = ~animal + Authors + Location + us(SE_r):units, 
                   data=metadata,  
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex3, trunc = T)
animal$DIC
#-6.929693

#Save the model for later
#save(animal, file = "nophylogeny.RDATA")

## Heterogenetiy
# Proportion of variance explained by random factors
rand <- animal$VCV/apply(animal$VCV,1,sum)
# Get median values (50%) and 95% quantiles
apply(rand,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))
Rand_Var <- apply(rand,2,mean)
Rand_Var
#    animal         Authors        Location SE_r:SE_r.units           units 
#0.03042971      0.37188837      0.14055656      0.25407999      0.20304538

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
#  37.18884  

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 3.042971

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
#animal 
#0.05026684

###### Run the model without Authors ####
prior.ex4<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

authors <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                    random = ~animal + Location + us(SE_r):units, 
                    data=metadata, pedigree = tree,
                    nitt = 600000, thin = 100, burnin = 400000, 
                    prior = prior.ex4, trunc = T)
authors$DIC
#62.69042

###### Run the model without Location ####
location <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                     random = ~animal + Authors + us(SE_r):units, 
                     data=metadata, pedigree = tree,
                     nitt = 600000, thin = 100, burnin = 400000, 
                     prior = prior.ex4, trunc = T)
location$DIC
#-7.703434

###### Run the model unweighted ####
weight <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                   random = ~animal + Authors + Location, 
                   data=metadata, pedigree = tree,
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex4, trunc = T)
weight$DIC
#-6.847557

###### Run the model using Season rather than Plasticity ####
#prior with expanded parameters (Cuachy Distribution close to a Fisher Z)
prior.ex2<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5),
                          G4 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model
season <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Season, 
                   random = ~animal + Authors + Location + us(SE_r):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2, trunc = T)
season$DIC
#-6.668026

###### Mixed Effects Model and strong prior using Plasticity ####
#prior with expanded parameters (Cuachy Distribution close to a Fisher Z)
prior.ex2<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5),
                          G4 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model
plasticity <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                   random = ~animal + Authors + Location + us(SE_r):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2, trunc = T)
plasticity$DIC
#-7.263151

###### Mixed Effects Model with Plasticity x Sex ####
#Run the model
pxs <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity
                   + Sex*Plasticity, 
                   random = ~animal + Authors + Location + us(SE_r):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2, trunc = T)
pxs$DIC
#-4.438337

###### Mixed Effects Model with Eu_Pheomelanin x Plasticity ####
#Run the model
ebyp <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity
                + Eu_Pheomelanin*Plasticity, 
                random = ~animal + Authors + Location + us(SE_r):units, 
                data=metadata, pedigree = tree, 
                nitt = 600000, thin = 100, burnin = 400000, 
                prior = prior.ex2, trunc = T)
ebyp$DIC
#-6.880187

###### Mixed Effects Model with Eu_Pheomelanin x Sex ####
#Run the model
ebys <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity
                 + Sex*Eu_Pheomelanin, 
                 random = ~animal + Authors + Location + us(SE_r):units, 
                 data=metadata, pedigree = tree, 
                 nitt = 600000, thin = 100, burnin = 400000, 
                 prior = prior.ex2, trunc = T)
ebys$DIC
#-5.001037

################## BEST MODEL ##############
load("R Files/mixed1.RDATA")
#Run the model
mixed1 <- MCMCglmm(rho ~ Eu_Pheomelanin + Vert_Invert + Sex + Plasticity, 
                   random = ~animal + Authors + Location + us(SE_r):units, 
                   data=metadata, pedigree = tree, 
                   nitt = 600000, thin = 100, burnin = 400000, 
                   prior = prior.ex2, trunc = T)
#Save the model for later
#save(mixed1, file = "mixed1.RDATA")

summary(mixed1)
#Iterations = 400001:599901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -7.138998 

#Model diagnostics
plot(mixed1$Sol)
autocorr(mixed1$Sol)
autocorr(mixed1$VCV)

xsim <- simulate(mixed1)

model_test <- data.frame(metadata$rho, metadata$SE_r, xsim)
colnames(model_test) <- c("rho", "invSE", "sim")
ggplot() +
  geom_point(data=model_test, aes(x=rho, y=invSE)) +
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
#animal         Authors        Location   SE_r:SE_r.units       units
#0.06795781     0.36324093     0.13195816      0.24363206  0.19321104

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
#  36.32409 

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 6.795781

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.1088353

# Proportion of variance explained by random factors
Sol <- mixed1$Sol/apply(mixed1$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(mixed1, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2680767 -0.5877061  1.1313426 

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(mixed1, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.26658276 -0.07468926  0.60112149

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$SE_r
MR<-metadata$rho-pred_matrix[1:149]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Call:
#glm(formula = zMR ~ Precision, family = "gaussian", data = metadata)

#Deviance Residuals: 
#  Min        1Q    Median        3Q       Max  
#-0.169997  -0.029655   0.003043   0.030062   0.191675  

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  -0.023670   0.009795  -2.417   0.0169 *
# Precision    0.144388   0.058625   2.463   0.0149 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 0.002782622)

#Null deviance: 0.42592  on 148  degrees of freedom
#Residual deviance: 0.40905  on 147  degrees of freedom
#AIC: -449.94

#Number of Fisher Scoring iterations: 2

Trim_Fill <- meta::trimfill(MR, metadata$SE_r)
Trim_Fill

#Number of studies combined: k = 218 (with 69 added studies)

#                                         95%-CI     z  p-value
#Random effects model -0.3426 [-0.4295; -0.2557] -7.73 < 0.0001

#Quantifying heterogeneity:
#tau^2 = 0.4029 [0.2603; 0.4116]; tau = 0.6347 [0.5102; 0.6416]
#I^2 = 97.9% [97.7%; 98.0%]; H = 6.85 [6.64; 7.07]

#Test of heterogeneity:
#  Q d.f. p-value
#10182.55  217       0

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry
cols <- rep(0,218)
for (i in 1:218){
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

shapes <- rep(0,218)
for (i in 1:218){
  if (i >= 150){
    shapes[i] <- 8
  } else if (metadata$Plasticity[i] == "Plastic"){
    shapes[i] <- 20
  } else if (metadata$Plasticity[i] == "No"){
    shapes[i] <- 18
  } 
}

FTF_plot <- meta::funnel(Trim_Fill, yaxis="invse", xlim = c(-3,3),
                         xlab = "rho", ylab = "Inverse SE", 
                         col = cols, pch = shapes,
                         level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
#legend
words <- c("Carotenoid", "Eumelanin", "Pheomelanin", "Unknown", 
           "Plastic", "Non-Plastic", "Trim and Fill Points")
Cols <- c("orange","black", "orangered3", "darkorchid4", 
          "grey50", "grey50", "dodgerblue")
points <- c(15,15,15,15,20,18,8)
ys <- c(40.1,39.3,38.5,37.7,36.9,36.1,35.3)
for(i in 1:7){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 6x6


#### Medians and 95% Credible Intervals ####
emmeans(mixed1, ~ Eu_Pheomelanin, data=metadata)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.176   -0.1296     0.434
#eumelanin       0.330    0.0427     0.624
#pheomelanin     0.295   -0.0238     0.613
#unknown         0.646    0.3057     1.027
Zcar <- 0.176; lcar <- -0.1296; ucar <- 0.434
Zeu <- 0.330;  leu <-   0.0427; ueu <- 0.624
Zph <- 0.295;  lph <-  -0.0238; uph <- 0.613
Zun <- 0.646;  lun <-   0.3057; uun <- 1.027

emmeans(mixed1, ~ Vert_Invert, data=metadata)
#Vert_Invert  emmean lower.HPD upper.HPD
#invertebrate  0.429   0.08942     0.817
#vertebrate    0.295  -0.00223     0.590
Zin <- 0.429;     lin <-  0.08942;   uin <- 0.817
Zvert <- 0.295; lvert <- -0.00223; uvert <- 0.590

emmeans(mixed1, ~ Sex, data=metadata)
#Sex     emmean lower.HPD upper.HPD
#both     0.360    0.0806     0.664
#females  0.303    0.0112     0.613
#males    0.424    0.1174     0.705
Zboth <- 0.360; lboth <- 0.0806; uboth <- 0.664
Zf <- 0.303;    lf <- 0.0112;    uf <- 0.613
Zm <- 0.298;    lm <- 0.1174;    um <- 0.705

emmeans(mixed1, ~ Plasticity, data=metadata)
#Plasticity emmean lower.HPD upper.HPD
#No          0.359    0.0775     0.632
#Plastic     0.364    0.0644     0.649
Znpl <- 0.359; lnpl <- 0.0775; unpl <- 0.632
Zpl <- 0.364;  lpl <- 0.0644;  upl <- 0.649

#Rho Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.9,1.3),axes=F,ann=F)
axis(1)
#Rho
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
title(xlab = "Correlation Coefficient (r)")
text(-2.1,1.2,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Eumelanin", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,0.8,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,0.6,"Unknown", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,0.4,"Vertebrates", cex = 0.9, adj = c(0,0))
text(-2.1,0.2,"Invertebrates", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,0,"Plastic", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.2,"Non-Plastic", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.4,"Both Sexes", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.6,"Females", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,-0.8,"Males", cex = 0.9, adj = c(0,0), font = 2)

#Export 6x8
#### Rho with Publication Bias correction ####
#Rho Plot
plot(NA,xlim=c(-2,2),ylim=c(-0.9,1.3),axes=F,ann=F)
axis(1)
#Rho
#Carotenoid
segments(lcar-0.3426,1.2,ucar-0.3426,1.2);
points(Zcar-0.3426,1.2,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu-0.3426,1,ueu-0.3426,1);
points(Zeu-0.3426,1,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph-0.3426,0.8,uph-0.3426,0.8);
points(Zph-0.3426,0.8,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun-0.3426,0.6,uun-0.3426,0.6);
points(Zun-0.3426,0.6,pch=16,col = "black",xpd=NA)
#Vertebrates
segments(lvert-0.3426,0.4,uvert-0.3426,0.4);
points(Zvert-0.3426,0.4,pch=16,col = "black",xpd=NA)
#Invertebrates
segments(lin-0.3426,0.2,uin-0.3426,0.2);
points(Zin-0.3426,0.2,pch=16,col = "black",xpd=NA)
#Plastic
segments(lpl-0.3426,0,upl-0.3426,0);
points(Zpl-0.3426,0,pch=16,col = "black",xpd=NA)
#Non-Plastic
segments(lnpl-0.3426,-0.2,unpl-0.3426,-0.2);
points(Znpl-0.3426,-0.2,pch=16,col = "black",xpd=NA)
#Sex Both 
segments(lboth-0.3426,-0.4,uboth-0.3426,-0.4);
points(Zboth-0.3426,-0.4,pch=16,col = "black",xpd=NA)
#Females
segments(lf-0.3426,-0.6,uf-0.3426,-0.6);
points(Zf-0.3426,-0.6,pch=16,col = "black",xpd=NA)
#Males
segments(lm-0.3426,-0.8,um-0.3426,-0.8);
points(Zm-0.3426,-0.8,pch=16,col = "black",xpd=NA)

#Add point for significant values


#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Correlation Coefficient (r)")
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
#### Rho for each Study #####
order <- read.csv("Excel Sheets/species_order.csv")
order <- cbind(order,rep(0,length(order$Study)),rep(0,length(order$Study)))
names(order)[41] <- "lci"
names(order)[42] <- "uci"

for (i in 1:length(order$Study)){
  order$uci[i] <- order$rho[i] + 1.96*order$SE_r[i]
  order$lci[i] <- order$rho[i] - 1.96*order$SE_r[i]
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
  points(order$rho[i],i*2,pch=shapesa[i],col = colsa[i],xpd=NA)
  text(-4.1, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
}

#Export 15x10

#### Rho by Phylogeny ####
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
  points(order$rho[i],i*2,pch=shapesa[i],col = colsa[i],xpd=NA)
  text(-5.1, i*2, order$Study[i], cex = 0.5, adj = c(0,0))
  text(3.5, i*2, order$Class[i], cex= 0.5, adj = c(0,0), font = 2)
}

#Export 15x10
