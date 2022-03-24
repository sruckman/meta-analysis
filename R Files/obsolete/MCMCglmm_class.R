################## CLASS ONLY MODEL ##############
load("R Files/class_only_Z.RDATA")

#prior with expanded parameters (Cuachy Distribution close to a Fisher Z)
prior.ex2<- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G2 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0.4, alpha.V = 0.5)), 
                 R = list(V=1, nu=0.02))

#Run the model
mixed.Z <- MCMCglmm(Fisher_Z ~ Eu_Pheomelanin, 
                    random = ~animal + Authors + us(SE_Z):units, 
                    data=metadata, pedigree = tree, 
                    nitt = 800000, thin = 100, burnin = 600000, 
                    prior = prior.ex2)
#Save the model for later
#save(mixed.Z, file = "class_only_Z.RDATA")

summary(mixed.Z)
#Iterations = 400001:599901
#Thinning interval  = 100
#Sample size  = 2000 

#DIC: -237.188 

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
#0.014855439     0.022474335     0.960898199     0.001772027 

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
#  2.247434

## species level heterogenity I^2s = varA/varT
I2u <- varA/varT
I2u*100
# animal 
# 1.485544  

## phylogenetic heritability, phylogenetic signal H2 = varA/varA + varS + varE
H2 = varA/(varA + varS + varE)
H2
# animal 
# 0.379917

# Proportion of variance explained by random factors
Sol <- mixed.Z$Sol/apply(mixed.Z$Sol,1,sum)
# Get median values (50%) and 95% quantiles
apply(Sol,2,function(c) quantile(c,probs = c(0.025,0.5,0.975)))

#### Get the prediction interval of the overall effect size ####
pred_matrix <- predict(mixed.Z, interval = "prediction")
pred_interval <- apply(pred_matrix, 2, mean)
pred_interval
#      fit        lwr        upr 
#0.2627427 -0.8443200  1.3678141 

#### Get the posterior mean and 95% CI of the overall effect size ####
pred_matrix <- predict(mixed.Z, interval = "confidence")
confidence <- apply(pred_matrix, 2, mean)
confidence
#      fit        lwr        upr 
#0.2588482 -0.1019163  0.6101872 

#### publication bias ####
# getting predictions (raw data - predictions = meta-analytic residuals)
Precision<-metadata$Weight
MR<-metadata$Fisher_Z-pred_matrix[1:150]
zMR<-MR*Precision
metadata[,c("zMR","Precision")]<-c(zMR,Precision)

# Egger's regression
Egger<-glm(zMR~Precision,family="gaussian",data=metadata)
summary(Egger)

#Call:
#glm(formula = zMR ~ Precision, family = "gaussian", data = metadata)

#Deviance Residuals: 
#  Min        1Q    Median        3Q       Max  
#-5.8825  -1.4309  -0.2376   1.2093  13.1286 

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.55791    0.36768   1.517   0.1313 
#Precision   -0.07533    0.04385  -1.718   0.0879 .
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 5.671266)

#Null deviance: 863.12  on 149  degrees of freedom
#Residual deviance: 846.24  on 148  degrees of freedom
#AIC: 691.21

#Number of Fisher Scoring iterations: 2

Trim_Fill <- meta::trimfill(MR, metadata$SE_Z)
Trim_Fill

#Number of studies combined: k = 150 (with 0 added studies)

#                                         95%-CI     z  p-value
#Random effects model 0.0152 [-0.0358; 0.0661] 0.58  0.5594

#Quantifying heterogeneity:
#tau^2 = 0.0688 [0.1131; 0.2041]; tau = 0.2622 [0.3363; 0.4518]
#I^2 = 82.7% [80.0%; 85.0%]; H = 2.40 [2.24; 2.58]

#Test of heterogeneity:
#      Q  d.f. p-value
#859.41  149 < 0.0001

#Details on meta-analytical method:
#- Inverse variance method
#- DerSimonian-Laird estimator for tau^2
#- Jackson method for confidence interval of tau^2 and tau
#- Trim-and-fill method to adjust for funnel plot asymmetry
cols <- rep(0,150)
for (i in 1:150){
  if (i >= 151){
    cols[i] <- "forestgreen"
  } else if (metadata$Classification[i] == "carotenoid"){
    cols[i] <- "orange"
  } else if (metadata$Eu_Pheomelanin[i] == "eumelanin"){
    cols[i] <- "black"
  } else if (metadata$Eu_Pheomelanin[i] == "pheomelanin"){
    cols[i] <- "orangered3"
  } else if (metadata$Classification[i] == "unknown"){
    cols[i] <- "darkorchid4"
  } else if (metadata$Classification[i] == "structural"){
    cols[i] <- "cornflowerblue"
  } else if (metadata$Classification[i] == "pteridine"){
    cols[i] <- "deeppink"
  } 
}

shapes <- rep(0,150)
for (i in 1:150){
  if (i >= 151){
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
           "Structural", "Pteridine", "Plastic", "Non-Plastic")
Cols <- c("orange","black", "orangered3", "darkorchid4", 
          "cornflowerblue", "deeppink", "grey50", "grey50")
points <- c(15,15,15,15,15,15,20,18)
ys <- c(20.1,19.3,18.5,17.7,16.9,16.1,15.3,14.5)
for(i in 1:9){
  points(x=-2.9, y=ys[i], pch=points[i], col=Cols[i])
  text(x=-2.9,y=ys[i], labels=words[i], pos=4,cex=.75)
}
#Export 6x6


#### Medians and 95% Credible Intervals ####
emmeans(mixed.Z, ~ Eu_Pheomelanin, data=metadata)
#Eu_Pheomelanin emmean lower.HPD upper.HPD
#carotenoid      0.1452   -0.2091     0.417
#eumelanin       0.3354    0.0310     0.656
#pheomelanin     0.2103   -0.1567     0.673
#pteridine      -0.0923   -1.0118     0.796
#structural      0.3385   -0.0792     0.802
#unknown         0.4199   -0.0087     0.831
Zcar <- 0.1452; lcar <- -0.2091; ucar <- 0.417
Zeu <- 0.3354;  leu <-   0.0310; ueu <- 0.656
Zph <- 0.2103;  lph <-  -0.1567; uph <- 0.673
Zpt <- -0.0923; lpt <-  -1.0118; upt <- 0.796
Zst <- 0.3385;  lst <-  -0.0792; ust <- 0.802
Zun <- 0.4199;  lun <-  -0.0087; uun <- 0.831


#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.4,1.8),axes=F,ann=F)
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
#Pteridine
segments(lpt,1,upt,1);
points(Zpt,1,pch=16,col = "black",xpd=NA)
#Structural
segments(lst,0.8,ust,0.8);
points(Zst,0.8,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun,0.6,uun,0.6);
points(Zun,0.6,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Eumelanin*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,1.2,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Pteridine", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Structural", cex = 0.9, adj = c(0,0))
text(-2.1,0.6,"Unknown", cex = 0.9, adj = c(0,0))

#Export 6x6
#### Fisher Z with Publication Bias correction ####
#Fisher Z Plot
plot(NA,xlim=c(-2,2),ylim=c(0.4,1.8),axes=F,ann=F)
axis(1)
#Fisher Z
#Carotenoid
segments(lcar-0.0152,1.6,ucar-0.0152,1.6);
points(Zcar-0.0152,1.6,pch=16,col = "black",xpd=NA)
#Eumelanin
segments(leu-0.0152,1.4,ueu-0.0152,1.4);
points(Zeu-0.0152,1.4,pch=16,col = "black",xpd=NA)
#Pheomelanin
segments(lph-0.0152,1.2,uph-0.0152,1.2);
points(Zph-0.0152,1.2,pch=16,col = "black",xpd=NA)
#Pteridine
segments(lpt-0.0152,1,upt-0.0152,1);
points(Zpt-0.0152,1,pch=16,col = "black",xpd=NA)
#Structural 
segments(lst-0.0152,0.8,ust-0.0152,0.8);
points(Zst-0.0152,0.8,pch=16,col = "black",xpd=NA)
#Unknown
segments(lun-0.0152,0.6,uun-0.0152,0.6);
points(Zun-0.0152,0.6,pch=16,col = "black",xpd=NA)

#Add dashed line at 0
abline(v = 0, lty = 1)
#Add axis labels
title(xlab = "Fisher Z")
text(-2.1,1.6,"Carotenoid", cex = 0.9, adj = c(0,0))
text(-2.1,1.4,"Eumelanin*", cex = 0.9, adj = c(0,0), font = 2)
text(-2.1,1.2,"Pheomelanin", cex = 0.9, adj = c(0,0))
text(-2.1,1,"Pteridine", cex = 0.9, adj = c(0,0))
text(-2.1,0.8,"Structural", cex = 0.9, adj = c(0,0))
text(-2.1,0.6,"Unknown", cex = 0.9, adj = c(0,0))

#Export 6x6
