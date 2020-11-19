#ANOVA Practice R Set
#this is the data that we are using
chickwts

#How to determine how many categories of feed there are
length(unique(chickwts$feed)) #6 levels

#Create the linear model
fit <- lm(chickwts$weight ~ chickwts$feed)
#Run the ANOVA
res <- aov(fit) #Need this form for Tukey Multiple Comparisons
#Can also use anova() but will not work with Tukey's

#See the results
summary(res)

#Since we have a significant result run multiple comparisons
TukeyHSD(res)

#summary(fit)
