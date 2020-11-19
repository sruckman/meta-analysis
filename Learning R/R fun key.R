#See where the working directory is
getwd()
#If not Learning R folder use setwd()
setwd("~/GitHub/meta-analysis/Learning R")

#Question 1 Dog Breeder t test
US <- c(12,12.5,11,10,13) 
UK <- c(12.5,13,12,12,11)
#H0: there is no difference between the corgi groups
#Ha: There is a difference between the corgi groups
t.test(US,UK)
#p-value = 0.5482
#Fail to reject
#There is insufficent evidence to conclude that the two
#groups of corgis have different heights

#Question 2 Treatment Effectiveness Chi-squared
#Read in the file
treat <- read.csv("treatment.csv", header =T)
#Add row names to the table
rownames(treat) <- treat[,1]
#Remove column 1 so you only have columns 2 and 3 (numeric data)
treat <- treat[,c(2,3)]
#H0: the drug is ineffective at treating migraines
#Ha: the drug is effective at treating migraines
chisq.test(treat)
#p-value = 0.03083
# Reject our null
#There is sufficient evidence to conclude that the drug
#is effective at treating migraines

#Question 3 Grasshopper paired t test
#Read in the data
grasshopper <- read.csv("grasshopper.csv", header=T)

#H0:there is no difference when the nerve is crushed
#Ha: there is a difference when the nerve is crushed
t.test(grasshopper[,2],grasshopper[,1], paired = T)
#p-value < 2.2 x 10^-16
#Reject null
#There is a difference in walking ability after the 
#nerve is crushed

#H0:there is no difference when the nerve is crushed (new nerve)
#Ha: there is a difference when the nerve is crushed (recovery)
t.test(grasshopper[,4],grasshopper[,3], paired = T)
#p-value = 1.103x10^-11
#Reject null
#There is a difference in walking ability after the nerve is recrushed
#This indicates that the grasshoppers recovered the nerve rather 
#than used a new nerve to walk
