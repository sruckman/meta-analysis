#Question 1 Linear model
#HP vs LP aggression time based on sex
#read in the data 
fish <- read.csv("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis\\Learning R\\fish_agg.csv")

#H0: There is no difference based on sex or predation on agg time
#Ha: There is a difference based on sex or predation on agg time
#create linear model y = time, x = sex and predation
res <- lm(fish$Time ~ fish$Predation*fish$Sex)
summary(res) #See the results
#Since the p-value is less than 0.05 for predationLP, we can conclude
#that there is a significant difference in agg time between LP and HP

#Question 2 For Loops 
#Basics:
# for repeats a block of code between the curly braces
# each time that it runs i will take a new value
for(i in 1:10){
  print(i)
}

for(i in 20:40){
  x <- i / 5
  print(x)
}


#Answer to question 2

for(i in 1:10){
  x <- i/2
  print(x)
}

#Question 3 Correlation (cor.test)
#read in the data
walk <- read.csv("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis\\Learning R\\walk.csv")

#first plot the data
plot(walk$Time, walk$Temp, pch = 16, xlab = "Time", ylab = "Temp (F)")

#since the relationship looks linear use cor.test to do the regression
cor.test(walk$Time, walk$Temp)
#There is a strong negative correlation between time spent on the walk
#and temp (-0.97)

