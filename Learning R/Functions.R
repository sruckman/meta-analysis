#Functions
#Read in the data (surveys.csv)
surveys <- read.csv("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis\\Learning R\\surveys.csv")

#We have a collaborator that realizes there are tons of issues with the
#data including that in 1984 the weights and hindfoot lengths have been
#manipulated. Both variables have been multiplied by 1.1245697375083747 
#and 10 was added. Like before we could use loops:

#first let's copy the data and like last time
surveys_adjusted <- surveys

#Create our loop to save and change the variables
for (i in 1:dim(surveys)[1]) {
  if (surveys$year[i] == 1984) {
    surveys_adjusted$weight[i] <- surveys$weight[i]*1.1245697375083747+ 10
    surveys_adjusted$hindfoot_length[i] <- 
      surveys$hindfoot_length[i]*1.1245697375093747 + 10
  } 
}

#Do you see a problem with the code?


#I accidently typed the long number wrong and there is a 9 instead of 8
#Typos happen all the time and the best way to control them if you are
#needing to type the same thing multiple times is to use functions

#Functions
#Basic Structure:
#function_name <- function(arguements){
#   what you want it to do
#}

#Let's create a function that converts the 1984 data
#We will give it the value and have it convert the data
convert_1984 <- function(myval) {
  myval_adjusted <- myval*1.1245697375083747 + 10
  return(myval_adjusted)
}

#The function takes values (myval) and converts it by multiplying by 
#1.1245697375083747 and adding ten. Let's try 1 to see if it works
convert_1984(1)

#Now let's add our function to the loop so we can go through the data
for (i in 1:dim(surveys)[1]) {
  if (surveys$year[i] == 1984) {
    surveys_adjusted$weight[i] <- convert_1984(surveys$weight[i])
    surveys_adjusted$hindfoot_length[i] <- 
      convert_1984(surveys$hindfoot_length[i])
  } 
}

#Now if our collaborator comes back and says that we need to change
#anything all we have to do is change the function. Super easy!

############################ Your Turn! ##############################
#Your collaborator tells you that you can use the length of the hindfoot
#to calculate brain volume. Apparently the hindfoot of these creatures 
#is equal to the diameter of their skulls. Write a function that will 
#calculate the volume of the animals skulls and apply it to this dataset.
#Hint: the volume of a sphere is (4/3)*pi*r^3




#Create volume function with r as your parameter
vol <- function(r){
  (4/3)*pi*r^3
}

#create empty vectors to store the data
volume <- c()
radius <- c()

#Using a for loop find the volume of each hindfoot length
#First, we need to divide the lengths by 2 to get the radius
#Then, we can use our function to calculate the volume

for (i in 1:dim(surveys_adjusted)[1]) {
  radius[i] <- surveys_adjusted$hindfoot_length[i]/2
  volume[i] <- vol(radius[i])
}

#Let's visualize the data
plot(radius, volume, pch = 16, col = rgb(1,0,1,0.25))
