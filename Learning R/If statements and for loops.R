#For Loops and If statements
#Read in the data (surveys.csv)
survey <- read.csv("C:\\Users\\sarah\\Documents\\GitHub\\meta-analysis\\Learning R\\surveys.csv")

#Let's visualize the data using a histogram of the weights
hist(survey$weight, main = "Distribution of weights", 
     xlab = "weight (g)", col = "cyan1")
#Take the mean of weight and ignore NA values
mean(survey$weight, na.rm = T)

#Our weights are between 0 and 250 g which sounds about right for birds, 
#rabbits, rodents, or small reptiles
#Now let's adjust all of the weights up by 10% if they were taken in 1984
#First we need to decide if any weights were taken in 1984 or not 
#this can be done two ways: subsetting or if/else statements

#subset method
#Basic structure:
#name <- subset(data, data$column == or != "value")
#== means equals and != means does not equal
#Pull out rows by the year 1984
survey1984 <- subset(survey, survey$year == "1984")
#Pull out all rows other than the year 1984
surveyother <- subset(survey, survey$year != "1984")
survey1984$weight <- survey1984$weight*1.1
#recombine the data to create the adjusted data
adjsurvey <- rbind(survey1984,surveyother)

#Take the mean and see difference with original value
mean(adjsurvey$weight, na.rm = T)

#visualize the new data and compare it to the first
hist(adjsurvey$weight, main = "Distribution of adjusted weights", 
     xlab = "weight (g)", col = "orchid")

#if/else statements
#Basic structure:
#if (condition is TRUE) {
#  do something
#} else {
#  do different thing
#}
#Let's try this for the first value to figure out if it is 1984 or not
if (survey$year[1] == 1984){
  print("Great Scott, it's 1984")
} else {
  print("it's not 1984")
}

############################# Your turn! #################################
#Let's say we're interested in knowing whether an animal is large or not, 
#with a cut-off of at least one ounce. Write an if/else statement that 
#evaluates whether the 523rd animal in our data is larger than an ounce. 
#(Hint: one ounce is 28.3g) 

#Using an if/else statement
if (survey$weight[523] >= 28.3){
  print("You're gonna need a bigger boat")
} else {
  print("This boat will do")
}

#For Loops
#Basic Structure:
#for (variable in vector) {
#  do something
#}

#Create a for loop to print numbers 1 through 10
for (i in 1:10){
  print(i)
}

#You can also do this with words
for (i in c("corgi", "shark", "bantha")) {
  print(i)
}

#dimensions and prints rows by columns
dim(survey)

for (i in 1:dim(survey)[1]) {
  if (survey$year[i] == 1984) {
    print("Great Scott, it's 1984!")
  } else {
    print("It's not 1984.")
  }
}
#It prints either command for all values
#Let's have R adjust the weights
#Make a copy of survey
survey_adjusted <- survey
for (i in 1:dim(survey_adjusted)[1]) {
  if (survey_adjusted$year[i] == 1984) {
    print(survey_adjusted$weight[i]*1.1)
  } else {
    print("It's not 1984.")
  }
}
#You should see some numbers pop up but 
#there are a lot of points outside 1984
#Only print data from 1984 that is adjusted
for (i in 1:dim(survey_adjusted)[1]) {
  if (survey_adjusted$year[i] == 1984) {
    print(survey_adjusted$weight[i]*1.1)
  } 
}

#This only prints the values and we really want to replace the values
for (i in 1:dim(survey_adjusted)[1]) {
  if (survey_adjusted$year[i] == 1984) {
    survey_adjusted$weight[i] <- survey_adjusted$weight[i]*1.1
  } 
}

#Check to see if we got the same mean as before
mean(survey_adjusted$weight, na.rm = T)
#Awesome they match!

#Check only 1984 for each method
original_1984_weight <- mean(survey$weight[survey$year == 1984], 
                             na.rm = TRUE)
original_1984_weight
subset_1984_weight <- mean(adjsurvey$weight[adjsurvey$year == 1984], 
                           na.rm = T)
subset_1984_weight
ifelse_1984_weight <- mean(survey_adjusted$weight[survey_adjusted$year 
                                                  == 1984], na.rm = T)
ifelse_1984_weight
original_1984_weight * 1.1

############################### Your Turn! ##############################
#Using a for loop and an if/else statement, tally the number of animals 
#that weigh over an ounce in our adjusted dataset. To get you started, 
#here is code to create a data.frame where all recrods with NA for the 
#weight are removed:
#  surveys_adj_no_na <- surveys_adjusted[!is.na(surveys_adjusted$weight),]
#For the animals that are not over an ounce in weight, how many of them 
#are female and how many of them are male?

surveys_adj_no_na <- survey_adjusted[!is.na(survey_adjusted$weight),]
large <- data.frame(weight = rep(0,32283),sex = rep(0,32283))
small <- data.frame(weight = rep(0,32283),sex = rep(0,32283))
#This will take a while so let it sit until you get > back
for (i in 1:length(surveys_adj_no_na$weight)) {
  if (surveys_adj_no_na$weight[i] >= 28.3) {
    large$weight[i] <- surveys_adj_no_na$weight[i]
    large$sex[i] <- surveys_adj_no_na$sex[i]
  } else {
    small$weight[i] <- surveys_adj_no_na$weight[i]
    small$sex[i] <- surveys_adj_no_na$sex[i]
  }
} 
#Number of large animals removing the original zeros
large_no0 <- subset(large,large$weight != 0)
small_no0 <- subset(small, small$weight != 0)
dim(large_no0)

#Number of males and females with less than one ounce in weight
females <- data.frame(weight = rep(0,length(small_no0$weight)),
                      sex = rep(0,length(small_no0$weight)))
males <- data.frame(weight = rep(0,length(small_no0$weight)),
                    sex = rep(0,length(small_no0$weight)))
for (i in 1:length(small_no0$weight)) {
  if (small_no0$sex[i] == "F") {
    females$weight[i] <- small_no0$weight[i]
    females$sex[i] <- small_no0$sex[i]
  } else if (small_no0$sex[i] == "M") {
    males$weight[i] <- small_no0$weight[i]
    males$sex[i] <- small_no0$sex[i]
  } else {
    print("no sex")
  }
}
females_no0 <- subset(females,females$weight != 0)
males_no0 <- subset(males, males$weight != 0)
dim(females_no0)
dim(males_no0)

#Next week we will talk about functions