# Moment Conditions #
## This file cleans up the data to give the necessary moment conditions.

## Loading Data ##
data <- read.csv("~/Desktop/NBER Data/SIPP93/smallclean.csv")
library(plyr)
tic <- proc.time()

# Putting people into wage brackets. The cutoffs were roughly determined from
# the 1993 income tax code.
famearn <- data$hus_ear1 + data$wife_ear1
taxbracket <- matrix(0)
for (i in 1:dim(data)[1]){
  if (famearn[i] >= 20000){
    taxbracket[i] <- 1
  }
  else if (famearn[i] >= 11000 & famearn[i] < 20000){
    taxbracket[i] <- 2
  }
  else if (famearn[i] >= 5000 & famearn[i] < 11000){
    taxbracket[i] <- 3
  }
  else if (famearn[i] >= 3000 & famearn[i] < 5000){
    taxbracket[i] <- 4
  }
  else if (famearn[i] >= 0 & famearn[i] < 3000){
    taxbracket[i] <- 5
  }
}

# Putting people in age brackets. Group 1 is from 18-29, Group 2 is from 30-49
# and Group 3 is > 50
agebracket <- matrix(0)
for (i in 1:dim(data)[1]){
  if (age_01[i] >= 18 & age_01[i] < 30){
    agebracket[i] <- 1
  }
  else if (age_01[i] >= 30 & age_01[i] < 50){
    agebracket[i] <- 2
  }
  else if (age_01[i] >= 50){
    agebracket[i] <-3
  }
}

## (1) Childbirth rate conditional on parents'age/wage ##
# I find number of births by subtracting final number of kids from initial
# number of kids and then I summarize the number of births by wage and age
births <- data$nkids36 - data$nkids1
data <- cbind(data, births)
ddply(data, .(taxbracket, agebracket), summarize, 
      mean=mean(births), sd=sd(births), freq=length(births))

## (2) Quitting rate (of male and female, respectively) conditional on childbirth
# (say within three months), conditional on parents'age/wage (if possible. i.e. 
# if we have enough sample size)
# NOTE: This makes no distinction between the dad quitting or losing his job.

# Creating a variable which takes the following values:
# 0 if they do not have a kid
# 1 if they have a kid and do not quit
# 2 if they have a kid and quit
quitm <- matrix(0, nrow = dim(data)[1], ncol = 1)
quitf <- matrix(0, nrow = dim(data)[1], ncol = 1)
start_a <- grep("nkids", colnames(data))[1]  # First col of kids
start_b <- grep("hus_emp", colnames(data))[1]  # First col of husband employment
start_c <- grep("wife_emp", colnames(data))[1]  # First col of wife employment
for (i in 1:dim(data)[1]){
  for (j in 1:36){
    if (data[, start_a + j][i] == data[, start_a + (j - 1)][i] + 1){ # Had a kid
      if (data[, start_b + (j - 1)][i] == 1 &  # Dad has a job
          (data[, start_b + j][i] == 0 |
           data[, start_b + max((j + 1), 36)][i] == 0 |
           data[, start_b + max((j + 2), 36)][i] == 0 )){  # Dad quits
      quitm[i] <- 2
      }
      else if (data[, start_b + (j - 1)][i] == 1 &  # Dad has a job
          (data[, start_b + j][i] == 1 &
           data[, start_b + max((j + 1), 36)][i] == 1 &
           data[, start_b + max((j + 2), 36)][i] == 1 )){  # Dad does not quit
      quitm[i] <- 1
      }
    }
  }
}

for (i in 1:dim(data)[1]){
  for (j in 1:36){
    if (data[, start_a + j][i] == data[, start_a + (j - 1)][i] + 1){ # Had a kid
      if (data[, start_c + (j - 1)][i] == 1 &  # Mom has a job
          (data[, start_c + j][i] == 0 |
           data[, start_c + max((j + 1), 36)][i] == 0 |
           data[, start_c + max((j + 2), 36)][i] == 0 )){  # Mom quits
      quitf[i] <- 2
      }
      else if (data[, start_c + (j - 1)][i] == 1 &  # Mom has a job
          (data[, start_c + j][i] == 1 &
           data[, start_c + max((j + 1), 36)][i] == 1 &
           data[, start_c + max((j + 2), 36)][i] == 1 )){  # Mom does not quit
      quitf[i] <- 1
      }
    }
  }
}

data <- cbind(data, quitm, quitf)

# Getting quitting rate conditional on having a kid
sum(quitm == 2) / (sum(quitm == 1) + sum(quitm == 2))
sum(quitf == 2) / (sum(quitf == 1) + sum(quitf == 2))

# (3) Mean leaving duration conditional on (2), and on parents'age/wage (again, 
# if possible)

# Getting subset of couples that have at least one person quitting upon having
# a kid.
quitters <- rbind(data[quitm == 2, ], data[quitf == 2, ])
quittersm <- quitters[quitters$quitm == 2, ]
quittersf <- quitters[quitters$quitf == 2, ]

durationm <- matrix(0, nrow = dim(quittersm)[1], ncol = 1)
durationf <- matrix(0, nrow = dim(quittersf)[1], ncol = 1)

# Going through and finding the period when they have a kid and then finding the
# earliest time that either the mom or the dad comes back to work.
for (i in 1:dim(quittersm)[1]){
  for (j in 1:36){
    if(quittersm[, start_a + j][i] == quittersm[, start_a + (j - 1)][i] + 1){
      durationm[i] <- min(which(quittersm[i, (start_b + j - 1):(start_b + 35)] == 1), 
                         36 - j)
    }
  }
}

for (i in 1:dim(quittersf)[1]){
  for (j in 1:36){
    if(quittersf[, start_a + j][i] == quittersf[, start_a + (j - 1)][i] + 1){
      durationf[i] <- min(which(quittersf[i, (start_c + j - 1):(start_c + 35)] == 1), 
                         36 - j)
    }
  }
}

mean(durationm)
mean(durationf)

# (4) Wage difference before and after the leave (male and female)
wagedifm <- matrix(0, nrow = dim(quittersm)[1], ncol = 1)
wagediff <- matrix(0, nrow = dim(quittersf)[1], ncol = 1)

start_d <- grep("hus_ear", colnames(data))[1]  # First col of husband earnings
start_e <- grep("wife_ear", colnames(data))[1]  # First col of wife earnings

for (i in 1:dim(quittersm)[1]){
  for (j in 1:36){
    if(quittersm[, start_a + j][i] == quittersm[, start_a + (j - 1)][i] + 1){
      startwage <- quittersm[, start_d + j - 1][i]
      returnwage <- min(quittersm[i, (start_d + j - 1):(start_d + 35)] * 
                        (quittersm[i, (start_d + j - 1):(start_d + 35)] > 0))
      wagedifm[i] <- returnwage - startwage
    }
  }
}

for (i in 1:dim(quittersf)[1]){
  for (j in 1:36){
    if(quittersf[, start_a + j][i] == quittersf[, start_a + (j - 1)][i] + 1){
      startwage <- quittersf[, start_e + j - 1][i]
      returnwage <- min(quittersf[i, (start_e + j - 1):(start_e + 35)] * 
                        (quittersf[i, (start_e + j - 1):(start_e + 35)] > 0))
      wagediff[i] <- returnwage - startwage
    }
  }
}

mean(wagedifm)
mean(wagediff)

# (5) Employment rate (of male and female respectively) conditional on parents' 
# age/youngest child's age (assume parents' age are the same)
# (6) Mean and standard deviation of wages (both male and female) conditional on 
# parents' age
# (7) Transition probability from ee to uu
# (8) Transition probability from eu to uu, and ee to ue without childbirth 
# (9) Transition probability from ue to uu, and ee to eu without childbirth
# (10) Transition probability from uu to ue, and eu to ee
# (11) Transition probability from uu to eu, and ue to ee
toc <- proc.time()
print(toc - tic)
