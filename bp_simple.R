set.seed(1)
n_trials<-1000

# create offspring distribution
n_offspring <- 0:4
#offspring_dist <- c(0.55, 0.2, 0.2, 0, 0.05) # c() combines values into a list
offspring_dist <- c(0.35, 0.35, 0.1, 0.15, 0.05)

# average n_offspring
sum(offspring_dist * n_offspring)
data <- data.frame(trial = 1:n_trials,
                   extinct = 0,
                   generations = 0,
                   n_population=0)

# simulate branching process
# loops through n_trials

for (i in 1:n_trials){
  # Start with initial population
  population <- 1
  generation <- 0
  # generate offspring until extinct
  while(TRUE){
    population <- sum(sample(x=n_offspring,
                             size=population,
                             p=offspring_dist,
                             replace=TRUE))
    # extinct condition
    if (population == 0){
      data$extinct[i] <- 1
      data$generations[i] <-generation
      data$n_population[i] <- population
      break
    }
    
    # increment generations
    generation <- generation + 1
    
    if (population > 30){
      data$extinct[i] <- 0
      data$generations[i] <-generation
      data$n_population[i] <- population
      break
    }
  }
}

#visualise
ggplot(data, aes(x=generations)) + geom_histogram() + xlab('generations') + ggtitle('Branching Process')

# analysis
# see which generation lives the longest
data[which.max(data$generations),]
# checkout a generation
data[data$generations == 31,]

#estimate extinction probability
mean(data$extinct)
