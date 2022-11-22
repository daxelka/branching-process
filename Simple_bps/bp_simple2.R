set.seed(0)
n_trials<-1000

# create offspring distribution
n_offspring <- 0:4
#offspring_dist <- c(0.55, 0.2, 0.2, 0, 0.05) # c() combines values into a list
offspring_dist <- c(0.35, 0.35, 0.1, 0.15, 0.05)
n_generations <-30

# average n_offspring
sum(offspring_dist * n_offspring)

# create data to store results
data <- data.frame(trial = 1:n_trials,
                   extinct = 0,
                   generations = 0,
                   n_population=0)

# create data to store simulation path
data_path <- matrix(data=0, nrow=n_generations, ncol=n_trials)
data_path[1,] <- 1


# simulate branching process
# loops through n_trials

for (i in 1:n_trials){
  # start with initial population
  population <- 1
  
  for (j in 2:n_generations){
    population <- sum(sample(x=n_offspring,
                             size=population,
                             p=offspring_dist,
                             replace=TRUE))
    
    data_path[j,i] <- population
    
    # extinct condition
    if (population == 0){
      data$extinct[i] <- 1
      data$generations[i] <-j
      data$n_population[i] <- population
      break
    }
  }
}

#visualise
data_plot <- data.table(generation =1:n_generations, data_path)
data_plot <- melt(data_plot, id.vars = "generation")
ggplot(data_plot, aes(x=generation, col=variable, y=value)) + 
  geom_line(alpha=1/10) + 
  xlab('generations') + ylab('generation size')+
  theme_bw() + theme(legend.position = "none")+
  ylim(c(0,100))+
  ggtitle('Branching Process')

# analysis
# numerical mean of generations
mean(data_path[10,])
# theoretical mean of generations
sum(offspring_dist*n_offspring)^9

# check extinction probability
mean(data$extinct)
mean(data$extinct[data_path[2,]==0])
mean(data$extinct[data_path[2,]==1])
mean(data$extinct[data_path[2,]==2])
mean(data$extinct[data_path[2,]==3])
mean(data$extinct[data_path[2,]==4])

mean(data$extinct[data_path[5,]==12])
sum(data_path[5,]==12)

# see which generation lives the longest
data[which.max(data$generations),]
# checkout a generation
data[data$generations == 31,]

#estimate extinction probability
mean(data$extinct)
