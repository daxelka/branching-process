###############################################################################################
## Project: cross community cascade diffusion
## Script purpose: create the multitype branching processes sim code
## Date: 01-05-2022
## Author: David JP O'Sullivan
###############################################################################################

rm(list = ls()) # tidy work space
gc()

# libraries, source files and data ---------------------------------------------

# load in packages and codde 
# source('./code/_project_setup.r')
source('./code/_code.R')
library(tidyverse) # load required packages


# bp simulations ----------------------------------------------------------


prob_ex_pois <- function(k, lambda)((k + 1)/ lambda) * dpois(k + 1, lambda)
# create the excess degree distrubion for a poi random var (we actually dont need
# this)

dpois(0:10, 9) %>% round(2) # see what the prob dist looks like

# create a function to sample n iid from the excess dist
rex_pois <- function(n, lambda) sample(x = 0:100, n, replace = T, prob = prob_ex_pois(0:100, lambda))

# create a function for the offspring distrubtion that gives the number of 
# of infections if we have  rpois(n, lambda) links 
g1 <- function(n, lambda, p) rbinom(1, rpois(n, lambda), p)
g <- function(n, lambda, p) rbinom(1, rex_pois(n, lambda), p)

# the excess degree dist for a sbm is also poi with the same parameter
# but its handy to build in the logic if we want to change it later

# multitype code ----------------------------------------------------------


r_mtbp_sim <- function(max_g = 5000, lambda_in = 6, lambda_out = 6, 
                       p_in = 0.05, p_out = 0.05){
  # create the data frame to hold the results of the sim
  sim_results <- tibble(n_gen = 0:max_g, n_s1 = NA, n_s2 = NA)
  # what generation do we start at; initial conditions for s1 and
  # s2 types
  n_gen <- 0; new_s1 <- 1; new_s2 <- 0;
  
  # do the simulation by:
  # first saving the initial coniditions
  sim_results$n_s1[n_gen + 1] <- new_s1
  sim_results$n_s2[n_gen + 1] <- new_s2
  
  repeat{ # and repeating the following processes 
    # until there are no new types or may gen is reached
    if((new_s1 + new_s2 == 0) | n_gen == max_g) break
    
    # save the number of nodes of each type in a temp var
    s1_temp <- new_s1
    s2_temp <- new_s2
    
    # if we are not beyond the seed generation
    if(n_gen == 0){
      
      # simuate the number of new s1 and s2 ndoes by
      # finding the number of of s1 of sprind that all s1 and s2 nodes have had
      new_s1 <- sum(g1(s1_temp, lambda_in, p_in), na.rm = TRUE) + sum(g1(s2_temp, lambda_out, p_out), na.rm = TRUE)
      # same as above but for s2 offspring that s1 and s2 have had
      new_s2 <- sum(g1(s1_temp, lambda_out, p_out), na.rm = TRUE) + sum(g1(s2_temp, lambda_in, p_in), na.rm = TRUE)
    } else { # but same as above but using the excess degree dist 
      new_s1 <- sum(g(s1_temp, lambda_in, p_in), na.rm = TRUE) + sum(g(s2_temp, lambda_out, p_out), na.rm = TRUE)
      new_s2 <- sum(g(s1_temp, lambda_out, p_out), na.rm = TRUE) + sum(g(s2_temp, lambda_in, p_in), na.rm = TRUE)
    }
    
    
    # save the results
    n_gen <- n_gen + 1
    sim_results$n_s1[n_gen + 1] <- new_s1
    sim_results$n_s2[n_gen + 1] <- new_s2
  }
  
  # once we have finished the simulation
  sim_results <- sim_results %>% # clean the results by
    filter(complete.cases(.)) %>% # making sure there is no NA values
    mutate( # adding col to the data
      # that are the total number of s1, s2 and over all cascades sizes
      t_s1 = sum(n_s1), t_s2 = sum(n_s2), t_s = t_s1 + t_s2,
      # along with the proportion of different types in the 
      # cascade
      p_s1 = t_s1/t_s, p_s2 = 1 - p_s1) %>% 
    slice(1:(n()-1))
  return(sim_results) # return the data frame
}


sim_res <- # run the above code by:
  tibble(
    sim = 1:5000, # setting up the number of monte carlo simulations 
    sim_res = list(NULL) # and a list to store the resutls of each sim
    ) %>%
  group_by(sim) %>% # for each sim id
  mutate(sim_res = list(r_mtbp_sim())) %>% # run the multi type bp sim and save the results
  unnest(sim_res) %>% 
  slice(n()) %>% 
  ungroup()


# plot some results

# sim_res %>% count(t_s) %>% mutate(p = n/sum(n)) %>%
#   ggplot(aes(x = t_s, y = p)) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10()
# 

sim_res %>% 
  # filter(p_s1 != 1) %>%
  # filter(t_s > 1) %>%
  ggplot(aes(x = p_s1)) + 
  geom_histogram(color = 'black', fill = 'gold', bins = 100)

sim_res %>% 
  select(t_s, t_s1, t_s2, p_s1) %>%
  mutate(p_s2 = 1 - p_s1) %>% 
  as.matrix() %>% 
  cov()

sim_res %>% # filter(n_gen > 0) %>%
  # select(t_s,t_s1, t_s2, p_s1) %>%
  mutate(p_s2 = 1 - p_s1) %>% 
  summarise(
    max(n_gen),
    mt_s = mean(t_s), mt_s1 = mean(t_s1), mt_s2 = mean(t_s2),
    p_s1 = mean(p_s1), p_s2 = mean(p_s2),
    cov12 = cov(t_s1,t_s2), var1 = var(t_s1), var2 = var(t_s2),
    cov1 = cov(t_s,t_s1), vart = var(t_s),
    ) %>% 
  mutate(
    est_s1moment = mt_s1/mt_s,
    est_2 = est_s1moment - (cov1/(mt_s^2)) + (vart * mt_s1)/(mt_s^3)
  ) %>%
  mutate(
    exp_ts12 = mt_s1/(mt_s1 + mt_s2) + (2*mt_s1*cov12)/(mt_s1 + mt_s2)^3 - cov12/(mt_s1 + mt_s2)^2 + (mt_s1*var1)/(mt_s1 + mt_s2)^3 - var1/(mt_s1 + mt_s2)^2 + (mt_s1*var2)/(mt_s1 + mt_s2)^3
  )

0.964/.881
# -------------------------------------------------------------------------

# the following code does a network base simulation of the above branching processes

library(igraph)

p_in <- 9/(1000)  # for the network size what should be the p_in and p_out for the SBM
p_out <- 4/(1000)
# create the matrix for how the communities are connected
pm <- matrix(c(p_in, p_out, p_out, p_in), byrow = TRUE, nrow = 2)

res_df <- tibble() # create where we are going to store the result
for(i in 1:5000){ # for each i in the number sim we want to run
  # sample a sbm network
  
  g_temp <- sample_sbm(2000, pm, c(1000,1000), directed = FALSE)
  
  # g_temp <- erdos.renyi.game(1000, p_in)
  # try({
  # g_temp <- sample_degseq(out.deg = rpois(n = 10000, lambda = 9))
  # g_temp <- k.regular.game(1000,k = 6)
  # degree(g_temp) %>% mean
  # degree_distribution(g_temp)
  # cluster.distribution(g_temp)
  
  # give the nodes names (the sim code need this)
  V(g_temp)$name <- 1:vcount(g_temp) %>% as.character()
  seed <- sample(V(g_temp)$name[1:1000], 1) # make sure the sim always selects 
  # a node in s1
  
  # then generate the sim by
  temp <- generate_cc_cascades(
    # using the sbm network we created
    follower.net = g_temp, follower.adj = get.adjacency(g_temp), 
    # using this probability of not adopting (1 - prob of adopt)
    q1 = 1 - 0.05, alpha = 0, total = 1, seed_ = seed) %>% 
    mutate(ID = i) # add a col that has the sim id
  
  # if at least one link was created add to results
  if(nrow(temp) > 0) res_df <- bind_rows(res_df, temp)
  
  # },silent = TRUE)
  
}

# once we have all the simulations from the network get some summary 
# statistics by 
comp_df <- 
  res_df %>% # taking the simulated results
  group_by(ID) %>% # and for each sim id 
  summarise(
    # find number of s1 childeren (they will have names in 1:1000)
    child_n_s1 = sum(child %in% 1:1000),
    # parent_n_s1 = sum(parent[1] %in% 1:1000),
    child_n_s2 = sum(!(child %in% 1:1000)),
    # parent_n_s2 = sum(!(parent[1] %in% 1:1000))
    n_gen = n()
  ) %>% 
  mutate( # calculate the same summary statistics 
    # as we had for the bp simulations
    t_s1 = child_n_s1 + 1,
    t_s2 = child_n_s2,
    t_s = t_s1 + t_s2, 
    p_s1 = t_s1/t_s, p_s2 = 1 - p_s1
  )

# comp_df %>% 
#   # filter(n_gen > 1) %>%
#   ggplot(aes(x = p_s1)) + 
#   geom_histogram(color = 'black', fill = 'gold', bins = 100)
# 
# comp_df %>% 
#   select(t_s1, t_s2, p_s1, p_s2) %>%
#   as.matrix() %>% 
#   cov()
# 
# 
# comp_df %>% count(t_s) %>% mutate(p = n/sum(n)) %>% 
#   ggplot(aes(x = t_s, y = p)) + 
#   geom_point() + 
#   scale_y_log10()

# create some overall simmaries
comp_df %>% # filter(n_gen == 1) %>% 
  summarise(
    mean(p_s1), mean(p_s2),
    mean(child_n_s1), 
    mean(t_s1), mean(t_s2), mean(t_s),
    mean(t_s1), mean(t_s2), mean(t_s))



# plot a couple of things. 

temp_1 <- sim_res %>% 
  filter(n_gen != 0) %>% count(t_s) %>% mutate(p = n/sum(n), type = 'bp')
temp_2 <- comp_df %>% count(t_s) %>% mutate(p = n/sum(n), type = 'net') 

temp <- bind_rows(temp_1, temp_2) %>% 
  group_by(type) %>% 
  mutate(ccdf = cumsum(p)) %>% ungroup()

temp %>% 
  ggplot(aes(x = t_s, y = p, color = type)) + 
  geom_point() +
  scale_y_log10()



# -------------------------------------------------------------------------

temp_1 <- sim_res %>% filter(n_gen != 0) %>% mutate(type = 'bp')
temp_2 <- comp_df %>% mutate(type = 'net') 

temp <- bind_rows(temp_1, temp_2) # %>% 
  group_by(type) %>% 
  mutate(ccdf = cumsum(p)) %>% ungroup()

temp %>% 
  ggplot(aes(x = p_s1, color = type, fill = type)) + 
  geom_histogram(position = "dodge", alpha = 0.8, bins = 100) 


ggplot(temp, aes(x = p_s1, fill = type)) + 
  geom_density(color = 'black', alpha = 0.2, bw = 0.01) 


# -------------------------------------------------------------------------


### ignore this down, just wanted to calculate a couple of things


library(mnormt)

cov_data <- sim_res %>% 
  select(t_s1, t_s2, p_s1) %>%
  mutate(p_s2 = 1 - p_s1) %>% 
  as.matrix() %>% 
  cov(); sigma <- cov_data[3:4,3:4] %>% as.matrix()

mu <- 
  sim_res %>% filter(n_gen > 0) %>% 
  select(t_s1, t_s2, p_s1) %>%
  mutate(p_s2 = 1 - p_s1) %>% 
  summarise(mean(p_s1), mean(p_s2)) %>% 
  as.numeric()

cov_df <- 
  expand.grid(p_s1 = seq(0, 1, 0.01), p_s2 = seq(0, 1, 0.01)) %>% 
  as_tibble %>% rowwise() %>% 
  mutate(z = dmnorm(x = c(p_s1, p_s2), mean = mu, varcov = sigma + 0.0001))

ggplot(cov_df, aes(x = p_s1, y = p_s2, z = z)) + 
  geom_contour_filled() + 
  geom_contour()




# simulate processes ------------------------------------------------------
M <- 100
res_m_df <- tibble(draw = 1:M, m_p_s1 = NA)
for(i in 1:M){
  sim_res <-
    tibble(sim = 1:500, sim_res = list(NULL)) %>%
    group_by(sim) %>%
    mutate(sim_res = list(r_mtbp_sim())) %>%
    unnest(sim_res) %>%
    slice(n()) %>% 
    ungroup()
  
  # sim_res %>% count(t_s) %>% mutate(p = n/sum(n)) %>%
  #   ggplot(aes(x = t_s, y = p)) +
  #   geom_point() +
  #   scale_x_log10() +
  #   scale_y_log10()
  # 
  
  
  res_m_df$m_p_s1[i] <- sim_res %>% filter(n_gen > 0) %>% 
    select(t_s1, t_s2, p_s1) %>%
    mutate(p_s2 = 1 - p_s1) %>% 
    summarise(m_p_s1 = mean(p_s1), mean(p_s2)) %>% 
    pull(m_p_s1)
  print(i)
}
  

res_m_df %>% ggplot(aes(x = m_p_s1)) + 
  geom_density(fill = 'gold', color = 'black', bw = 0.0025) + 
  stat_function(fun = dnorm, args = list(mean = mu[1], sd = sigma[1,1])) +
  xlim(c(0.7, 0.9))

res_m_df %>% ggplot(aes(x = m_p_s1)) + 
  geom_histogram(fill = 'gold', color = 'black')



# cascade size ------------------------------------------------------------

mu1 <- 5 * 0.05
mu2 <- 3 * 0.05

M = matrix(
  ncol = 4, byrow = TRUE, 
  c(mu1, mu2, 1,1,
    mu2, mu1, 1,0,
    0, 0, 1, 0,
    0, 0, 0, 1)
  )

res <- c(1,0,0,0) %*% round(M %^% 2, 10)

res[1,4]/res[1,3]

1./1.65