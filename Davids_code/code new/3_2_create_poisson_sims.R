
rm(list = ls())

# Include the Rcpp package
library(tidyverse)
library(Rcpp)
library(igraph)

# Load the C++ implementation
sourceCpp("./code/cascade_estimation_new.cpp")

p_in <- 8/(2500)  # for the network size what should be the p_in and p_out for the SBM
p_out <- 2/(2500)
pinf <- 0.08
M <- 10

# create the matrix for how the communities are connected
pm <- matrix(c(p_in, p_out, p_out, p_in), byrow = TRUE, nrow = 2)

res_df <- tibble() # create where we are going to store the result

g_temp <- sample_sbm(5000, pm, c(2500,2500), directed = FALSE)
g_temp |> degree() |> mean()
adj_comb <- get.adjacency(g_temp)

sim_res_com <- run_simulations_cas(adj_comb, pinf, 1:2500, 0, 10^4,100) |> as_tibble()

sim_res_df <- sim_res_com |> 
  group_by(simulation_num, generation) |> 
  summarise(
    S1 = sum(child %in% 1:2500),
    S2 = sum(child %in% 2501:5000), 
    ST = S1 + S2
  ) |> 
  rename(
    t = generation
  ) |> 
  # add_row(.before = 1, t = 0, S1 = 1, S2 = 0, ST = 1)
  ungroup()



sim_summ_q <- sim_res_df |> #unnest(sim_list) |> 
  group_by(t) |> # for each generation t
  summarise(
    q_1_sim  = 1 - sum(S1 > 0)/M, # find the prob that it is dead for the different types
    q_2_sim  = 1 - sum(S2 > 0)/M,
    q_b_sim  = 1 - sum(ST > 0)/M,
    
    e_1_sim = sum(S1 == 0)/M,
    e_2_sim = sum(S2 == 0)/M,
    e_b_sim = sum(ST == 0)/M,
    
  )



ggplot(sum_res_com, aes(x = cascade_size + 1, y = prob)) + 
  geom_point() +
  scale_y_log10()



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

