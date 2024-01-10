###########################################################################
## Project: CCCD
## Script purpose: test mode results vs simulations
## Date: 02-03-2023
## Author: David JP O'Sullivan
###########################################################################

# Include the Rcpp package
library(Rcpp)

# Load the C++ implementation
sourceCpp("./code/cascade_estimation.cpp")

mean(degree(g_combined))
# mean(degree(g_simple))

1/mean(degree(g_combined)) # what should be the critical number for the simple case
# 1/mean(degree(g_simple))

2 * ecount(g_combined)/(vcount(g_combined) * (vcount(g_combined) - 1))
# 2 * ecount(g_simple)/(vcount(g_simple) * (vcount(g_simple) - 1))

M <- 1*10^6
pinf <- 0.1

## simulate on the network with the community structure
adj_comb <- as_adj(g_combined)
sim_res_com <- run_simulations(adj_comb, pinf, 0, 0, M) %>% as_tibble()
sum_res_com <- 
  sim_res_com %>% 
  count(num_active) %>% 
  mutate(
    cascade_size = num_active,
    prob = n/sum(n), 
    type = "sim com")

## simulate from a network with the same degree dists but with no community str
# Number of nodes in the graph
# n <- 5000
# # Calculate the degree sequence from the probability vector
# degree_sequence <- sample(0:(length(X) - 1), size = n, replace = TRUE, prob = X)
# # Create a graph with the sampled degree sequence
# is_graphical(degree_sequence)
# g_simple <- degree.sequence.game(degree_sequence, method = "simple.no.multiple")

## rewire the network to destroy the community structure
g_simple <- g_combined %>% rewire(keeping_degseq(niter = 10000))

## what does this look like
adj_simple <- as_adj(g_simple)
sim_res_simple <- run_simulations(adj_simple, pinf, 0, 0, M) %>% as_tibble()
sum_res_simple <- 
  sim_res_simple %>% 
  count(num_active) %>% 
  mutate(
    cascade_size = num_active,
    prob = n/sum(n), 
    type = "sim simple")

# group together both the results
sum_res <- bind_rows(sum_res_com, sum_res_simple)

sum_res %>% 
  ggplot(aes(x = cascade_size + 1, y = prob, color = type)) + 
  geom_line()+ 
  scale_x_log10() + 
  scale_y_log10()

# invert pgf's and check results
e_df3 <- 
  bind_rows(
    invert_pgf_via_ifft(t = 10, pin = pinf, pout = pinf, M = 1000) %>% mutate(type = "community"),
    # invert_pgf_via_ifft(t = 10, pin = 0.05, pout = 0.0, M = 1000) %>% mutate(type = "w comm 1 (off)"),
    invert_pgf_via_ifft_s(t = 10, pinf = pinf, M = 1000) %>% mutate(type = "no comm"),
    
  ) 

# cut pgf results off at the same point as the simulations
min_prob <- summarise(sum_res, min = min(prob))
e_df4 <- e_df3 %>% filter(prob >= min_prob$min)

# group them together
comp <- 
  bind_rows(
    e_df4, 
    sum_res
  ) %>% 
  group_by(type) %>% 
  mutate(
    prob = prob/sum(prob),
    ccdf = 1 - cumsum(prob)) %>% 
  ungroup()

# plots the results
# comp %>% 
#   ggplot(aes(x = cascade_size + 1, y = prob, color = type)) + 
#   geom_line()+ 
#   scale_x_log10() + 
#   scale_y_log10()

# comp %>%
#   ggplot(aes(x = cascade_size + 1, y = ccdf, color = type)) +
#   geom_line()+
#   scale_x_log10() +
# scale_y_log10()



# nice plots --------------------------------------------------------------

library(latex2exp)

e_df4$type %>% table
sum_res$type %>% table

theory_df <- e_df4 %>% 
  mutate(
    type = case_when(
      type == "community" ~ "With Communities (MTBP)",
      type == "no comm" ~ "No Communities (SBP)",
    )
  )

simulation_df <- sum_res %>% 
  mutate(
    type = case_when(
      type == "sim com" ~ "With Communities (MTBP)",
      type == "sim simple" ~ "No Communities (SBP)",
    )
  )

theory_df %>% mutate(theo = type) %>% 
  ggplot() + 
  geom_point(data = simulation_df, 
             aes(x = cascade_size + 1, y = prob, color = type), alpha = 0.5) +
  geom_line(aes(x = cascade_size + 1, y = prob, linetype = theo), size = 1.4) + 
  scale_x_log10() +
  scale_y_log10() + 
  xlab("Cascade size + 1") + 
  ylab("Prob of Cascade") + 
  labs(linetype = "Theory", color = "Simulation")
