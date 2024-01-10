###########################################################################
## Project: CCCD
## Script purpose: create hazard prob and test vs sims
## Date: 02-03-2023
## Author: David JP O'Sullivan
###########################################################################

rm(list = ls()) # tidy work space
gc()

# libraries, source files and data ----------------------------------------

# source('./code/_project_setup.r')


library(tidyverse)
library(cowplot)
theme_set(cowplot::theme_cowplot())
# functions

library(purrr)
library(tidyverse)
library(igraph)



# Assuming prob_vector and exponents are complex vectors
map2_cplx <- function(prob_vector, exponents, pinf, x) {
  mapply(function(p, e) p * (1 - pinf + pinf * x) ^ e, prob_vector, exponents)
}


# Define the function to create a PGF
create_pgf <- function(prob_vector) {
  pgf <- function(s1, s2, s1_t, s2_t, sT, pinf) {
    exponents <- seq_along(prob_vector) - 1
    # result <- map2_dbl(prob_vector, exponents, function(p, e) p * (1 - pinf + pinf * x) ^ e)
    # Assuming prob_vector and exponents are of the same length
    result <- numeric(length(prob_vector))
    
    for (i in seq_along(prob_vector)) {
      p <- prob_vector[i]
      e <- exponents[i]
      result[i] <- p * (1 - pinf + pinf * s1 * s2 * s1_t * s2_t  * sT) ^ e
    }
    
    return(sum(result))
  }
  return(pgf)
}



# Define the function to compute the excess degree distribution
excess_degree_distribution <- function(degree_dist) {
  # k_values <- 0:(length(degree_dist) - 1)
  # exp_value <- sum(k_values * degree_dist)
  # 
  # numerators <- (k_values + 1)[-1] * degree_dist[-1]
  # excess_degrees <- numerators / exp_value
  
  # Assuming prob_vector is the degree distribution P(k)
  degree_values <- seq_along(degree_dist) - 1
  
  # Calculate the average degree <k>
  avg_degree <- sum(degree_values * degree_dist)
  
  # Calculate the excess degree distribution Q(k)
  excess_degree <- numeric(length(degree_dist) - 1)
  
  for (i in seq_along(excess_degree)) {
    k <- degree_values[i]
    excess_degree[i] <- (k + 1) * degree_dist[i + 1] / avg_degree
  }
  return(excess_degree) # Remove the first element (corresponding to k = 0)
}


n <- 50000
p <- 0.025/n
g <- igraph::ba.game(n, 1)
# g <- igraph::erdos.renyi.game(n, p = 0.01)
g_combined <- disjoint_union(g, g)

# Add edges between communities with probability p
# V1_indices <- V(g)
# V2_indices <- n + V(g)
# for (v1 in V1_indices) {
#   for (v2 in V2_indices) {
#     if (runif(1) < p) {
#       g_combined <- g_combined + edge(v1, v2)
#     }
#   }
#   if(v1 %% 250 == 0)print(v1)
# }

clusters(g_combined)
Xin <- degree_distribution(g)
Xin %>% plot(log = "xy")
Xin_t <- excess_degree_distribution(Xin)

Xin %>% sum
Xin_t %>% sum

g_Xin <- create_pgf(Xin)
g_Xin_t <- create_pgf(Xin_t)

# Xout <- prop.table(table(apply(get.adjacency(g_combined)[1:n, (n+1):(n+n)],1, sum)))
# Xout_t <- excess_degree_distribution(Xout)
# Xout %>% plot
Xout <- dpois(0:5, lambda = 2)/sum(dpois(0:5, lambda = 2))
Xout_t <- excess_degree_distribution(Xout)

g_Xout <- create_pgf(Xout)
g_Xout_t <- create_pgf(Xout_t)

g_T <- function(sT) sT
g_ic <- function(s1, s2, s1_t, s2_t, sT) s1 * sT

# create a function to iterate it -----------------------------------------



# pgf for the number of nodes at time t
G_N_t <- function(s1, s2, s1_t, s2_t, sT, t, pin, pout){
  if(t > 0){ # have to iterate the pgf at least one
    for(t_i in t:1){ # iterate the function t times
      
      new_s1 <- g_Xin(1, 1, s1_t, 1, sT, pin) * g_Xout_t(1, s2, 1, 1, sT, pout)
      new_s2 <- g_Xin(1, 1, 1, s2_t, sT, pin) * g_Xout_t(s1, 1, 1, 1, sT, pout)
      
      new_s1_t <- g_Xin_t(1, 1, s1_t, 1, sT, pin) * g_Xout(1, s2, s1_t, 1, sT, pout)
      new_s2_t <- g_Xin(1, 1, 1, s2_t, sT, pin) * g_Xout_t(1, 1, 1, s2_t, sT, pout)
      
      new_sT <- g_T(sT)
      
      # new_s1 <- g_Xout(s2, pinf) * g_Xout_t(s2_t, pinf)
      # new_s2 <- g_Xout(s1, pinf) * g_Xout_t(s1_t, pinf)
      # 
      # new_s1_t <- g_Xin(s1, pinf) * g_Xin_t(s1_t, pinf)
      # new_s2_t <- g_Xin(s2, pinf) * g_Xin_t(s2_t, pinf)
      
      # new_sT <- g_T(sT) * new_s1 * new_s2 * new_s1_t * new_s2_t
      
      s1 <- new_s1
      s2 <- new_s2
      
      s1_t <- new_s1_t
      s2_t <- new_s2_t
      
      sT <- new_sT
    }
  }

  ans = g_ic(s1, s2, s1_t, s2_t, sT) # apply initial conditions
  return(ans)
}



G_N_tV <- Vectorize(G_N_t)

invert_pgf_via_ifft <- function(t, pin, pout, M = 10^4){
  x <- exp(2*pi*1i*(0:(M-1))/M)
  
  
  pdf <- Re(fft(G_N_tV(s1 = 1, s2 = 1, s1_t = 1, s2_t = 1, sT = x, 
                       t, pin, pout) %>% Conj(),inverse = TRUE))
  pdf[pdf < 0] <- 0
  pdf <- pdf/sum(pdf)
  
  cas_dist <- tibble(cascade_size = 0:(length(pdf) - 1), prob = pdf) 
  
  return(cas_dist)
}
# G_N_tV(s1 = x[1:4], s2 = 1, s1_t = 1, s2_t = 1, sT = 1, t = 1, pinf = 0.1)

# invert_pgf_via_ifft(t = 1, pinf = .001, M = 50) %>% plot
# test --------------------------------------------------------------------


# G_N_t(0, 0, 0, 0, 1, 5, 0.1)
# G_N_t(1, 1, 1, 1, 1, 4, 0.001)

e_df <- tibble(t = 0:50) %>% 
  rowwise() %>% 
  mutate(
    e_prob_1 = G_N_t(0, 0, 0, 0, 1, t, 0.08, 0.08),
    e_prob_2 = G_N_t(0, 1, 0, 1, 1, t, 0.05, 0.0)
  ) %>% ungroup()

e_df %>% summarise(val=max(e_prob_1))

e_df %>% ggplot(aes(x = t, y = e_prob_1)) + 
  geom_point() + 
  ylim(c(0,1)) + 
  scale_y_log10() + 
  scale_x_log10()


e_df2 <- 
  bind_rows(
    invert_pgf_via_ifft(t = 30, pin = 0.08, pout = 0.08, M = 1000) %>% mutate(type = "1"), 
    invert_pgf_via_ifft(t = 30, pin = 0.08, pout = 0.05, M = 1000) %>% mutate(type = "2"), 
    invert_pgf_via_ifft(t = 30, pin = 0.08, pout = 0., M = 1000) %>% mutate(type = "3")
    ) 

e_df2 %>% group_by(type) %>% 
  summarise(sum(cascade_size * prob))

e_df2 %>% 
  group_by(type) %>% 
  mutate(
    ccdf = 1- cumsum(prob)
  ) %>% 
  ggplot(aes(x=cascade_size+1, y = prob, color = type)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_log10()


e_df2 %>% 
  group_by(type) %>% 
  mutate(
    ccdf = 1- cumsum(prob)
  ) %>% 
  ggplot(aes(x=cascade_size+1, y = ccdf, color = type)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_log10()


