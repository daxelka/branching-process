###########################################################################
## Project: CCCD
## Script purpose: create general mtbp for cccd
## Date: 02-03-2023
## Author: David JP O'Sullivan
###########################################################################

rm(list = ls()) # tidy work space
gc()

# libraries, source files and data ----------------------------------------

# source('./code/_project_setup.r')

# load libraries
library(tidyverse)
library(cowplot)
library(purrr)
library(tidyverse)
library(igraph)
library(Rcpp)
# set ploting theme
theme_set(cowplot::theme_cowplot())


# create pgf builder functions --------------------------------------------

# Define the function to create a PGF (Probability Generating Function) 
# given a vectors for probabilities
create_pgf <- function(prob_vector) { 
  pgf <- function(s, sT, pinf) { # return a function that takes these parameters
    exponents <- seq_along(prob_vector) - 1 # the powers of the pgf P[x=i] s ^ exponents
    result <- numeric(length(prob_vector))
    
    for (i in seq_along(prob_vector)) { # for every prob in prob_vectors
      p <- prob_vector[i] # get the p_i
      e <- exponents[i] # and the i for the p_i in the p_i s ^ i
      # calculate number of infected nodes which will give you a s type and a sT particle
      result[i] <- p * (1 - pinf + pinf * s * sT)^e
    }

    return(sum(result)) # sum the results
  }
  return(pgf)
}

# Define the function to compute the excess degree distribution
excess_degree_distribution <- function(degree_dist) { 

  # Assuming prob_vector is the degree distribution P(k)
  degree_values <- seq_along(degree_dist) - 1

  # Calculate the average degree E(k)
  avg_degree <- sum(degree_values * degree_dist)

  # Calculate the excess degree distribution q(k)
  excess_degree <- numeric(length(degree_dist) - 1)

  for (i in seq_along(excess_degree)) { # for each excess degree prob
    k <- degree_values[i] # grab the degree value
    excess_degree[i] <- (k + 1) * degree_dist[i + 1] / avg_degree
  }
  return(excess_degree) 
}


# create the network that we want to test ---------------------------------

n <- 5000 # number of nodes 
p_edge_between <- 0.05 / n # probability of edges between communities

g <- igraph::ba.game(n, 1, directed = F) # generate a ba network
(degree_distribution(g)) |> plot(log = "xy")

# now create the two communities, with no edges between them 
g_combined <- disjoint_union(g, g)
g_combined |>
  degree.distribution() |>
  plot(log = "xy")
# 
# # Add edges between communities with probability p
# V1_indices <- V(g)
# V2_indices <- n + V(g)
# for (v1 in V1_indices) {
#   for (v2 in V2_indices) {
#     if (runif(1) < p_edge_between) {
#       g_combined <- g_combined + edge(v1, v2)
#     }
#   }
#   if (v1 %% 250 == 0) print(v1)
# }

# for a large graph, there is a c++ function in the following file that will speed
# up the adding of these edges
sourceCpp("./code/cascade_estimation.cpp")
combined_adj <- as_adj(g_combined) # get the sparse adj matrix
# use the C++ function to add links between nodes
combined_adj <- add_edges_between_communities(combined_adj, n, p_edge_between)
# turn it back into a igraph
g_combined <- graph.adjacency(combined_adj, mode = "undirected")

(degree_distribution(g_combined)) |> plot(log = "xy")

# check if network is in one connected comp
clusters(g_combined)
# create the degree dist and excess degree dist for inside and outside the communities
Xin <- degree_distribution(g)
Xin_t <- excess_degree_distribution(Xin)

# Xout <- Xin/1.1
# Xout[1] <- 1 - sum(Xout)
# Xout_t <- excess_degree_distribution(Xout)

Xout <- prop.table(table(apply(get.adjacency(g_combined)[1:n, (n + 1):(n + n)], 1, sum)))
Xout_t <- excess_degree_distribution(Xout)

# create the pgfs for each
g_Xin <- create_pgf(Xin)
g_Xin_t <- create_pgf(Xin_t)

g_Xout <- create_pgf(Xout)
g_Xout_t <- create_pgf(Xout_t)

g_T <- function(sT) sT
g_ic <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin(s1, sT, pin) * g_Xout(s2, sT, pout) * sT #### <------ maybe this should be tilda s2_t

## orginal formulation of offspring function for differnet types
# g_1_in <- function(s1_in, s2_in, s1_out, s2_out, sT, pin, pout) g_Xin_t(s1_in, sT, pin) * g_Xout(s2_out, sT, pout)
# g_2_in <- function(s1_in, s2_in, s1_out, s2_out, sT, pin, pout) g_Xin_t(s2_in, sT, pin) * g_Xout(s1_out, sT, pout)
# g_1_out <- function(s1_in, s2_in, s1_out, s2_out, sT, pin, pout) g_Xin(s1_in, sT, pin) * g_Xout_t(s2_out, sT, pout)
# g_2_out <- function(s1_in, s2_in, s1_out, s2_out, sT, pin, pout) g_Xin(s2_in, sT, pin) * g_Xout_t(s1_out, sT, pout)

## this one seems to be the correct ones
g_1_in <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin_t(s1_t, sT, pin) * g_Xout(s2, sT, pout)
g_2_in <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin_t(s2_t, sT, pin) * g_Xout(s1, sT, pout)
g_1_out <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin(s1, sT, pin) * g_Xout_t(s2_t, sT, pout)
g_2_out <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin(s2, sT, pin) * g_Xout_t(s1_t, sT, pout)


# create a function to iterate it -----------------------------------------


# pgf for the number of nodes at time t
G_N_t <- function(s1, s2, s1_t, s2_t, sT, t, pin, pout) {
  if (t > 0) { # have to iterate the pgf at least one
    for (t_i in t:1) { # iterate the function t times
      
      # find the offspring of each type from
      # nodes that were arrived at inside each community
      new_s1 <- g_1_in(s1, s2, s1_t, s2_t, sT, pin, pout)
      new_s2 <- g_2_in(s1, s2, s1_t, s2_t, sT, pin, pout)

      # nodes that were arrived at from outside each communitiy
      new_s1_t <- g_1_out(s1, s2, s1_t, s2_t, sT, pin, pout)
      new_s2_t <- g_2_out(s1, s2, s1_t, s2_t, sT, pin, pout)
      
      # total number of nodes
      new_sT <- g_T(sT)
      
      # update values
      s1 <- new_s1
      s2 <- new_s2

      s1_t <- new_s1_t
      s2_t <- new_s2_t

      sT <- new_sT
    }
  }

  ans <- g_ic(s1, s2, s1_t, s2_t, sT, pin, pout) # apply initial conditions
  return(ans)
}
# create a vectorised version of the function
G_N_tV <- Vectorize(G_N_t)

# invert the pgf to get the total cascade size
invert_pgf_via_ifft <- function(t, pin, pout, M = 10^4) {
  # sample point on the complex unit cuircle
  x <- exp(2 * pi * 1i * (0:(M - 1)) / M) 
  pdf <- Re(fft(G_N_tV( # get the pdf via ifft
    s1 = 1, s2 = 1, s1_t = 1, s2_t = 1, sT = x,
    t, pin, pout
  ) |> Conj(), inverse = TRUE))
  # remove negative points and normalise
  pdf[pdf < 0] <- 0
  pdf <- pdf / sum(pdf)
  
  # create the cascade size
  cas_dist <- tibble(cascade_size = 0:(length(pdf) - 1), prob = pdf)
  return(cas_dist)
}
# G_N_tV(s1 = x[1:4], s2 = 1, s1_t = 1, s2_t = 1, sT = 1, t = 1, pinf = 0.1)

# plot an example
invert_pgf_via_ifft(t = 10, pin = 0.05, pout = .05, M = 1000) |>
  ggplot(aes(x = cascade_size, y = prob)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10()


# do the simple branching processes ---------------------------------------

X <- degree_distribution(g_combined)
X_t <- excess_degree_distribution(X)

g_X <- create_pgf(X)
g_X_t <- create_pgf(X_t)

# pgf for the number of nodes at time t
G_N_t_s <- function(s, sT, t, pinf) {
  if (t > 0) { # have to iterate the pgf at least one
    for (t_i in t:1) { # iterate the function t times

      new_s <- g_X_t(s, sT, pinf)
      new_sT <- g_T(sT)

      s <- new_s
      sT <- new_sT
    }
  }

  ans <- g_X(s, sT, pinf) * sT # apply initial conditions
  return(ans)
}

G_N_t_sV <- Vectorize(G_N_t_s)

invert_pgf_via_ifft_s <- function(t, pinf, M = 10^4) {
  x <- exp(2 * pi * 1i * (0:(M - 1)) / M)
  pdf <- Re(fft(G_N_t_sV(
    s = 1, sT = x,
    t, pinf
  ) |> Conj(), inverse = TRUE))
  pdf[pdf < 0] <- 0
  pdf <- pdf / sum(pdf)
  cas_dist <- tibble(cascade_size = 0:(length(pdf) - 1), prob = pdf)
  return(cas_dist)
}



# test for some values ----------------------------------------------------

# e_df <- tibble(t = 0:50) |>
#   rowwise() |>
#   mutate(
#     e_prob_1 = G_N_t(0, 0, 0, 0, 1, t, 0.1, 0.1),
#     e_prob_2 = G_N_t(0, 1, 0, 1, 1, t, 0.1, 0.0),
#     e_prob_s = G_N_t_s(0, 1, t, 0.1)
#   ) |>
#   ungroup()
# 
# e_df |> summarise(val = max(e_prob_1))
# 
# e_df |> ggplot(aes(x = t, y = e_prob_1)) +
#   geom_point() +
#   ylim(c(0, 1)) +
#   scale_y_log10() +
#   scale_x_log10()
# 
# 
# e_df2 <-
#   bind_rows(
#     invert_pgf_via_ifft(t = 10, pin = 0.1, pout = 0.1, M = 1000) |> mutate(type = "w comm 1", par = "pin=0.1"),
#     # invert_pgf_via_ifft(t = 10, pin = 0.1, pout = 0.0, M = 1000) |> mutate(type = "w comm 1 (off)", par = "pin=0.1"),
#     invert_pgf_via_ifft_s(t = 10, pinf = 0.1, M = 1000) |> mutate(type = "simple 1", par = "pin=0.1"),
#     
#     invert_pgf_via_ifft(t = 10, pin = 0.05, pout = 0.05, M = 1000) |> mutate(type = "w comm 2", par = "pin=0.05"),
#     # invert_pgf_via_ifft(t = 10, pin = 0.05, pout = 0.00, M = 1000) |> mutate(type = "w comm 2 (off)", par = "pin=0.05"),
#     invert_pgf_via_ifft_s(t = 10, pinf = 0.05, M = 1000) |> mutate(type = "simple 2", par = "pin=0.05")
#   )
# 
# e_df2 |>
#   group_by(type) |>
#   summarise(sum(cascade_size * prob))
# 
# e_df2 |>
#   group_by(type) |>
#   filter(prob >= 1e-15) |>
#   ggplot(aes(x = cascade_size + 1, y = prob, color = type)) +
#   # geom_point() +
#   geom_line() +
#   scale_y_log10() +
#   scale_x_log10() + 
#   facet_wrap(~par)

# e_df2 |>
#   group_by(type) |>
#   mutate(
#     ccdf = 1 - cumsum(prob)
#   ) |>
#   filter(prob >= 1e-14) |>
#   ggplot(aes(x = cascade_size + 1, y = ccdf, color = type)) +
#   geom_line() +
#   scale_y_log10() +
#   scale_x_log10()



# nice plots --------------------------------------------------------------

Xin <- c(3/4, 1/8, 1/8)
Xin_t <- excess_degree_distribution(Xin)

Xout <- c(3/4,1/4)
Xout_t <- excess_degree_distribution(Xout)

# create the pgfs for each
g_Xin <- create_pgf(Xin)
g_Xin_t <- create_pgf(Xin_t)

g_Xout <- create_pgf(Xout)
g_Xout_t <- create_pgf(Xout_t)

g_T <- function(sT) sT
g_ic <- function(s1, s2, s1_t, s2_t, sT, pin, pout) (s1^0)*g_Xout(s2_t, sT, 1)*g_Xin(s1_t, sT, 1)*(s2_t^0)*sT #### <------ maybe this should be tilda s2_t

s <- 0.5
g_ic(s1 = 0.5,s2 = 0.5,s1_t = 0.5, s2_t = 0.5,sT = 1)


## orginal formulation of offspring function for differnet types
# g_1_in <- function(s1_in, s2_in, s1_out, s2_out, sT, pin, pout) g_Xin_t(s1_in, sT, pin) * g_Xout(s2_out, sT, pout)
# g_2_in <- function(s1_in, s2_in, s1_out, s2_out, sT, pin, pout) g_Xin_t(s2_in, sT, pin) * g_Xout(s1_out, sT, pout)
# g_1_out <- function(s1_in, s2_in, s1_out, s2_out, sT, pin, pout) g_Xin(s1_in, sT, pin) * g_Xout_t(s2_out, sT, pout)
# g_2_out <- function(s1_in, s2_in, s1_out, s2_out, sT, pin, pout) g_Xin(s2_in, sT, pin) * g_Xout_t(s1_out, sT, pout)

## this one seems to be the correct ones
g_1_in <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin_t(s1_t, sT, pin) * g_Xout(s2, sT, pout)
g_2_in <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin_t(s2_t, sT, pin) * g_Xout(s1, sT, pout)
g_1_out <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin(s1_t, sT, pin) * g_Xout_t(s2_t, sT, pout)
g_2_out <- function(s1, s2, s1_t, s2_t, sT, pin, pout) g_Xin(s2_t, sT, pin) * g_Xout_t(s1_t, sT, pout)

G_N_t <- function(s1, s2, s1_t, s2_t, sT, t, pin, pout) {
  if (t > 0) { # have to iterate the pgf at least one
    for (t_i in 1:t) { # iterate the function t times
 
      # find the offspring of each type from
      # nodes that were arrived at inside each community
      new_s1 <- g_1_out(s1, s2, s1_t, s2_t, sT, pin, pout)
      new_s2 <- g_2_out(s1, s2, s1_t, s2_t, sT, pin, pout)
      
      # nodes that were arrived at from outside each community
      new_s1_t <- g_1_in(s1, s2, s1_t, s2_t, sT, pin, pout)
      new_s2_t <- g_2_in(s1, s2, s1_t, s2_t, sT, pin, pout)
      
      # total number of nodes
      new_sT <- g_T(sT)
      
      # update values
      s1 <- new_s1
      s2 <- new_s2
      
      s1_t <- new_s1_t
      s2_t <- new_s2_t
      
      sT <- new_sT
    }
  }
  
  ans <- g_ic(s1, s2, s1_t, s2_t, sT, pin, pout) # apply initial conditions
  return(ans)
}
# create a vectorised version of the function
G_N_tV <- Vectorize(G_N_t)

# invert the pgf to get the total cascade size
invert_pgf_via_ifft <- function(t, pin, pout, M = 20) {
  # sample point on the complex unit cuircle
  x <- exp(2 * pi * 1i * (0:(M - 1)) / M) 
  pdf <- Re(fft(G_N_tV( # get the pdf via ifft
    s1 = x, s2 = x, s1_t = x, s2_t = x, sT = 1,
    t, pin, pout
  ) |> Conj(), inverse = TRUE))
  # remove negative points and normalise
  pdf[pdf < 0] <- 0
  pdf <- pdf / sum(pdf)
  
  # create the cascade size
  cas_dist <- tibble(cascade_size = 0:(length(pdf) - 1), prob = pdf)
  return(cas_dist)
}
# G_N_tV(s1 = x[1:4], s2 = 1, s1_t = 1, s2_t = 1, sT = 1, t = 1, pinf = 0.1)

# plot an example
invert_pgf_via_ifft(t = 0, pin = 1, pout = 1, M = 10) |> mutate(prob = round(prob,10))
invert_pgf_via_ifft(t = 1, pin = 1, pout = 1, M = 10) |> mutate(prob = round(prob,10))
invert_pgf_via_ifft(t = 2, pin = 1, pout = 1, M = 20) |> mutate(prob = round(prob,15))
invert_pgf_via_ifft(t = 3, pin = 1, pout = 1, M = 20) |> mutate(prob = round(prob,15))


G_N_tV(s1 = 0, s2 = 0, s1_t = 0, s2_t = 0, sT = 1, t = 0, pin = 1, pout = 1)
G_N_tV(s1 = 0, s2 = 0, s1_t = 0, s2_t = 0, sT = 1, t = 1, pin = 1, pout = 1)
G_N_tV(s1 = 0, s2 = 0, s1_t = 0, s2_t = 0, sT = 1, t = 2, pin = 1, pout = 1)
G_N_tV(s1 = 0, s2 = 0, s1_t = 0, s2_t = 0, sT = 1, t = 3, pin = 1, pout = 1)
