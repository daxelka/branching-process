

###############################################################################################
## Project: cross community branching process
## Script purpose: example of inverting bp
## Date: 03/07/2022
## Author: David JP O'Sullivan
###############################################################################################

library(tidyverse)

# set up the branching processes' offspring distribution

# poisson distribution offspring
f_x <- function(s1,s2, l = 0.9) exp(l * (s1 * s2 - 1))
# type 2 (s2) can not die
f_T <- function(s) s
# initial conditions: one of each type
fic <- function(s1, s2) s1 * s2

# function to itterate each type
itt_gen <- function(s1, s2, lam = 0.9, n_max = 100) {
  # s1 and s2: dummy var of for each type
  # lam par for poisson dist
  # n_max: number of generation 
  for (n in n_max:1) { # iterate the generating
    new_s1 <- f_x(s1, s2, l = lam)
    new_s2 <- f_T(s2)
    
    s1 <- new_s1
    s2 <- new_s2
  }
  
  # apply initial conditions
  ans <- fic(s1, s2)
  
  return(ans)
}


# create point on the complex unit circle
# that are M point, equally spaced apart
M <- 100
x <- exp(2*pi*1i*(0:(M-1))/M)


plot(x)
itt_gen(1,x) # check results

# coefficints are mixed up
pdf <- Re(fft(itt_gen(1,x) %>% Conj(), inverse = TRUE))
pdf[pdf <= 0] <- 0 # replace any negatative coeff
pdf <- pdf / length(pdf) # require to sub to one
length(pdf) 
pdf %>% sum # check that it summs to one

# create the probabilties
dist_poi <- tibble(x = 0:99, prob = pdf)

dist_poi %>% ggplot(aes(x = x, y = prob)) + 
  geom_line() +
  geom_point() + 
  scale_x_log10()


# try a couple of different parameters ------------------------------------

# vec of lam
pos_vec <- seq(from = 0.1, to = 0.9, by = 0.2)
dist_poi_2 <- tibble() # create where to store the poisson

M <- 10000 # how many probability do we want? 
x <- exp(2 * pi * 1i * (0:(M - 1)) / M) 

for(i in seq_along(pos_vec)){
  
  # coefficints are mixed up
  pdf <- Re(fft(itt_gen(1, x, pos_vec[i], 1000) %>% Conj(), inverse = TRUE))
  pdf[pdf <= 0] <- 0 # replace any negatative coeff
  pdf <- pdf / length(pdf) # require to sub to one
  
  # create the probabilties
  dist_poi <- tibble(x = 0:(M-1), prob = pdf)
  dist_poi$lam <- pos_vec[i]
  
  dist_poi_2 <- bind_rows(dist_poi_2, dist_poi)
  
  print(glue::glue("Working on lam = {pos_vec[i]}, which is {i/length(pos_vec)}%"))
}

dist_poi_2 %>% filter(prob > 1e-16) %>% 
  ggplot(aes(x = x + 1, y = prob, color = factor(lam))) + 
  geom_line() + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10() + 
  scale_colour_viridis_d() + 
  theme_classic()
  
