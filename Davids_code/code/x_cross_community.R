##################################################################
## Project: what project are we working on?
## Script purpose: what does this code do?
## Date: what date is it?
## Author: David JP O'Sullivan
##################################################################

rm(list = ls()) # tidy work space
gc()

# libraries, source files and data ----------------------------------------

# want to calculate the cascade size distribution from a single clique
# pgf for each clique type


f_1 <- function(s1, s2, p = 0.2, l = 4) exp(l * (p * s1 * s2 - p))
f_2 <- function(s1, s2) s2

f_1(1,1)
f_2(1,1)

fic <- function(s1,s2) s1 * s2

multi_composite <- function(f, f1, f2) function(s1,s2) f(f1(s1,s2), f2(s1,s2))

## first method
# take a function and iterate n times
compose_itt_n <- function(fn, n) {
  if (n == 0) return(fn) # once its done stop and return the fn
  Recall(multi_composite(fn, f_1, f_2), n - 1) # call the function this is inside
  # n times each times apply the generation fn
}

# test how far we can push the cstack
compose_itt_n(fic, 1)(1,1)
compose_itt_n(fic, 900)(1,1)



invert_pgf_via_ifft <- function(gen_fn, M = 10^5){
  x <- exp(2*pi*1i*(0:(M-1))/M)
  
  pdf <- Re(fft(gen_fn(x) %>% Conj(),inverse = TRUE))
  pdf[pdf < 0] <- 0
  pdf <- pdf/length(pdf)
  
  cas_dist <- tibble(cascade_size = 0:(length(pdf) - 1), prob = pdf) 
  
  return(cas_dist)
}

library(tidyverse)

cas_after_100_gen <- function(x) compose_itt_n(fic, 900)(1,x)
cas_dist_100_gen <- invert_pgf_via_ifft(gen_fn = cas_after_100_gen, M = 10^6)

cas_dist_100_gen %>% summarise(prob = sum(prob))

cas_dist_100_gen %>% # mutate(prob = round(prob,10)) %>% 
  filter(prob != 0) %>%
  ggplot(aes(x = cascade_size, y = prob)) + 
  geom_line() +
  # xlim(c(1,1000)) +
  scale_x_log10() + 
  scale_y_log10()







# clean up work space --------------------------------------------------------
rm(list = ls()) # tidy work space
gc()