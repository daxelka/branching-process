
# r code ------------------------------------------------------------------


library(tidyverse)
library(cowplot)
theme_set(cowplot::theme_cowplot())
# functions

# offspring functions for type 1 and type 2 with poisson 
# degree dist
g_X_1 <- function(s1, s2, lin, lout, pin) {
  ans <- exp(lin * pin * (s1-1)) * exp(lout * pin * (s2-1))
  return(ans)
}

g_X_2 <- function(s1, s2, lin, lout, pin) {
  ans <- exp(lin * pin * (s2-1)) * exp(lout * pin * (s1-1))
  return(ans)
}

# initial conditions, we always start with a seed from com 1
gic_1 <- function(s1, s2) s1

# pgf for the number of nodes at time t
G_N_t <- function(s1, s2, t, lin, lout, pin){
  if(t > 0){ # have to iterate the pgf at least one
    for(t_i in t:1){ # iterate the function t times
      new_s1 <- g_X_1(s1, s2, lin, lout, pin)
      new_s2 <- g_X_2(s1, s2, lin, lout, pin)
      s1 <- new_s1
      s2 <- new_s2
    }
  }
  
  ans = gic_1(s1, s2) # apply initial conditions
  return(ans)
}

# pgf for for N(t) > 0
G_Y_t <- function(s1, s2, t, lin = 8, lout = 2, pin = 0.05){
  if(t>=0){
    # find the prob of process ending at t
    G_N_t_num <- G_N_t(0,0,t, lin, lout, pin)
    
    # create pgf for G_Y_t by subtracting G_N_t_num and normalising it
    ans <- ((G_N_t(s1, s2, t, lin, lout, pin) - G_N_t_num) / (1 - G_N_t_num))
  } else{
    ans <- 0
  }
  return(ans)
}

# pgf for P(N_1(t) = k_1, N_2(t) = k_2 | N(t-1) > 0)
G_H_t <- function(s1, s2, t, lin = 8, lout = 2, pin = 0.05){
  if(t>=0){
    
    G_Y_t_num <- G_Y_t( # iterate G_Y_t again t - 1 times then
      # roll the number of offspring forward one generation
      s1 = g_X_1(s1, s2, lin, lout, pin), 
      s2 = g_X_2(s1, s2, lin, lout, pin),
      t-1, lin, lout, pin)
    ans <- G_Y_t_num
  } else{
    ans <- 0
  }
  return(ans)
}

G_C1_t <- function(s1 = 0, s2 = 1, t, lin = 8, lout = 2, pin = 0.05){
  if(t == 0){
    ans <- NA
  } else if(t>=1) {
    ans <- G_N_t( # iterate G_Y_t again t - 1 times then
      # roll the number of offspring forward one generation
      s1 = s1,
      s2 = g_X_2(s1, s2, lin, lout, pin),
      t-1, lin, lout, pin)
  } else{
    ans <- 0
  }
  return(ans)
}

G_C2_t <- function(s1 = 1, s2 = 0, t, lin = 8, lout = 2, pin = 0.05){
  if(t == 0){
    ans <- NA
  } else if(t>=1) {
    ans <- G_N_t( # iterate G_Y_t again t - 1 times then
      # roll the number of offspring forward one generation
      s1 = g_X_1(s1, s2, lin, lout, pin),
      s2 = s2,
      t-1, lin, lout, pin)
  } else{
    ans <- 0
  }
  return(ans)
}

# create hazards ----------------------------------------------------------

pin <- 0.09
lin <- 8
lout <- 2
(lout + lin) * pin


full_sur_dis_df <- tibble(t = 0:10) %>% 
  rowwise() %>%
  mutate(
    
    hazard_1 = G_H_t(0,1,t,lin, lout, pin),
    hazard_2 = G_H_t(1,0,t,lin, lout, pin),
    hazard_b = G_H_t(0,0,t,lin, lout, pin),
    
    q_1 = G_N_t(0,1,t,lin, lout, pin),
    q_2 = G_N_t(1,0,t,lin, lout, pin),
    q_b = G_N_t(0,0,t,lin, lout, pin),
    
  ) %>% 
  ungroup() %>% 
  mutate(
    
    e_1 = 0,
    e_2 = c(1,rep(0,n()-1)),
    ce_1 = 0,
    ce_2 = c(1,rep(0,n()-1))
  )

for(i in 2:nrow(full_sur_dis_df)){
  
  etm1_1 <- full_sur_dis_df$e_1[i-1]
  hp_1 <- full_sur_dis_df$hazard_1[i]
  
  ce_1 <- G_C1_t(t = full_sur_dis_df$t[i],lin = lin, lout = lout, pin = pin)
  full_sur_dis_df$ce_1[i] <- ce_1
  full_sur_dis_df$e_1[i] <- ce_1 * etm1_1 + hp_1 * (1 - etm1_1)
  
  
  etm1_2 <- full_sur_dis_df$e_2[i-1]
  hp_2 <- full_sur_dis_df$hazard_2[i]
  
  ce_2 <- G_C2_t(t = full_sur_dis_df$t[i],lin = lin, lout = lout, pin = pin)
  full_sur_dis_df$ce_2[i] <- ce_2
  full_sur_dis_df$e_2[i] <- ce_2 * etm1_2 + hp_2 * (1 - etm1_2)
}


# simulation code ---------------------------------------------------------


# function to sim bp
bp_sim <- function(Z_mat, lin = 8, lout = 2, 
                   pin = 0.09, pout = 0.09){
  # initial seed
  Z_mat[1,] <-c(1,0) # save the results the matrix
  for(i in 2:nrow(Z_mat)){ # simulate the process
    Z_new <- c(0,0)
    Z_old <- Z_mat[i-1,]
    # generate offspring from type 1
    Z_new[1] <- rpois(Z_old[1], lin * pin) %>% sum
    Z_new[2] <- rpois(Z_old[1], lout * pout) %>% sum
    
    # generate offspring from type 2 and add to the total
    Z_new[2] <- Z_new[2] + rpois(Z_old[2], lin * pin) %>% sum
    Z_new[1] <- Z_new[1] + rpois(Z_old[2], lout * pout) %>% sum
    
    Z_mat[i,] <- Z_new
    if(Z_new[1] == 0 & Z_new[2] == 0) break
  }
  
  return(Z_mat)
  # return(Z_mat[1:i,])
}

# max 5 gen
res_mat <- matrix(0, ncol = 2, nrow = 10)

no_sims <- 2500
sim_res_df <- tibble(sim = 1:no_sims, sim_list = list(NULL))


for(i in 1:nrow(sim_res_df)){
  res <- 
    bp_sim(res_mat, lin, lout, pin = pin, pout = pin) %>%  # take the emtpy matrix, add the sims
    `colnames<-`(c('S1','S2')) %>% # give the results col names
    as_tibble() %>%  # convert to a tibble
    mutate(t = 0:(n()-1), ST = S1 + S2) %>% # add a generation col, and total
    select(t, everything()) # reorder cols
  
  
  sim_res_df$sim_list[[i]] <- res
  
  if(i %% 500 == 0) print(glue::glue("Currently finished {i} of {nrow(sim_res_df)} ({(i/nrow(sim_res_df)) * 100}%)."))
  
}



full_sur_dis_df

sim_summ_e <-
  sim_res_df %>% unnest(sim_list) %>%
  group_by(sim) %>%
  mutate(
    lag_1 = lag(S1), lag_2 = lag(S2),
  ) %>% ungroup() %>%
  group_by(t) %>% # for each generation t
  summarise(
    # q_1_sim  = 1 - sum(S1 > 0)/no_sims, # find the prob that it is dead for the different types
    # q_2_sim  = 1 - sum(S2 > 0)/no_sims,
    # q_b_sim  = 1 - sum(ST > 0)/no_sims,

    e_1_sim = sum(S1 == 0)/n(),
    e_2_sim = sum(S2 == 0)/n(),

    ce_1 = sum(S1 ==0 & lag_1 == 0)/sum(lag_1 == 0),
    ce_2 = sum(S2 ==0 & lag_2 == 0)/sum(lag_2 == 0),

    h_1 = sum(S1 ==0 & lag_1 > 0)/sum(lag_1 > 0),
    h_2 = sum(S2 ==0 & lag_2 > 0)/sum(lag_2 > 0),

  )
  
  # 
  # %>% 
  # pivot_longer( # bring all the cols that start with q_ into two cols, one
  #   # with the values and the other with the names
  #   starts_with("e_")
  # )

sim_summ_h <- 
  sim_res_df %>% unnest(sim_list) %>% 
  group_by(sim) %>% 
  mutate(
    S1_at_tm1 = lag(S1),
    S2_at_tm1 = lag(S2),
    ST_at_tm1 = lag(ST)
  ) %>% 
  filter(ST_at_tm1 > 0) %>% 
  mutate(
    haz_1 = if_else(S1 == 0, TRUE, FALSE),
    haz_2 = if_else(S2 == 0, TRUE, FALSE),
    haz_b = if_else(ST == 0, TRUE, FALSE),
  ) %>% 
  group_by(t) %>% 
  summarise(
    hazard_1_sim = mean(haz_1),
    hazard_2_sim = mean(haz_2),
    hazard_b_sim = mean(haz_b),
  )

sim_summ_h <- 
  sim_res_df %>% unnest(sim_list) %>% 
  group_by(sim) %>% 
  mutate(
    S1_at_tm1 = lag(S1),
    S2_at_tm1 = lag(S2),
    ST_at_tm1 = lag(ST)
  ) %>% 
  filter(S2_at_tm1 == 0) %>% 
  mutate(
    ce_2 = if_else(S2 == 0, TRUE, FALSE),
    # ce_2 = if_else(S2 == 0, TRUE, FALSE),
    # ce_b = if_else(ST == 0, TRUE, FALSE),
  ) %>% 
  group_by(t) %>% 
  summarise(
    ce_2_sim = mean(ce_2),
    # ce_2_sim = mean(ce_2),
    # ce_b_sim = mean(ce_b),
  )

# plot --------------------------------------------------------------------

com_df <- left_join(full_sur_dis_df, sim_summ_e %>% select(1:3), by = c("t")) %>% 
  pivot_longer(cols = starts_with("e_")) %>% 
  select(t,name,value)


ggplot(com_df, aes(x = t, y = value, color = name)) + 
  geom_point() + 
  geom_line()
