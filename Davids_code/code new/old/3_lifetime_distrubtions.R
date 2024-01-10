###########################################################################
## Project: ccbp
## Script purpose: simulate lifetimes
## Date: 05-12-22
## Author: David JP O'Sullivan
###########################################################################

rm(list = ls()) # tidy work space
gc()

# libraries, source files and data ----------------------------------------

# source('./code/_project_setup.r')

library(tidyverse)
library(cowplot)
theme_set(cowplot::theme_cowplot())


M_mat <- function(lin = 9, lout = 4, 
                  pin = 0.01, pout = 0.01) {
  matrix( 
    c(lin * pin, lout * pout,
    lout * pout, lin * pin), 
    nrow = 2, ncol = 2, byrow = TRUE)
}

find_max_eigen <- function(lin = 9, lout = 4, 
                           pin = 0.01, pout = 0.01){
  mat <- M_mat(lin, lout, pin, pout)
  lam <- eigen(mat)$values
  lam_m <- lam[which.max(lam)]
  return(lam_m)
}

M_mat(pin = 0.075, pout = 0.075) %>% eigen()

g1 <- function(s1,s2, lin = 9, lout = 4, 
               pin = 0.1, pout = 0.1) exp(lin * pin * (s1 - 1)) * exp(lout * pout * (s2 - 1))
g2 <- function(s1,s2, lin = 9, lout = 4, 
               pin = 0.1, pout = 0.1) exp(lin * pin * (s2 - 1)) * exp(lout *pout* (s1 - 1))


g_ic <- function(s1,s2) s1 

g1(0,0)

g1(1,1)
g2(1,1)

library(tidyverse)

lifetime_dist <- function(pin = 0.075, pout = 0.075){
  lifetime_df <- tibble(gen = 0:50, e1 = NA, e2 = NA, et = NA)
  lifetime_df$e1[1] <- 0
  lifetime_df$e2[1] <- 1
  lifetime_df$et[1] <- 1
  for(i in 2:nrow(lifetime_df)){
    lifetime_df$e1[i] <- g1(lifetime_df$e1[i-1], lifetime_df$e2[i-1], 
                            pin = pin, pout = pout)
    lifetime_df$e2[i] <- g2(lifetime_df$e1[i-1], lifetime_df$e2[i-1],
                            pin = pin, pout = pout)
    lifetime_df$et[i] <- g2(lifetime_df$et[i-1], lifetime_df$et[i-1],
                            pin = pin, pout = pout)
    # lifetime_df$et[i] <- lifetime_df$e1[i] * lifetime_df$e2[i]
  }
  # lifetime_df_long <- 
  #   lifetime_df %>% 
  #   # mutate(et = e1 * e2) %>%
  #   pivot_longer(2:3, names_to = "type", values_to = "prob") 
  
  return(lifetime_df)  
}


lifetime_df <- lifetime_dist()

lifetime_df %>% slice(-1) %>% 
  mutate(
  S1 = cumprod(1- e1),
  S2 = cumprod(1- e2)
) %>% filter(gen < 10) %>% 
  pivot_longer(5:6) %>% 
  ggplot(aes(x = gen, y =value, color = name)) + geom_point() + 
  # scale_x_log10() + 
  scale_y_log10()



lifetime_df_long %>% filter(gen>=0) %>% ggplot(aes(x = gen, y = prob, color = type, linetype = type)) + 
  geom_line(size = 1.2) + 
  scale_y_log10() + 
  scale_x_log10() + 
  cowplot::theme_cowplot()
  

# cycle throught ----------------------------------------------------------

find_max_eigen(pin = 0.08,pout = 0.08)


par_sweep <- expand_grid(
  pin = seq(from = 0.00, to = 1, by = 0.005),
  pout = seq(from = 0.00, to = 1, by = 0.005),
  ) %>% 
  mutate(ld = list(NULL))


for(i in 1:nrow(par_sweep)){
  par_sweep$ld[[i]] <- lifetime_dist(pin = par_sweep$pin[i], pout = par_sweep$pout[i])
  if(i %% 100 == 0)print(glue::glue("Finished {i} of {nrow(par_sweep)}."))
}

par_sweep <-
  par_sweep %>% rowwise() %>% 
  mutate(l_max = find_max_eigen(pin = pin, pout = pout),
         is_sub = l_max < 1) %>% 
    unnest(ld) 


ggplot(par_sweep, aes(x = pin, y = pout, z = l_max)) + 
  geom_contour_filled() +
  geom_contour(breaks = 1, color = "red", size = 1)

# %>% filter(is_sub, pin == pout) %>% filter(max(pout) == pout)

par_sweep_small %>% count(is_sub)

par_sweep_small <- par_sweep %>% filter(
  round(pin,2) == 0.08,
  between(pout, 0.06, 0.08) | round(pout,2) == 0 | round(pout,2) == 0.01, 
  ) %>% 
  mutate(par_set = glue::glue("pin = {pin}, pout = {pout}") %>% as.character())

par_sweep_small %>% count(pin, is_sub)

par_sweep_small %>% # filter(gen>=0) %>% 
  ggplot(aes(x = gen, y = prob, color = par_set)) + 
  geom_line(size = 1.0, linetype = "solid") + 
  geom_point() + 
  scale_y_log10() +
  scale_x_log10() +
  cowplot::theme_cowplot() +
  facet_wrap(~type + pin) + 
  scale_color_viridis_d() + 
  # theme(legend.position = "") + 
  xlab("Generation") + 
  ylab(latex2exp::TeX("$P(e)$"))


par_sweep_small %>% filter(round(pout,2) == 0 | round(pout,2) == 0.01) %>% 
  ggplot(aes(x = gen, y = prob, color = par_set)) + 
  geom_line(size = 1.0, linetype = "solid") + 
  geom_point() + 
  scale_y_log10() +
  scale_x_log10() +
  cowplot::theme_cowplot() +
  facet_wrap(~type + pin) + 
  scale_color_viridis_d() + 
  # theme(legend.position = "") + 
  xlab("Generation") + 
  ylab(latex2exp::TeX("$P(e)$"))

par_sweep_small %>% filter(pin == 0.08, pout == 0.08)

# multitype version -------------------------------------------------------


g1 <- function(s1, s2, lin = 9, lout = 4, pin = 0.01) exp(lin * pin * (s1-1)) * exp(lout * pin * (s2-1))
g2 <- function(s1, s2, lin = 9, lout = 4, pin = 0.01) exp(lin * pin * (s2-1)) * exp(lout * pin * (s1-1))
gic_1 <- function(s1, s2) s1


Gn <- function(s1, s2, t, lin = 9, lout = 4, pin = 0.06){
  
  if(t > 0){
    for(i in 1:t){
      new_s1 <- g1(s1, s2, lin, lout, pin)
      new_s2 <- g2(s1, s2, lin, lout, pin)
      
      s1 <- new_s1
      s2 <- new_s2
    }
  }
  
  ans = gic_1(s1, s2)
  
  return(ans)
}

test_df <- tibble(t = 0:15) %>% rowwise() %>% 
  mutate(q = Gn(0,0,t)) %>% 
  ungroup()

test_df %>% ggplot(aes(x = t, y = q)) + geom_point()

Gh <- function(s1, s2, t, lin = 8, lout = 2, pin = 0.01){
  if(t>=1){
    gh <- Gn(g1(s1,s2), g2(s1,s2), t-1, lin, lout, pin)
    gn <- Gn(0,0,t-1, lin, lout, pin)
    
    ans <- ((gh - gn) / (1 - gn))
  } else{
    print("now we have a problem.")
  }
  return(ans)
}

pin <- 0.09
lin <- 8
lout <- 2
(lout + lin) * pin

Gh(0,0,1,lin, lout, pin)

full_sur_dis_df <- tibble(t = 1:25) %>% 
  rowwise() %>%
  mutate(
    hazard_1 = Gh(0,1,t,lin, lout, pin),
    hazard_2 = Gh(1,0,t,lin, lout, pin),
    hazard_b = Gh(0,0,t,lin, lout, pin),
    
    q_1 = Gn(0,1,t,lin, lout, pin),
    q_2 = Gn(1,0,t,lin, lout, pin),
    q_b = Gn(0,0,t,lin, lout, pin),
    
    # qp_1 = lag(q_1, default = 1) - q_1,
    # qp_2 = lag(q_2, default = 1) - q_2,
    # qp_b = lag(q_b, default = 1) - q_b
    # hazard_2 = if_else(hazard_2 > 0, hazard_2, 0)
  ) %>% 
  ungroup() %>% 
  add_row(t = 0, hazard_1 = 0, hazard_2 = 0, hazard_b = 0, 
          .before = 1) %>% 
  mutate(
    S_1 = cumprod(1 - hazard_1),
    S_2 = cumprod(1 - hazard_2),
    S_b = cumprod(1 - hazard_b),
    
    p_1 = lag(S_1, default = 1) - S_1,
    p_2 = lag(S_2, default = 1) - S_2,
    p_b = lag(S_b, default = 1) - S_b,
    
    # F_1 = cumsum(p_1),
    # F_2 = cumsum(p_2),
    # F_b = cumsum(p_b)
    # p_b
  ) 

full_sur_dis_df %>% summarise(across(starts_with("p"), .fns = ~ sum(.x)))

full_sur_dis_df %>% 
  pivot_longer(
    starts_with("q_")
  ) %>% 
  mutate(
    name = case_when(
      name == "q_1" ~ "type 1", 
      name == "q_2" ~ "type 2",
      name == "q_b" ~ "both"
    )
  ) %>% 
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line(size = 1) + 
  geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$q^{(i)}(t)$")) + 
  labs(color ='type') 


full_sur_dis_df %>% 
  pivot_longer(
    starts_with("S_")
  ) %>%
  mutate(
    name = case_when(
      name == "S_1" ~ "type 1",
      name == "S_2" ~ "type 2",
      name == "S_b" ~ "both"
    )
  ) %>%
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line(size = 1) + 
  geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$S^{(i)}(t)$")) + 
  scale_y_log10() + 
  annotation_logticks(sides = "l") + 
  labs(color ='type') 

full_sur_dis_df %>% pivot_longer(starts_with("p_")) %>% 
  mutate(
    name = case_when(
      name == "p_1" ~ "type 1",
      name == "p_2" ~ "type 2",
      name == "p_b" ~ "both"
    )
  ) %>%
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line(size = 1) + 
  geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$p^{(i)}(t)$")) + 
  scale_y_log10() + 
  annotation_logticks(sides = "l") + 
  labs(color ='type') 




# sim ---------------------------------------------------------------------

(9 + 3) * 0.076

bp_sim <- function(Z_mat, lin = 9, lout = 3, 
                   pin = 0.076, pout = 0.076){
  Z_vec <- c(1,0)
  Z_mat[1,] <- Z_vec
  for(i in 2:nrow(Z_mat)){
    Z_new <- c(0,0)
    Z_new[1] <- rpois(Z_vec[1], lin * pin) %>% sum
    Z_new[2] <- rpois(Z_vec[1], lout * pout) %>% sum
    
    Z_new[2] <- Z_new[2] + rpois(Z_vec[2], lin * pin) %>% sum
    Z_new[1] <- Z_new[1] + rpois(Z_vec[2], lout * pout) %>% sum
    
    Z_vec <- Z_new
    Z_mat[i,] <- Z_new
    if(Z_vec[1] == 0 & Z_vec[2] == 0)break
  }
  
  return(Z_mat[1:i,])
}


res_mat <- matrix(NA, ncol = 2, nrow = 50)

sim_res_df <- tibble(sim = 1:5000, sim_res = list(NULL))

for(i in 1:nrow(sim_res_df)){
  sim_res_df$sim_res[[i]] <- 
    bp_sim(res_mat) %>% 
    `colnames<-`(c('S1','S2')) %>% 
    as_tibble() %>% 
    mutate(gen = 0:(n()-1)) %>% 
    select(gen, everything())
  
  
  if(i %% 100 == 0) print(i)
  
}


survival_prob <- 
  sim_res_df %>% # filter(sim == 2) %>% 
  unnest(sim_res) %>%
  filter(gen != 0) %>%  
  mutate(
    is_S1_0 = S1 == 0,
    is_S2_0 = S2 == 0, 
    is_Sb_0
  ) %>% 
  group_by(sim) %>% 
  summarise(
    S1_diff = which(is_S1_0 == TRUE)[1] - 1,
    S2_diff = which(is_S2_0 == TRUE)[1] - 1
    ) %>% ungroup() %>% 
  filter(complete.cases(.))
  
survival_prob %>% count(S1_diff) %>% mutate(p = n/sum(n))
survival_prob %>% count(S2_diff) %>% mutate(p = n/sum(n))


# clean up work space -----------------------------------------------------

rm(list = ls()) # tidy work space
gc()