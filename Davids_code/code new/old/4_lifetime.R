
# multitype version -------------------------------------------------------
rm(list = ls())
gc()

library(tidyverse)
library(cowplot)
theme_set(cowplot::theme_cowplot())

g1 <- function(s1, s2, lin, lout, pin) {
  ans <- exp(lin * pin * (s1-1)) * exp(lout * pin * (s2-1))
  return(ans)
}

g2 <- function(s1, s2, lin, lout, pin) {
  ans <- exp(lin * pin * (s2-1)) * exp(lout * pin * (s1-1))
  return(ans)
  }
gic_1 <- function(s1, s2) s1


Gn <- function(s1, s2, t, lin, lout, pin){
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

Gh <- function(s1, s2, t, lin = 8, lout = 2, pin = 0.06){
  if(t>=1){
    s1 <- g1(s1,s2,lin, lout, pin)
    s2 <- g2(s1,s2,lin, lout, pin)
    
    gh <- Gn(s1, s2, t-1, lin, lout, pin)
    gn <- Gn(0,0,t-1, lin, lout, pin)
    
    ans <- ((gh - gn) / (1 - gn))
  } else{
    ans <- 0
  }
  return(ans)
}

pin <- 0.095
lin <- 8
lout <- 2
(lout + lin) * pin


full_sur_dis_df <- tibble(t = 1:30) %>% 
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
    
    F_1 = cumsum(p_1),
    F_2 = cumsum(p_2),
    F_b = cumsum(p_b)
    # p_b
  ) 

full_sur_dis_df %>% summarise(across(starts_with("p"), .fns = ~ sum(.x)))

p_q <- full_sur_dis_df %>% 
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

ggsave(filename = "./plots/mt_q.png", plot = p_q, width = 6, height = 5)


p_h <- full_sur_dis_df %>% 
  pivot_longer(
    starts_with("hazard_")
  ) %>% 
  mutate(
    name = case_when(
      name == "hazard_1" ~ "type 1", 
      name == "hazard_2" ~ "type 2",
      name == "hazard_b" ~ "both"
    )
  ) %>% 
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line(size = 1) + 
  geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$h^{(i)}(t)$")) + 
  labs(color ='type') 

ggsave(filename = "./plots/mt_h.png", plot = p_h, width = 6, height = 5)


p_S <- full_sur_dis_df %>% filter(t < 10) %>% 
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
  ggplot(aes(x = t, value, color = name)) + 
  geom_line(size = 1) + 
  geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$S^{(i)}(t)$")) + 
  # scale_y_log10() + 
  annotation_logticks(sides = "l") + 
  labs(color ='type') + 
  ylim(c(0,1))

ggsave(filename = "./plots/mt_S.png", plot = p_S, width = 6, height = 5)


p_pe <- full_sur_dis_df %>% 
  filter(t !=0, t < 10) %>% 
  pivot_longer(starts_with("p_")) %>% 
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



ggsave(filename = "./plots/mt_pe.png", plot = p_pe, width = 6, height = 5)



# find this ---------------------------------------------------------------


bp_sim <- function(Z_mat, lin = 8, lout = 2, 
                   pin = 0.07, pout = 0.07){
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
    if(Z_vec[1] == 0 & Z_vec[2] == 0) break
  }
  
  return(Z_mat[1:i,])
}


res_mat <- matrix(NA, ncol = 2, nrow = 100)

sim_res_df <- tibble(sim = 1:10000, sim_res = list(NULL))

for(i in 1:nrow(sim_res_df)){
  sim_res_df$sim_res[[i]] <- 
    bp_sim(res_mat) %>% 
    `colnames<-`(c('S1','S2')) %>% 
    as_tibble() %>% 
    mutate(gen = 0:(n()-1)) %>% 
    select(gen, everything())
  
  
  if(i %% 500 == 0) print(glue::glue("Currently finished {i} of {nrow(sim_res_df)}."))
  
}


survival_prob <- 
  sim_res_df %>% # filter(sim == 2) %>% 
  unnest(sim_res) %>%
  # filter(gen != 0) %>%  
  mutate(
    is_S1_0 = S1 == 0,
    is_S2_0 = S2 == 0,
    is_Sb_0 = S2 == 0 & S2 == 0,
  ) %>% 
  group_by(sim) %>% 
  summarise(
    S1_diff = list(which(is_S1_0 == TRUE)[1] - 1),
    S2_diff = list(which(is_S2_0 == TRUE)[2] - 1)
  ) %>% ungroup() 

sum_1 <- survival_prob %>%
  unnest(S1_diff) %>%
  count(t = S1_diff) %>%
  mutate(p_1_sim = n/sum(n)) %>%
  select(-n)

sum_2 <- survival_prob %>%
  unnest(S2_diff) %>%
  count(t = S2_diff) %>%
  mutate(p_2_sim = n/sum(n)) %>%
  select(-n)


# sum_1 <- tibble(S1_diff = gen_1 - 1) %>% 
#   count(t = S1_diff) %>% 
#   mutate(p_1_sim = n/sum(n)) %>% 
#   select(-n)
# 
# sum_2 <- survival_prob %>% 
#   unnest(S2_diff) %>% 
#   count(t = S2_diff) %>% 
#   mutate(p_2_sim = n/sum(n)) %>% 
#   select(-n)


all_sim_df <- left_join(sum_1,sum_2, by = 't') %>% pivot_longer(2:3)


p_ep_vs_sim <- 
  full_sur_dis_df %>% 
  filter(t !=0, t < 10) %>% 
  pivot_longer(starts_with("p_")) %>% 
  mutate(
    name = case_when(
      name == "p_1" ~ "type 1",
      name == "p_2" ~ "type 2",
      name == "p_b" ~ "both"
    )
  ) %>% filter(name != "both") %>% 
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line(size = 1) + 
  # geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$p^{(i)}(t)$")) + 
  scale_y_log10() +
  annotation_logticks(sides = "l") + 
  labs(color ='type') + 
  geom_line(data = all_sim_df %>% filter(t<10), 
            mapping = aes(x = t, y = value), 
            size = 1, linetype = "dashed") + 
  geom_point(data = all_sim_df %>% filter(t<10), 
             mapping = aes(x = t, y = value), size = 3)

p_ep_vs_sim
ggsave(filename = "./plots/mt_ep_vs_sim_l.png", p_ep_vs_sim, width = 6, height = 5)


# create larger plot ------------------------------------------------------

pin <- 0.05
lin <- 8
lout <- 2
(lout + lin) * pin

par_set <- expand_grid(
  pin = seq(0.04, 0.08, 0.02),
  t = 0:20
  ) %>% 
  mutate(lin = 8, lout = 2)


full_sur_dis_df <- par_set %>% 
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
  # group_by(pin) %>% arrange(t) %>%  
  # add_row(t = 0, hazard_1 = 0, hazard_2 = 0, hazard_b = 0, 
  #         .before = 1) %>% 
  # ungroup() %>%
  group_by(pin, lout, lin) %>% 
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



full_sur_dis_df %>% write_csv(file = "full_extin_dist_v3.csv")


