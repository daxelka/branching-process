

s1 <- g1(s1,s2,lin, lout, pin)
s2 <- g2(s1,s2,lin, lout, pin)

g1

Ge <- function(s1, s2, t,l, lin = 8, lout = 2, pin = 0.06){
  if(t>=1){
    for(i in 1:l){
      s1_new <- g1(s1,s2,lin, lout, pin)
      s2_new <- g2(s1,s2,lin, lout, pin)
      
      s1 <- s1_new
      s2 <- s2_new
    }
        
    gh <- Gn(s1, s2, t-1, lin, lout, pin)
    gn <- Gn(0,0,t-1, lin, lout, pin)
    
    ans <- ((gh - gn) / (1 - gn))
  } else{
    ans <- 0
  }
  return(ans)
}


Ge(1,1,4,4)

le_df <- 
  expand_grid(t = 1:5, l = 1:20, lin = c(8,16)) %>% rowwise() %>% 
  mutate(
    prob_1 = 1- Ge(0,1,t,l, lin),
    prob_2 = 1- Ge(1,0,t,l, lin),
    prob_b = 1- Ge(0,0,t,l, lin)
  ) %>% 
  pivot_longer(starts_with("prob")) %>% 
  mutate(
    par_set = 
      glue::glue("t={t}, lam_I = {lin}")
  )

p_dur_e <- le_df %>%  
  ggplot(aes(x = l, y = value, color = par_set)) + 
  geom_line() + 
  geom_point() +
  # scale_x_log10() + 
  # scale_y_log10() + 
  facet_wrap(~name) + 
  xlab("duration of ext") + 
  ylab("prob")

ggsave(filename = "./plots/mt_dur_e.png", plot = p_dur_e,  width = 10, height = 4)


le_df %>%  
  ggplot(aes(x = l, y = value, color = par_set)) + 
  geom_line() + 
  geom_point() +
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap(~name)

le_df %>% filter(t == 3)

le_df %>% filter(l == 15, p)

test <- sim_res_df %>% # filter(sim == 2) %>% 
  unnest(sim_res) %>%
  # filter(gen != 0) %>%  
  mutate(
    is_S1_0 = S1 == 0,
    is_S2_0 = S2 == 0,
    is_Sb_0 = S1 == 0 & S2 == 0,
  ) 


ids <- test$sim %>% unique()
gen_1 <- numeric()
gen_2 <- numeric()
for(i in 1:length(ids)){ # i <- 1
  f_df <- test %>% filter(sim == ids[i])
  max_gen <- max(f_df$gen)
  
  for(j in 1:max_gen){ # j <- 1
    is_alive_p <- f_df$is_Sb_0[j] == FALSE
    is_dead_1 <- f_df$is_S1_0[j+1]
    is_dead_2 <- f_df$is_S2_0[j+1]
    
    if(is_alive_p == TRUE & is_dead_1) gen_1 <- c(gen_1, j + 1)
    
    if(is_alive_p == TRUE & is_dead_2) gen_2 <- c(gen_2, j + 1)
    
    
    
  }
  if(i %% 500 == 0) print(glue::glue("Finished {i} of {length(ids)}."))
}

gen_1 %>% table

  group_by(sim) %>% 
  summarise(
    S1_diff = list(which(is_S1_0 == TRUE)[1] - 1),
    S2_diff = list(which(is_S2_0 == TRUE)[2] - 1)
  ) %>% ungroup() 