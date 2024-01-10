rm(list = ls())
gc()

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



# test --------------------------------------------------------------------


par_set <- expand_grid(
  lout = seq(0,2,0.5),
  pin = seq(0.02, 0.099, 0.001),
  t = 0:40
) %>% 
  mutate(lin = 8)


full_sur_dis_df <- 
  par_set %>% 
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
  # group_by(pin, lout, lin) %>% #arrange(pin, lout, lin, t) %>%  
  # group_modify(~ .x %>% add_row(t = 0, hazard_1 = 0, hazard_2 = 0, hazard_b = 0, 
  #         .before = 1)) %>%
  # ungroup() %>% 
  # arrange(pin, lout, lin) %>% 
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
  ) %>% ungroup()



full_sur_dis_df %>% 
  # pivot_longer(starts_with("S_1")) %>% 
  filter(pin == 0.099) %>% 
  ggplot(aes(x = t, y = 1 - S_b, color = factor(lout))) + 
  geom_line(size = 1) + 
  geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$S^{(i)}(t)$")) + 
  # scale_y_log10() +
  # scale_x_log10() + 
  annotation_logticks(sides = "l") + 
  labs(color ='type') + 
  ylim(c(0,1))


full_sur_dis_df %>% 
  pivot_longer(starts_with("p_")) %>%
  filter(pin == 0.099) %>% 
  ggplot(aes(x = t, y = value, color = factor(lout))) + 
  geom_line(size = 1) + 
  geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$p^{(i)}(t)$")) + 
  # scale_y_log10() +
  # scale_x_log10() +
  annotation_logticks(sides = "l") + 
  labs(color ='type') + 
  ylim(c(0,0.5)) + 
  facet_wrap(~name)
  



full_sur_dis_df %>% 
  filter(lout == 0, between(pin, 0.01,0.02)) %>% 
  mutate_all( ~ replace(.,is.nan(.), 1)) %>% print(n = Inf)

sum <- full_sur_dis_df %>% 
  mutate_all( ~ replace(.,is.nan(.), 0)) %>% 
  group_by(pin, lout, lin) %>% 
  summarise(
    exp_1 = sum(t * p_1),
    exp_2 = sum(t * p_2),
    exp_b = sum(t * p_b)
  )


sum %>% 
  ggplot(aes(x = pin, y = exp_b0, color = factor(lout))) + 
  geom_line()