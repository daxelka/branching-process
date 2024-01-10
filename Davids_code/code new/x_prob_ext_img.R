


g_X <- function(s, p = 0.9) exp(p *( s - 1))


g_I(0)

g_I <- function(s, p = 0.1) exp(p * (s-1))
g_ic <- function(s) s * g_I(0)

G_N <- function(s,t) {
  s_old <- g_ic(s)
  s_new <- s_old

  if(t > 0){
    for(i in 1:t){
      s_new <- g_X(s_old) * g_I(0)
      s_old <- s_new
    }
  }
  return(s_old)
}

G_N(0,100000)


G_N(0,1) * (1 - G_N(0,2))

# duration of death
D <- function(t, t_len){
  
  prob <- 1
  if(t_len > 1){
    for(i in 0:(t_len-1)) prob <- G_N(0,t + i) * prob 
  }
  prob = prob * (1 - G_N(0,t + t_len))
  
  return(prob)
}


D_df <- expand_grid(t = 1:100, t_len = 1:10) %>% 
  rowwise() %>% 
  mutate(
    p_duration = D(t,t_len)
  )


E_df <- D_df %>% group_by(t) %>% 
  summarise(
    exp_duration = sum(t_len * p_duration)
  )


ggplot(D_df, aes(x = t, y = p_duration, color = t_len %>% factor)) + 
  geom_line() + 
  geom_point(size = 1) + 
  scale_color_viridis_d() + 
  scale_y_log10()


ggplot(E_df, aes(x = t, y = exp_duration)) + 
  geom_line() + 
  geom_point(size = 1) + 
  scale_color_viridis_d()






















