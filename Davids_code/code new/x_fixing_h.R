d1 <- sim_res_df %>% unnest(sim_list) %>% 
  group_by(sim) %>% 
  mutate(
    S1_at_tm1 = lag(S1),
  ) %>% 
  filter(S1_at_tm1 > 0) %>% 
  mutate(
    haz_1 = if_else(S1 == 0, TRUE, FALSE),
  ) %>% 
  group_by(t) %>% 
  summarise(
    value = mean(haz_1),
  ) %>% mutate(name = 'type 1')

d2 <- sim_res_df %>% unnest(sim_list) %>% 
  group_by(sim) %>% 
  mutate(
    S2_at_tm1 = lag(S2),
  ) %>% 
  filter(S2_at_tm1 > 0) %>% 
  mutate(
    haz_2 = if_else(S2 == 0, TRUE, FALSE),
  ) %>% 
  group_by(t) %>% 
  summarise(
    value = mean(haz_2),
  ) %>% mutate(name = "type 2")

d3 <- sim_res_df %>% unnest(sim_list) %>% 
  group_by(sim) %>% 
  mutate(
    ST_at_tm1 = lag(ST)
  ) %>% 
  filter(ST_at_tm1 > 0) %>% 
  mutate(
    haz_b = if_else(ST == 0, TRUE, FALSE),
  ) %>% 
  group_by(t) %>% 
  summarise(
    value = mean(haz_b),
  ) %>% mutate(name = "both")


alt <- bind_rows(d1,d2,d3)


full_sur_dis_df %>% 
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
  # geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$h^{(i)}(t)$")) + 
  labs(color ='type') + 
  xlim(c(1,10)) + 
  geom_point(
    data = alt, aes(y = value, shape = name), 
    size = 2,
    color = "black")




# check -------------------------------------------------------------------

sim_res_df %>% unnest() %>% group_by(sim) %>% 
  mutate(
    S2_at_tm1 = lag(S2),
  ) %>% 
  filter(S2_at_tm1 == 0) %>% 
  mutate(
    haz_2 = if_else(S2 == 0, TRUE, FALSE),
  ) %>% 
  group_by(t) %>% 
  summarise(
    value = mean(haz_2),
  ) %>% mutate(name = "type 2")




# check the c probls ------------------------------------------------------


c1 <- sim_res_df %>% unnest() %>% group_by(sim) %>% 
  mutate(
    S1_at_tm1 = lag(S1),
  ) %>% 
  filter(S1_at_tm1 == 0) %>% 
  mutate(
    haz_1 = if_else(S1 == 0, TRUE, FALSE),
  ) %>% 
  group_by(t) %>% 
  summarise(
    value = mean(haz_1),
  ) %>% mutate(name = "type 1")



c2 <- sim_res_df %>% unnest() %>% group_by(sim) %>% 
  mutate(
    S2_at_tm1 = lag(S2),
  ) %>% 
  filter(S2_at_tm1 == 0) %>% 
  mutate(
    haz_2 = if_else(S1 == 0, TRUE, FALSE),
  ) %>% 
  group_by(t) %>% 
  summarise(
    value = mean(haz_2),
  ) %>% mutate(name = "type 2")

c_all <- bind_rows(c1, c2)


full_sur_dis_df %>% 
  pivot_longer(
    starts_with("c")
  ) %>% 
  mutate(
    name = case_when(
      name == "c1" ~ "type 1", 
      name == "c2" ~ "type 2"
    )
  ) %>% 
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line(size = 1) + 
  # geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab(latex2exp::TeX("$h^{(i)}(t)$")) + 
  labs(color ='type') + 
  xlim(c(1,10)) + 
  geom_point(
    data = c_all, aes(y = value, shape = name), 
    size = 2,
    color = "black")















