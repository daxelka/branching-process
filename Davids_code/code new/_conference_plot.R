ggsave(filename = "./plots/mt_pe.png", plot = p_pe, width = 6, height = 5)


ggsave(filename = "./plots/mt_pe.png", plot = p_pe, width = 5, height = 4)

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
  # geom_point(size = 1.8) + 
  xlab("generation t") + 
  ylab("probablity") + 
  scale_y_log10() +
  scale_x_continuous(breaks = c(1:10)) +
  annotation_logticks(sides = "l") + 
  labs(color ='type') 



p_cas_size <- all_df %>% filter(x != 0, x < 50) %>% 
  mutate(
    name = case_when(
      type == "1" ~ "type 1",
      type == "2" ~ "type 2",
      type == "both" ~ "both"
    )
  ) %>% 
  ggplot(aes(x = x, y = prob, color = type)) +
  geom_line(size = 1) +
  # geom_point(size = 2) +
  # geom_point() + 
  scale_y_log10() + 
  xlab("cascade size") + 
  ylab("probablity") + 
  scale_y_log10() +
  # scale_x_continuous(breaks = c(1:10)) +
  annotation_logticks(sides = "l") + 
  labs(color ='type') +
  theme(legend.position = "none")


pg_both <- cowplot::plot_grid(p_cas_size, p_pe, labels = c("(a)","(b)"), 
                   rel_widths = c(0.8,1))


ggsave(filename = './plots/netsci.png', plot = pg_both, width = 10)
