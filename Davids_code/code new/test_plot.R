



cascade_id_yes <- rep("yes",100)
cascade_id_no <- rep("no",10)

cascade_ids <- rep(1:100,each = 10) # create 100 cascade ids
seeded_by <- c(rep("yes",each = 500), rep("no",each = 500)) # split them between yes and no seeded, 500/500

yes_com <- sample(c("yes","no"),size = 500, replace = TRUE, prob = c(0.70, 0.3)) # first 500: sample users in them
no_com <- sample(c("yes","no"),size = 500, replace = TRUE, prob = c(0.3, 0.7)) # second 500: sample users in them

test_df <- tibble(cascade_ids, seeded_by, user_com = c(yes_com, no_com)) # put into test data frame

summ_test <- test_df %>% 
  group_by(cascade_ids, seeded_by) %>% # for each cascade id and seed_by
  summarise(
    prop_yes = sum(user_com == 'yes')/n() # count the proportion of yes and no 
    ) %>% ungroup()

# plot them, we can see the no cascades have a much lower proportion of no suporters
# 
summ_test %>% 
  ggplot(aes(x = prop_yes, fill = seeded_by))  +
  geom_density(alpha = 0.5) 
  



