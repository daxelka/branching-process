
node <- ends(g_temp,sample(E(g_temp),5000, replace = T))[,1]
table(degree(g_temp, node)-1) %>% prop.table %>% round(.,2)

mean(degree(g_temp, node)-1)

degree(g_temp) %>% table %>% prop.table 

degree(g_temp) %>% mean