g1 <- function(s1, s2,c1,c2, lin, lout, pin) {
  ans <- exp(lin * pin * (s1*c1-1)) * exp(lout * pin * (s2*c2-1))
  return(ans)
}

g2 <- function(s1, s2,c1,c2, lin, lout, pin) {
  ans <- exp(lin * pin * (s2 * c2 -1)) * exp(lout * pin * (s1 * c1-1))
  return(ans)
}
gic <- function(s1, s2,c1,c2) s1 * c1

gc <- function(c) c

# function to itterate each type
itt_gen <- function(s1, s2,c1,c2, lin = 8, lout = 2, pin = 0.03, n_max = 10) {
  # s1 and s2: dummy var of for each type
  # lam par for poisson dist
  # n_max: number of generation 
  for (n in n_max:1) { # iterate the generating
    new_s1 <- g1(s1, s2,c1,c2, lin, lout, pin)
    new_s2 <- g2(s1, s2,c1,c2, lin, lout, pin)
    
    new_c1 <- gc(c1)
    new_c2 <- gc(c2)
    
    s1 <- new_s1
    s2 <- new_s2
    c1 <- new_c1
    c2 <- new_c2
  }
  
  # apply initial conditions
  ans <- gic(s1, s2,c1,c2)
  
  return(ans)
}

itt_gen(1,1,1,1,10,2,0.01, 100)


M <- 1000
x <- exp(-2*pi*1i*(0:(M-1))/M)
get_inverse <- function(c1,c2,lin = 8, lout = 2, pin = 0.095){
  # coefficints are mixed up
  pdf <- Re(fft(itt_gen(s1 = 1,s2 = 1, c1, c2,lin, lout, pin), inverse = TRUE))
  pdf[pdf <= 0] <- 0 # replace any negatative coeff
  pdf <- pdf / length(pdf) # require to sub to one
  dist_poi <- tibble(x = 0:(length(x)-1), prob = pdf)
  return(dist_poi)
}


all_df <- get_inverse(x,x) %>% mutate(type = "both") %>% 
  bind_rows(., get_inverse(x,1) %>% mutate(type = "1")) %>% 
  bind_rows(., get_inverse(1,x) %>% mutate(type = "2")) %>% 
  filter()


all_df %>% 
  ggplot(aes(x = x, y = prob, color = type)) +
  geom_line() +
  # geom_point() + 
  scale_y_log10() + 
  xlab("cascade size") + 
  ylab("probobablity")


# sim for comp ------------------------------------------------------------

bp_sim <- function(Z_mat, lin = 8, lout = 2, 
                   pin = 0.03, pout = 0.03){
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

sim_res_df <- tibble(sim = 1:20000, C1 = NA, C2 = NA, Cb = NA)

for(i in 1:nrow(sim_res_df)){
    
    res <- 
      bp_sim(res_mat) %>% 
      `colnames<-`(c('S1','S2')) %>% 
      as_tibble() %>% 
      mutate(gen = 0:(n()-1)) %>% 
      select(gen, everything())
    
    sim_res_df$C1[i] <- res$S1 %>% sum
    sim_res_df$C2[i] <- res$S2 %>% sum
    sim_res_df$Cb[i] <- sim_res_df$C1[i] + sim_res_df$C2[i]
    
  
  if(i %% 500 == 0) print(glue::glue("Currently finished {i} of {nrow(sim_res_df)}."))
  
}




pdf_c1 <- sim_res_df %>% count(C1) %>% mutate(p = n/sum(n), sim = "sim_1") %>% rename(x = C1)
pdf_c2 <- sim_res_df %>% count(C2) %>% mutate(p = n/sum(n), sim = "sim_2") %>% rename(x = C2)
pdf_cb <- sim_res_df %>% count(Cb) %>% mutate(p = n/sum(n), sim = "sim_b") %>% rename(x = Cb)

sim_pdf <- bind_rows(pdf_c1, pdf_c2, pdf_cb) %>% filter()

p_sim_cascade_size <- all_df %>% filter(x < 13) %>% 
  ggplot(aes(x = x, y = prob, color = type)) +
  geom_line(size = 1) +
  # geom_point() + 
  scale_y_log10() + 
  xlab("cascade size") + 
  ylab("probobablity") + 
  geom_point(data = sim_pdf, aes(y = p, shape = sim), size = 2.2, color = "black") + 
  geom_line(data = sim_pdf, aes(y = p, shape = sim, group = sim), color = "grey", linetype = "dashed")

ggsave(filename = "./plots/mt_pe.png", plot = p_pe, width = 5, height = 4)

ggsave(filename = "./plots/mt_cas_size_vs_sim.png", plot = p_sim_cascade_size)



# weakly connected comms --------------------------------------------------

all_df <- get_inverse(x, x, lin = 5, lout = 5, pin = 0.02) %>% mutate(type = "equal") %>% 
  bind_rows(., get_inverse(x, x,lin = 6, lout = 4, pin = 0.02) %>% mutate(type = "all most")) %>% 
  bind_rows(., get_inverse(x, x,lin = 9, lout = 1, pin = 0.02) %>% mutate(type = "ext")) %>% 
  filter()


all_df %>% filter(x<50) %>% 
  ggplot(aes(x = x, y = prob, color = type)) +
  geom_line(size = 1.2, aes(linetype = type)) +
  # geom_point() + 
  scale_y_log10() + 
  xlab("cascade size") + 
  ylab("probability")
