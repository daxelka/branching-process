###########################################################################
## Project: CCCD
## Script purpose: test mode results vs simulations
## Date: 02-03-2023
## Author: David JP O'Sullivan
###########################################################################

# Include the Rcpp package
library(tidyverse)
library(Rcpp)
library(igraph)

# Load the C++ implementation
sourceCpp("./code/cascade_estimation.cpp")

set.seed(123)
network_simulation_w_com <- function(g_combined, pinf, M){
  ## simulate on the network with the community structure
  adj_comb <- as_adj(g_combined)
  sim_res_com <- run_simulations(adj_comb, pinf, 0, 0, M) |> as_tibble()
  sum_res_com <- 
    sim_res_com |> 
    count(num_active) |> 
    mutate(
      cascade_size = num_active,
      prob = n/sum(n), 
      type = "sim com")
  
  return(sum_res_com)
}

network_simulation_wo_com <- function(g_combined, pinf, M){
  ## simulate on the network with the community structure
  g_simple <- g_combined |> rewire(keeping_degseq(niter = 10000))
  adj_comb <- as_adj(g_simple)
  sim_res_com <- run_simulations(adj_comb, pinf, 0, 0, M) |> as_tibble()
  sum_res_com <- 
    sim_res_com |> 
    count(num_active) |> 
    mutate(
      cascade_size = num_active,
      prob = n/sum(n), 
      type = "sim simple")
  
  return(sum_res_com)
}


network_theory_curves <- function(pinf){
  ## simulate on the network with the community structure
  # invert pgf's and check results
  e_df3 <- 
    bind_rows(
      invert_pgf_via_ifft(t = 10, pin = pinf, pout = pinf, M = 1000) |> mutate(type = "community"),
      # invert_pgf_via_ifft(t = 10, pin = 0.05, pout = 0.0, M = 1000) |> mutate(type = "w comm 1 (off)"),
      invert_pgf_via_ifft_s(t = 10, pinf = pinf, M = 1000) |> mutate(type = "no comm"),
    ) 
  return(e_df3)
}

# test parameter values ---------------------------------------------------
1000000
M <-10^6
# pinf_vec <- seq(from = 0.05, to = 0.2, by = 0.05)
# pinf_vec <- c(0.05,0.10, 0.15, 0.25)
pinf_vec <- c(0.2)



simulation_res <- tibble()
theory_res <- tibble()

for(i in 1:length(pinf_vec)){ # i <- 1
  pinf <- pinf_vec[i]
  sim_res <- bind_rows(
    network_simulation_w_com(g_combined,pinf, M), 
    network_simulation_wo_com(g_combined,pinf, M)
  )
  
  print("Sim finished.")
  
  th_df <- network_theory_curves(pinf)
  # cut pgf results off at the same point as the simulations
  min_prob <- summarise(sim_res, min = min(prob))
  th_small_df <- 
    th_df |>  
    filter(prob >= min_prob$min)
    
  simulation_res <- bind_rows(simulation_res, sim_res |> mutate(p = pinf))
  theory_res <- bind_rows(theory_res, th_small_df |> mutate(p = pinf))
  
  print(glue::glue("Finished {i} of {length(pinf_vec)}."))
  
}


# nice plots --------------------------------------------------------------

library(latex2exp)

theory_df <- theory_res |> 
  mutate(
    type = case_when(
      type == "community" ~ "With Communities (MTBP)",
      type == "no comm" ~ "No Communities (SBP)",
    )
  ) |> 
  group_by(type) |> 
  mutate(theo = type, ccdf = 1 - cumsum(prob)) |> 
  ungroup()

simulation_df <- simulation_res |> 
  mutate(
    type = case_when(
      type == "sim com" ~ "With Communities (MTBP)",
      type == "sim simple" ~ "No Communities (SBP)",
    )
  ) |> 
  group_by(type) |> 
  mutate(ccdf = 1 - cumsum(prob)) |> 
  ungroup()

pm1 <- 
  theory_df |> 
  ggplot() + 
  geom_point(data = simulation_df,
             aes(x = cascade_size + 1, y = prob, color = type), alpha = 0.5) +
  geom_line(aes(x = cascade_size + 1, y = prob, linetype = theo), size = 1.4) + 
  scale_x_log10() +
  scale_y_log10(
    breaks = scales::breaks_log(base = 10),
    labels = scales::label_log(base = 10)
  ) + 
  xlab("Cascade size + 1") + 
  ylab("Probability") + 
  labs(linetype = "Theory", color = "Simulation") + 
  facet_wrap(~p, labeller = labeller(p = function(x) glue::glue("prob of infection = {x}"))) + 
  theme(
    legend.position = "bottom", 
    
    )
pm1

# ggsave(filename = "./plots/_nice_theory_multi_2.png", plot = pm1, bg = "white", 
#        width = 12, height = 7)


pm2 <- 
  theory_df |> 
  ggplot() + 
  geom_point(data = simulation_df |> filter(ccdf >= min(theory_df$ccdf)),
             aes(x = cascade_size + 1, y = ccdf, color = type), alpha = 0.5) +
  geom_line(aes(x = cascade_size + 1, y = ccdf, linetype = theo), size = 1.4) + 
  scale_x_log10() +
  scale_y_log10(
    breaks = scales::breaks_log(base = 10),
    labels = scales::label_log(base = 10)
  ) + 
  xlab("Cascade size + 1") + 
  ylab("CCDF") + 
  labs(linetype = "Theory", color = "Simulation") + 
  facet_wrap(~p, labeller = labeller(p = function(x) glue::glue("prob of infection = {x}"))) + 
  theme(
    legend.position = "bottom", 
    
  )
pm2

# create multiple plots of p = 0.25 ---------------------------------------
# 
op_df <-
  theory_df |>
  mutate(theo = type) |>
  filter(p == pinf_vec)

simulation_op_df <- simulation_df |> 
  filter(p == pinf_vec) |> 
  mutate(si = 
           case_when(
             type == "No Communities (SBP)" ~ "white",
             type == "With Communities (MTBP)"~ "blue")
  )

op_df |> count(type)

op0 <- op_df |>
  ggplot() +
  geom_point(data = simulation_op_df,
             aes(x = cascade_size + 1, y = prob, color = type), alpha = 0.5) +
  # geom_line(aes(x = cascade_size + 1, y = prob, linetype = theo), size = 1.4) +
  scale_x_log10() +
  scale_y_log10(
    breaks = scales::breaks_log(base = 10),
    labels = scales::label_log(base = 10)
  ) +
  xlab("Cascade size + 1") +
  ylab("Probability") +
  labs(linetype = "Theory", color = "Simulation") +
  theme(
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) + 
  scale_color_manual(values = c("white", "#00BFC4"))

op1 <- op_df |>
  ggplot() +
  geom_point(data = simulation_op_df |> filter(type == "No Communities (SBP)"),
             aes(x = cascade_size + 1, y = prob, color = type), alpha = 0.5) +
  # geom_line(aes(x = cascade_size + 1, y = prob, linetype = theo), size = 1.4) +
  scale_x_log10() +
  scale_y_log10(
    breaks = scales::breaks_log(base = 10),
    labels = scales::label_log(base = 10)
  ) +
  xlab("Cascade size + 1") +
  ylab("Probability") +
  labs(linetype = "Theory", color = "Simulation") +
  theme(
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

op2 <- op_df |>
  ggplot() +
  geom_point(data = simulation_op_df,
             aes(x = cascade_size + 1, y = prob, color = type), alpha = 0.5) +
  # geom_line(aes(x = cascade_size + 1, y = prob, linetype = theo), size = 1.4) +
  scale_x_log10() +
  scale_y_log10(
    breaks = scales::breaks_log(base = 10),
    labels = scales::label_log(base = 10)
  ) +
  xlab("Cascade size + 1") +
  ylab("Probability") +
  labs(linetype = "Theory", color = "Simulation") +
  theme(
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

op3 <- op_df |> filter(type == "No Communities (SBP)") |>
  ggplot() +
  geom_point(data = simulation_op_df,
             aes(x = cascade_size + 1, y = prob, color = type), alpha = 0.5) +
  geom_line(aes(x = cascade_size + 1, y = prob, linetype = theo), size = 1.4) +
  scale_x_log10() +
  scale_y_log10(
    breaks = scales::breaks_log(base = 10),
    labels = scales::label_log(base = 10)
  ) +
  xlab("Cascade size + 1") +
  ylab("Probability") +
  labs(linetype = "Theory", color = "Simulation") +
  theme(
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))


op4 <- op_df |>
  ggplot() +
  geom_point(data = simulation_op_df,
             aes(x = cascade_size + 1, y = prob, color = type), alpha = 0.5) +
  geom_line(aes(x = cascade_size + 1, y = prob, linetype = theo), size = 1.4) +
  scale_x_log10() +
  scale_y_log10(
    breaks = scales::breaks_log(base = 10),
    labels = scales::label_log(base = 10)
  ) +
  xlab("Cascade size + 1") +
  ylab("Probability") +
  labs(linetype = "Theory", color = "Simulation") +
  theme(
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

op0
op1
op2
op3
op4

ggsave(filename = "./plots/_nice_opening_plot_0_v1.png", plot = op0,
       bg = "white", width = 6, height = 4.5)
ggsave(filename = "./plots/_nice_opening_plot_1_v1.png", plot = op1,
       bg = "white", width = 6, height = 4.5)
ggsave(filename = "./plots/_nice_opening_plot_2_v1.png", plot = op2,
       bg = "white", width = 6, height = 4.5)
ggsave(filename = "./plots/_nice_opening_plot_3_v1.png", plot = op3,
       bg = "white", width = 6, height = 4.5)
ggsave(filename = "./plots/_nice_opening_plot_4_v1.png", plot = op4,
       bg = "white", width = 6, height = 4.5)
# 
# write_csv(x = simulation_df, file = './data/2023_08_01_simulations.csv')
# write_csv(x = theory_df, file = './data/2023_08_01_theory.csv')
