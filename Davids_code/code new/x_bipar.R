# Loading necessary libraries
library(tidyverse)
library(igraph)

# This function creates a bipartite network
create_bi_network <- function(com_1_size = 10, com_2_size = 10, p = 0.4, p_rewire = 0.05, 
                              CREATE_LAYOUT = FALSE){
  # Generate random bipartite graph
  g <- sample_bipartite(n1 = com_1_size, n2 = com_2_size, p = p, directed = FALSE)
  
  if(CREATE_LAYOUT == TRUE){
    spacing <- 4
    y_pos_f <- 1
    y_pos_t <- 1
    for(i in 1:vcount(g)){
      # Check the type attribute and assign positions
      if(V(g)[i]$type == FALSE){
        V(g)[i]$x <- 1; V(g)[i]$y <- y_pos_f * spacing
        y_pos_f <- y_pos_f + 0.25
      } else {
        V(g)[i]$x <- 2; V(g)[i]$y <- y_pos_t * spacing
        y_pos_t <- y_pos_t + 0.25
      }
    }
  }
  
  # Assigning colors to nodes based on type
  V(g)[V(g)$type == FALSE]$color <- "red"
  V(g)[V(g)$type == TRUE]$color <- "blue"
  
  # Rewire the graph according to the probability p_rewire
  g_re <- rewire(g, each_edge(p = p_rewire, loops = FALSE))
  
  # Set edge colors
  E(g_re)$color = "green"
  E(g_re)[V(g_re)[type == TRUE] %--% V(g_re)[type == TRUE]]$color <- 'red'
  E(g_re)[V(g_re)[type == FALSE] %--% V(g_re)[type == FALSE]]$color <- 'red'
  
  # Return a list containing the adjacency matrices and graphs, both for the original and rewired networks
  return(list(A_re =  get.adjacency(g_re), g_re = g_re, A = get.adjacency(g), g = g))
}

# Running the create_bi_network function
bi_obj <- create_bi_network(com_1_size = 10, com_2_size = 5, CREATE_LAYOUT = TRUE)
g <- bi_obj$g
g_re <- bi_obj$g_re
A <- bi_obj$A
A_re <- bi_obj$A_re

# Plot the original and rewired bipartite networks
plot(g, vertex.label = "")
plot(g_re, vertex.label = "")

# Multiply adjacency matrices
A %*% A # strict block form
A_re %*% A_re # Almost strict block form

# Compute the Manhattan distance for the squared adjacency matrix and perform hierarchical clustering
d <- dist(A_re %*% A_re, method = "manhattan")
hc <- hclust(d)
plot(hc)

# Creating a result tibble
res <- tibble(mem = c(rep(1,10), rep(2, 5)),hc_res = cutree(hc, k=2))

res |> count(mem,hc_res)

# Simulations ----------------------------------------------------------------

# Setup for the simulations
res_df <- expand_grid(
  sim_id = 1:200, # do 200 sim 
  p = seq(0.05,0.15, 0.005) # for this list of rewiring
  )  |>  
  mutate(prop_cor = NA) # save results in here

# look at group sizes of
com_1_size <- 100
com_2_size <- 100

# Run simulations
for(i in 1:nrow(res_df)){
  p_i <- res_df$p[i]
  
  # Create the bipartite network
  adj <- create_bi_network(com_1_size, com_2_size, p = p_i)$A_re
  
  # Compute the Manhattan distance and perform hierarchical clustering
  d <- dist(adj %*% adj, method = "manhattan")
  hc <- hclust(d)
  res <- tibble(mem = c(rep(1,com_1_size), rep(2,com_2_size)),hc_res = cutree(hc, k=2))
  
  # Calculate the proportion of correct classifications by:
  prop_cor <- res |> # taking the results
    summarise(prop_cor = sum(mem == hc_res)/n()) |> # find the prop correct
    pull(prop_cor) # pull the number into a scaler object
  
  res_df$prop_cor[i] <- prop_cor # save the number
  
  # print progress
  if(i %% 100 == 0) print(glue::glue("Finished {round(i/nrow(res_df),2)*100}%"))
}

# Visualizing the distribution of prop_cor for different values of p
res_df |> 
  ggplot(aes(prop_cor, fill = as.character(p))) + 
  geom_density(alpha = 0.8) + 
  xlim(c(0.95,1))

# Grouping by p and plotting the mean proportion correct
res_df |> group_by(p) |> 
  summarise(m_p = mean(prop_cor)) |> 
  ggplot(aes(x = p, y = m_p)) + 
  geom_point() + 
  geom_line()