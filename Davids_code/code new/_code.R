
# code that simulate cascades for a network.

# 
p_cc_basic <- function(q1 = 0.995, alpha, n){
  return(1-q1*(1-alpha)**(n-1))
}

# same as before but with a threshold m
p_cc_threshold <- function(q1 = 0.998, alpha, n, m){
  vals <- numeric(length(n))
  vals[which(n<m)] <- 1-q1*(1-alpha[which(n<m)])**(n[which(n<m)]-1)
  vals[which(n>=m)] <- 1-q1*(1-alpha[which(n>=m)])**(m-1)
  return(vals)
}

# this is from when we were looking at the probability of adoption decreasing after
# a certain number of exposures
p_cc_threshold_drop2 <- function(q1 = 0.998, alpha, n, m = 5, beta = 0.5){
  vals <- numeric(length(n))
  vals[which(n<m)] <- 1-q1*(1-alpha[which(n<m)])**(n[which(n<m)]-1)
  vals[which(n>=m)] <- beta*(1-q1*(1-alpha[which(n>=m)])**(m-1))
  return(vals)
}

# has been altered to regular threshold version, fix seed
# you will probably need to specify a different seed, or I usually like to pick a random one each time
generate_cc_cascades <- function(follower.net = follower.net, follower.adj = follower.adj, q1 = 0.998, alpha = 0.002, total = 1000, seed_ = NULL){
  all_cascades.df <- tibble(parent = character(), child = character(), generation = numeric(), ID = numeric(), exposures = numeric())
  for (j in 1:total) {
    active <- numeric()
    inactive <- numeric()
    removed <- numeric()
    
    vertex_names <- vertex_attr(follower.net)$name
    if(is.null(seed_)){
      seed <- sample(vertex_names,1)
      active <- seed
    } else {
      active <- seed_
    }
        
    
    inactive <- vertex_names[! vertex_names %in% seed]
    
    exposures <- numeric(gorder(follower.net))
    names(exposures) <- vertex_names
    cascade.df <- tibble(parent = character(), child = character(), generation = numeric())
    generation <- 1
    while (length(active)>0) {
      new_active <- character()
      # shuffle active
      if(length(active)>1){
        active <- sample(active)
      }
      for (i in active) {
        followers <- vertex_names[which(follower.adj[,i]==1)]
        potential_adopters <- followers[followers %in% inactive]
        exposures[potential_adopters] <- exposures[potential_adopters] + 1
        if(length(potential_adopters)>0){
          # fix this, problem is with having n a vector
          adopters <- potential_adopters[runif(length(potential_adopters)) < p_cc_threshold(q1 = q1, alpha = rep(alpha, length(potential_adopters)), n = exposures[potential_adopters], m = 5)]
          if(length(adopters)>0){
            new_active <- c(new_active, adopters)
            inactive <- inactive[! inactive %in% new_active]
            cascade.df <- cascade.df %>% add_row(parent = rep(i, length(adopters)), child = adopters, generation = rep(generation, length(adopters)))
          }
        }
      }
      generation <- generation + 1
      removed <- c(removed, active)
      active <- new_active
    }
    if(nrow(cascade.df)>0){
      all_cascades.df <- all_cascades.df %>% add_row(cascade.df %>% mutate(ID = rep(j, nrow(cascade.df)), exposures = exposures[cascade.df$child]))
    }
  }
  return(all_cascades.df)
}
