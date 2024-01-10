#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Calculates the probability of adoption for a node based on the number of exposures,
// adoption probability p, and social reinforcement.
// [[Rcpp::export]]
double complex_contagion(int n_exposures, double p, double social_reinforcement) {
  double prob = p * std::pow(1 - social_reinforcement, n_exposures);
  return prob;
}

// [[Rcpp::export]]
arma::sp_mat add_edges_between_communities(arma::sp_mat& g_combined, arma::uword n, double p_edge_between) {
  for (arma::uword v1 = 0; v1 < n; v1++) {
    for (arma::uword v2 = n; v2 < g_combined.n_rows; v2++) {
      if (R::runif(0, 1) < p_edge_between) {
        g_combined(v1, v2) = 1;
        g_combined(v2, v1) = 1;
      }
    }
    if (v1 % 250 == 0) Rcout << "v1: " << v1 << std::endl;
  }
  return(g_combined);
}




// new functions|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

// [[Rcpp::export]]
// Simulates the independent cascade model with a pre-quarantine step using a
// complex contagion function to determine a node's probability of adoption.
DataFrame cascade(arma::sp_mat adj_sp_mat, double p, int seed_node, double alpha, 
                  double social_reinforcement, int n_iter = 100) {
  
  // Rcout << "Starting simulation" << std::endl;  // Debugging
  
  // Initialize the sets of activated nodes, new activated nodes, and quarantined nodes.
  std::set<int> activated_nodes {seed_node};
  std::set<int> new_activated_nodes {seed_node};
  std::set<int> quarantined_nodes;
  
  // Output vectors
  std::vector<int> parent_node;
  std::vector<int> child_node;
  std::vector<int> generation;
  
  
  // Run the simulation for n_iter iterations.
  for (int iter = 0; iter < n_iter; ++iter) {
    // Rcout << "Iteration: " << iter << std::endl;  // Debugging
    
    // If there are no newly activated nodes, stop the simulation.
    if (new_activated_nodes.empty()) {
      // Rcout << "No new activated nodes. Ending simulation." << std::endl;  // Debugging
      break;
    }
    
    std::set<int> next_new_activated_nodes;
    
    // Iterate over the newly activated nodes.
    for (int node : new_activated_nodes) {
      // Rcout << "Active Node: " << node << std::endl;  // Debugging
      
      // Iterate over the node's neighbors.
      for (arma::sp_mat::const_col_iterator it = adj_sp_mat.begin_col(node); it != adj_sp_mat.end_col(node); ++it) {
        int neighbor = it.row();
        
        // If the neighbor is not already activated...
        if (activated_nodes.count(neighbor) == 0) {
          // Rcout << "Neighbor not activated: " << neighbor << std::endl;  // Debugging
          
          // Quarantine the neighbor with probability alpha.
          if (R::runif(0, 1) < alpha) {
            // Rcout << "Neighbor quarantined: " << neighbor << std::endl;  // Debugging
            quarantined_nodes.insert(neighbor);
          }
          // If the neighbor is not quarantined, attempt to activate it.
          else if (quarantined_nodes.count(neighbor) == 0) {
            // Rcout << "Attempting to activate neighbor: " << neighbor << std::endl;  // Debugging
            
            // Calculate the number of exposures the neighbor has to activated nodes.
            int n_exposures = 0;
            for (int active_node : activated_nodes) {
              if (adj_sp_mat(active_node, neighbor) != 0) {
                n_exposures++;
              }
            }
            
            // Calculate the adoption probability using the complex contagion function.
            // double adoption_probability = complex_contagion(n_exposures, p, social_reinforcement);
            double adoption_probability = p;
            // Rcout << "exp: " << n_exposures << std::endl;  // Debugging
            // Rcout << "Adoption probability: " << adoption_probability << std::endl;  // Debugging
            
            
            // Activate the neighbor with the calculated adoption probability.
            if (R::runif(0, 1) < adoption_probability) {
              // Rcout << "Neighbor activated: " << neighbor << std::endl;  // Debugging
              next_new_activated_nodes.insert(neighbor);
              
              parent_node.push_back(node);
              child_node.push_back(neighbor);
              generation.push_back(iter + 1);
            }
          }
        }
      }
    }
    
    // Update the sets of activated nodes and new activated nodes for the next iteration.
    activated_nodes.insert(next_new_activated_nodes.begin(), next_new_activated_nodes.end());
    new_activated_nodes.swap(next_new_activated_nodes);
    next_new_activated_nodes.clear();
  }
  
  // Return the DataFrame with parent, child, and generation data.
  return DataFrame::create(Named("parent") = parent_node,
                           Named("child") = child_node,
                           Named("generation") = generation);
}

// [[Rcpp::export]]
DataFrame run_simulations_cas(arma::sp_mat adj_matrix, double p, IntegerVector seed_nodes = -1, 
                              double alpha = 0, double social_reinforcement = 0, 
                              int M = 1000, int n_iter = 100) {
  
  // We will store the infection data in these vectors
  std::vector<int> simulation_num;
  std::vector<int> parents;
  std::vector<int> children;
  std::vector<int> generations;
  
  // Rcout << "start 1" << std::endl;  // Debugging
  
  for (int i = 0; i < M; ++i) {
    int seed_node;
    // Rcout << "start 2" << std::endl;  // Debugging
    // If -1 is provided, select a random node as the seed
    if (seed_nodes[0] == -1) {
      IntegerVector all_nodes = seq_len(adj_matrix.n_rows);  // Create a sequence of all node indices
      seed_node = Rcpp::sample(all_nodes, 1, false)[0] - 1; // Randomly select a node index, subtract 1 for 0-based indexing
      // Rcout << "start 3" << std::endl;  // Debugging
    } else {
      // Otherwise, randomly select a seed node from the provided vector
      seed_node = Rcpp::sample(seed_nodes, 1, false)[0] - 1; // Subtract 1 for 0-based indexing
      // Rcout << "start 4" << std::endl;  // Debugging
    }
    // Rcout << "seed node: " << seed_node << std::endl;  // Debugging
    // Run the simulation
    DataFrame result = cascade(adj_matrix, p, seed_node, alpha, social_reinforcement, n_iter);
    // Rcout << "test: " << seed_node << std::endl;  // Debugging
    // Append the results to the vectors
    int num_infections = result.nrows();
    if(num_infections > 0){
      for (int j = 0; j < num_infections; ++j) {
        simulation_num.push_back(i+1);  // Simulation number, add 1 to convert to 1-based index
        
        IntegerVector parents_col = result["parent"];
        IntegerVector children_col = result["child"];
        IntegerVector generations_col = result["generation"];
        
        parents.push_back(parents_col[j]);
        children.push_back(children_col[j]);
        generations.push_back(generations_col[j]);
      }
    }
    
  }
  
  // Create a data frame from the vectors
  return DataFrame::create(Named("simulation_num") = simulation_num,
                           Named("parent") = parents,
                           Named("child") = children,
                           Named("generation") = generations);
}
