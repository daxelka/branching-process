#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Calculates the probability of adoption for a node based on the number of exposures,
// adoption probability p, and social reinforcement.
double complex_contagion(int n_exposures, double p, double social_reinforcement) {
  double prob = p * std::pow(1 - social_reinforcement, n_exposures);
  return prob;
}

// Simulates the independent cascade model with a pre-quarantine step using a
// complex contagion function to determine a node's probability of adoption.
List independent_cascade_pre_quarantine_cpp(arma::sp_mat adj_sp_mat, double p, int seed_node, double alpha, 
                                            double social_reinforcement, int n_iter = 100) {
  
  // arma::sp_mat adj_sp_mat = as<arma::sp_mat>(adj_matrix);
  // int n = adj_sp_mat.n_rows;
  
  // Initialize the sets of activated nodes, new activated nodes, and quarantined nodes.
  std::set<int> activated_nodes {seed_node};
  std::set<int> new_activated_nodes {seed_node};
  std::set<int> quarantined_nodes;
  
  // Run the simulation for n_iter iterations.
  for (int iter = 0; iter < n_iter; ++iter) {
    // If there are no newly activated nodes, stop the simulation.
    if (new_activated_nodes.empty()) {
      break;
    }
    
    std::set<int> next_new_activated_nodes;
    
    // Iterate over the newly activated nodes.
    for (int node : new_activated_nodes) {
      // Iterate over the node's neighbors.
      for (arma::sp_mat::const_col_iterator it = adj_sp_mat.begin_col(node); it != adj_sp_mat.end_col(node); ++it) {
        int neighbor = it.row();
        
        // If the neighbor is not already activated...
        if (activated_nodes.count(neighbor) == 0) {
          // Quarantine the neighbor with probability alpha.
          if (R::runif(0, 1) < alpha) {
            quarantined_nodes.insert(neighbor);
          }
          // If the neighbor is not quarantined, attempt to activate it.
          else if (quarantined_nodes.count(neighbor) == 0) {
            // Calculate the number of exposures the neighbor has to activated nodes.
            int n_exposures = 0;
            for (int active_node : activated_nodes) {
              if (adj_sp_mat(active_node, neighbor) != 0) {
                n_exposures++;
              }
            }
            // Calculate the adoption probability using the complex contagion function.
            double adoption_probability = complex_contagion(n_exposures, p, social_reinforcement);
            // Activate the neighbor with the calculated adoption probability.
            if (R::runif(0, 1) < adoption_probability) {
              next_new_activated_nodes.insert(neighbor);
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
  
  // Return the final number of active and quarantined nodes.
  return List::create(Named("num_active") = activated_nodes.size(),
                      Named("num_quarantine") = quarantined_nodes.size());
}

// [[Rcpp::export]]
DataFrame run_simulations(arma::sp_mat adj_matrix, double p, double alpha, 
                          double social_reinforcement, int M, int n_iter = 100) {
  // Convert the S4 adjacency matrix to an Armadillo sparse matrix
  // arma::sp_mat adj_matrix = Rcpp::as<arma::sp_mat>(adj_matrix_R);
  int n = adj_matrix.n_rows;
  
  IntegerVector num_active(M);
  IntegerVector num_quarantine(M);
  
  for (int i = 0; i < M; ++i) {
    // Choose a random seed node
    int seed_node = R::runif(0, n);
    
    // Run the simulation
    List result = independent_cascade_pre_quarantine_cpp(adj_matrix, p, seed_node, alpha, social_reinforcement, n_iter);
    
    // Store the number of active and quarantine nodes
    num_active[i] = result["num_active"];
    num_quarantine[i] = result["num_quarantine"];
  }
  
  return DataFrame::create(Named("num_active") = num_active,
                           Named("num_quarantine") = num_quarantine);
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


