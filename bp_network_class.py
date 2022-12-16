import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import pandas as pd


def get_max_generation(data):
    return [len(data[str(j)]) - 1 for j in range(len(data))]


def get_lifetime_distribution(data, extend_to=False):
    generations = np.array(get_max_generation(data))
    max_generation = max(generations)
    n_sim = len(generations)
    lifetime_probs = [len(generations[generations >= n]) / n_sim for n in
                      range(1, max(generations) + 1)]
    if extend_to and extend_to > max_generation:
        lifetime_distribution = pd.DataFrame({'gens': list(range(1, extend_to + 1)),
                                              'probs': lifetime_probs + [0]*(extend_to-max_generation)})
    else:
        lifetime_distribution = pd.DataFrame({'gens': list(range(1, max(generations) + 1)),
                                              'probs': lifetime_probs})
    # max_generation = get_max_generation(data)
    # # lifetime distribution
    # vals, frequencies = np.unique(max_generation, return_counts=True)
    # lifetime_vals = vals + 1  # counting first generation
    # lifetime_frequencies = frequencies / len(max_generation)  # normalising
    # lifetime_distribution = pd.DataFrame({'gens': list(lifetime_vals),
    #                                       'probs': list(lifetime_frequencies)})
    return lifetime_distribution


def get_average_number_offspring(results):
    offspring = {}
    for sim in range(len(results)):
        results_i = results[str(sim)]
        for gen in range(len(results_i)):
            if str(gen) in offspring.keys():
                offspring[str(gen)].append(results_i[str(gen)]['average_n_offspring'])
            else:
                offspring[str(gen)] = []
                offspring[str(gen)].append(results_i[str(gen)]['average_n_offspring'])

    # calculate an average number of offspring
    average_n_offspring = dict.fromkeys(offspring.keys())
    for gen_str in offspring.keys():
        average_n_offspring[gen_str] = sum(offspring[gen_str]) / len(offspring[gen_str])
    return average_n_offspring


def get_total_infections(results, community_size, block='total'):
    total_infections_multi_runs = {}
    for sim in range(len(results)):
        results_i = results[str(sim)]
        total_infections = []
        # infections_block0 = []
        # infections_block1 = []

        for gen in range(len(results_i)):
            # take list of active nodes in the generation
            active_total = np.array(results_i[str(gen)]['active_nodes'])

            if block == 'total':
                total_infections.append(len(active_total))
            elif block == 'block0':
                active_bloc0 = active_total[active_total < community_size]
                total_infections.append(len(active_bloc0))
            elif block == 'block1':
                active_bloc1 = active_total[active_total >= community_size]
                total_infections.append(len(active_bloc1))
            # If an exact match is not confirmed
            else:
                raise ValueError('unknown value for parameter block')

        total_infections_multi_runs[str(sim)] = np.cumsum(total_infections)

    return total_infections_multi_runs


class BranchingProcessNetwork:

    def __init__(self, p_in=9/1000, p_out=4/1000, sizes=[500, 500], prob_infection=0.05):
        # setting branching process parameters
        self.G = nx.stochastic_block_model(sizes, [[p_in, p_out], [p_out, p_in]], seed=0)
        self.edges = self.G.edges()
        self.nodes = self.G.nodes(data=True)
        self.prob_infection = prob_infection
        self.nodes_block0 = [x for x, y in self.G.nodes(data=True) if y['block'] == 0]
        # internal parameters
        self.max_generation = 100
        self.rng = np.random.default_rng()

    def run(self):
        # make a copy of the network for this simulation
        G_copy = nx.create_empty_copy(self.G, with_data=True)
        G_copy.add_edges_from(self.edges)
        # pick a seed node from block0
        seed = random.sample(self.nodes_block0, 1)
        active_nodes = seed
        results_dic = {'0': {}}
        generation = 0
        while active_nodes and generation <= self.max_generation:
            results = {}
            infected_nodes_total = []
            # print('active nodes', active_nodes)
            for active_node in active_nodes:
                # print('active node', active_node)
                # take IDs of all nodes that are connected by a edge which received a possitive outcome in Bernoulli trials
                neighbours = list(G_copy.neighbors(active_node))
                # exclude active nodes from list of neighbours
                neighbours_eligible = list(set(neighbours).difference(active_nodes))
                neighbours_eligible_np = np.array(neighbours_eligible)
                if neighbours_eligible:
                    infections = self.rng.binomial(1, self.prob_infection, (1, len(neighbours_eligible)))
                    infections_np = np.array(infections[0])
                    infected_nodes = list(neighbours_eligible_np[infections_np > 0])
                    # record all infections for this generation (with possibility of a node being infected twice)
                    infected_nodes_total = infected_nodes_total + infected_nodes

                # remove current active node
                G_copy.remove_node(active_node)

            # recording results
            results['active_nodes'] = active_nodes
            results['offspring'] = list(set(infected_nodes_total))  # remove dublicates
            results['n_active'] = len(active_nodes)
            results['average_n_offspring'] = len(set(infected_nodes_total)) / len(active_nodes)
            results_dic[str(generation)] = results
            # print('results_dic', results_dic[str(generation)])
            active_nodes = list(set(infected_nodes_total))
            # print('active nodes for next generation', active_nodes)
            generation += 1
        return results_dic

    def simulations(self, n_simulations=10):
        simulations = {}
        for i in range(n_simulations):
            results = self.run()
            simulations[str(i)] = results
        return simulations
