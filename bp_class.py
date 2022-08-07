import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import pandas as pd


def get_max_generation(data):
    return [len(data[str(j)]) - 1 for j in range(len(data))]


def get_average_number_offspring(results):
    offspring = {}
    for sim in range(len(results)):
        results_i = results[str(sim)]
        for gen in range(len(results_i)):
            # print(results_i[str(gen)])
            # print('gen', gen, 'results', results_i[str(gen)]['average_n_offspring'])
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

#     print('number of active nodes', len(my_dic[str(i)]['active_nodes']))
#     print('number of offspring', sum([len(k) for k in my_dic[str(i)]['offspring']]))
#     print('average number of offspring per node',
#           sum([len(k) for k in my_dic[str(i)]['offspring']])/len(my_dic[str(i)]['active_nodes']))


class BranchingProcess:

    def __init__(self, p_in=9/1000, p_out=4/1000, sizes=[500, 500], prob_infection=0.05):
        # setting branching process parameters
        self.G = nx.stochastic_block_model(sizes, [[p_in, p_out], [p_out, p_in]], seed=0)
        self.edges = self.G.edges()
        self.nodes = self.G.nodes(data=True)
        self.prob_infection = prob_infection
        self.nodes_block0 = [x for x, y in self.G.nodes(data=True) if y['block'] == 0]
        # internal parameters
        self.max_generation = 5

    def run(self):
        G_copy = nx.create_empty_copy(self.G, with_data=True)
        G_copy.add_edges_from(self.edges)
        seed = random.sample(self.nodes_block0, 1)
        active_nodes = seed
        rng = np.random.default_rng()
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
                # eligible neighbours
                neighbours_eligible = list(set(neighbours).difference(active_nodes))
                neighbours_eligible_np = np.array(neighbours_eligible)
                if neighbours:
                    infections = rng.binomial(1, self.prob_infection, (1, len(neighbours_eligible)))
                    infections_np = np.array(infections[0])
                    infected_nodes = list(neighbours_eligible_np[infections_np > 0])
                    # print('infected_nodes for active node', infected_nodes)
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
