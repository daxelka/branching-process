import numpy as np
import pandas as pd


class MTBPAnalysis:

    def hazard_function_in_communities(results_df):
        max_generation = max(results_df.generation)
        hazard_prob_1 = [0]
        hazard_prob_2 = [0]
        hazard_prob_both = [0]
        for gen in range(1, max_generation):
            data_current_gen = results_df[results_df.generation == gen]
            data_previous_gen = results_df[results_df.generation == gen - 1]

            nonzeros_previous_gen_both_communities = np.count_nonzero(
                data_previous_gen.new_infections_1 + data_previous_gen.new_infections_2)
            nonzeros_current_gen_both_communities = np.count_nonzero(
                data_current_gen.new_infections_1 + data_current_gen.new_infections_2)
            nonzeros_current_gen_1 = np.count_nonzero(data_current_gen.new_infections_1)
            nonzeros_current_gen_2 = np.count_nonzero(data_current_gen.new_infections_2)

            hazard_prob_1.append(1 - nonzeros_current_gen_1 / nonzeros_previous_gen_both_communities)
            hazard_prob_2.append(1 - nonzeros_current_gen_2 / nonzeros_previous_gen_both_communities)
            hazard_prob_both.append(1 - nonzeros_current_gen_both_communities / nonzeros_previous_gen_both_communities)

        hazard_function = pd.DataFrame({'gens': list(range(max_generation)),
                                        'probs_1': hazard_prob_1,
                                        'probs_2': hazard_prob_2,
                                        'probs_both': hazard_prob_both})
        return hazard_function

    def get_max_simulation_id(self, results_df):
        return max(results_df.simulation_id)

    def cascade_distribution(self, results_df):
        max_simulation_id = self.get_max_simulation_id(results_df)
        cascades_1 = []
        cascades_2 = []
        cascades_both = []
        for simulation in range(max_simulation_id+1):
            data_slice = results_df[results_df.simulation_id == simulation]
            cascades_1.append(max(data_slice.total_infections_1))
            cascades_2.append(max(data_slice.total_infections_2))
            cascades_both.append(max(data_slice.total_infections))

        # print(cascades_both)
        # print(pd.Series(cascades_both))
        cascade_size_dist_1 = pd.Series(cascades_1).value_counts().sort_index()
        cascade_size_dist_2 = pd.Series(cascades_2).value_counts().sort_index()
        cascade_size_dist_both = pd.Series(cascades_both).value_counts().sort_index()

        return cascade_size_dist_both, cascade_size_dist_1, cascade_size_dist_2

    def duration_extinction(results_df):
        max_simulation_id = max(results_df.simulation_id)
        for simulation in range(max_simulation_id):
            data_current = results_df[results_df.simulation_id == simulation]
            print(np.count_nonzero(data_current.new_infections_1 == 0))

        return 1

    def get_lifetime_distribution(results):
        # to count the total frequencies for each values of total_infections
        max_generation_bp = results.generation.value_counts() / results.simulation.size
        max_generation_bp = max_generation_bp.sort_values()
        lifetime_distribution = pd.DataFrame({'gens': np.array(max_generation_bp.index),
                                              'probs': list(max_generation_bp.values)})
        return lifetime_distribution
