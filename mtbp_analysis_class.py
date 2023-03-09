import numpy as np
import pandas as pd


class MTBPAnalysis:

    def hazard_function_in_communities(self,results_df):
        max_generation = max(results_df.generation)
        hazard_prob_1 = [0]
        hazard_prob_2 = [0]
        hazard_prob_both = [0]
        for gen in range(1, max_generation):
            current_gen = results_df[results_df.generation == gen]
            previous_gen = results_df[results_df.generation == gen - 1]

            nonzeros_previous_gen_both_communities = np.count_nonzero(
                previous_gen.new_infections_1 + previous_gen.new_infections_2)
            nonzeros_current_gen_both_communities = np.count_nonzero(
                current_gen.new_infections_1 + current_gen.new_infections_2)
            nonzeros_current_gen_1 = np.count_nonzero(current_gen.new_infections_1)
            nonzeros_current_gen_2 = np.count_nonzero(current_gen.new_infections_2)

            hazard_prob_1.append(1 - nonzeros_current_gen_1 / nonzeros_previous_gen_both_communities)
            hazard_prob_2.append(1 - nonzeros_current_gen_2 / nonzeros_previous_gen_both_communities)
            hazard_prob_both.append(1 - nonzeros_current_gen_both_communities / nonzeros_previous_gen_both_communities)

        hazard_function = pd.DataFrame({'gens': list(range(max_generation)),
                                        'probs_1': hazard_prob_1,
                                        'probs_2': hazard_prob_2,
                                        'probs_both': hazard_prob_both})
        return hazard_function


    def hazard_function_in_communities_v2(self, results_df):
        results_df['new_infections_both'] = results_df['new_infections_1'] + results_df['new_infections_2']
        hazard_function = (results_df
                          .assign(new_infections_1_prev_gen=lambda x: x['new_infections_1'].shift(1),
                                  new_infections_2_prev_gen=lambda x: x['new_infections_2'].shift(1),
                                  new_infections_both_prev_gen=lambda x: x['new_infections_both'].shift(1))
                          .query('new_infections_both_prev_gen > 0')
                          .assign(new_infections_1_true=lambda x: (x['new_infections_1'] == 0).replace({True: 1, False: 0}),
                                  new_infections_2_true=lambda x: (x['new_infections_2'] == 0).replace({True: 1, False: 0}),
                                  new_infections_both_true=lambda x: (x['new_infections_both'] == 0).replace({True: 1, False: 0}))
                          .groupby('generation')
                          .agg(hazard_1=('new_infections_1_true', 'mean'),
                               hazard_2=('new_infections_2_true', 'mean'),
                               hazard_both=('new_infections_both_true', 'mean'))
                          )

        # adding a row on top for the o generation
        row_generation_0 = pd.DataFrame({'hazard_1': 0, 'hazard_2': 0, 'hazard_3': 0}, index=[0])
        hazard_function = pd.concat([row_generation_0, hazard_function]).reset_index(drop=True)

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


    def duration_extinction(self, results_df):
        max_simulation_id = max(results_df.simulation_id)
        extinction_arr = []
        for simulation in range(max_simulation_id + 1):
            data_current = results_df[results_df.simulation_id == simulation]
            # print(np.count_nonzero(data_current.new_infections_1 == 0))
            print('data', np.array(data_current.new_infections_1))
            print('extint',self.zero_sequences(np.array(data_current.new_infections_1)))

            extinction_arr = extinction_arr + self.zero_sequences(np.array(data_current.new_infections_1))

        return extinction_arr


    def zero_sequences(self, arr):
        first_non_zero = None
        zero_sequences = []
        current_sequence = 0
        for i in arr:
            if first_non_zero is None:
                if i != 0:
                    first_non_zero = i
            else:
                if i == 0:
                    current_sequence += 1
                else:
                    if current_sequence != 0:
                        zero_sequences.append(current_sequence)
                    current_sequence = 0
        if current_sequence != 0:
            zero_sequences.append(current_sequence)
        return zero_sequences

    def get_lifetime_distribution(results):
        # to count the total frequencies for each values of total_infections
        max_generation_bp = results.generation.value_counts() / results.simulation.size
        max_generation_bp = max_generation_bp.sort_values()
        lifetime_distribution = pd.DataFrame({'gens': np.array(max_generation_bp.index),
                                              'probs': list(max_generation_bp.values)})
        return lifetime_distribution