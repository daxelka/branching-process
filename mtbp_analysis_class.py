import numpy as np
import pandas as pd


class MTBPAnalysis:

    def hazard_function_in_communities(self,results_df):
        """
           Calculates the hazard function (P[N(t)=0|N(t-1)>0) for a multi-type branching process
           from simulation data stored in a Pandas DataFrame called results_df.
           The DataFrame should have columns 'generation', 'new_infections_1', 'new_infections_2',
           and 'simulation_id'.
           :return hazard function for community 1, hazard function for community 2, hazard function for both communities
           Note: faster method
        """
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
        """
          Calculates the hazard function (P[N(t)=0|N(t-1)>0) for a multi-type branching process
           from simulation data stored in a Pandas DataFrame called results_df.
           The DataFrame should have columns 'generation', 'new_infections_1', 'new_infections_2',
           and 'simulation_id'.
           :return hazard function for community 1, hazard function for community 2, hazard function for both communities
           Note: Slower method
        """
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
                          .agg(probs_1=('new_infections_1_true', 'mean'),
                               probs_2=('new_infections_2_true', 'mean'),
                               probs_both=('new_infections_both_true', 'mean'))
                          )

        # adding a row on top for the o generation
        row_generation_0 = pd.DataFrame({'probs_1': 0, 'probs_2': 0, 'probs_both': 0}, index=[0])
        hazard_function = pd.concat([row_generation_0, hazard_function]).reset_index(drop=True)
        hazard_function['gens'] = list(range(len(hazard_function['probs_both'])))

        # deleting the last row as it always corresponds to zero offsprings
        hazard_function = hazard_function.drop(hazard_function.tail(1).index)

        return hazard_function

    # def extinction_probability(self, results_df):
    #     results_df['new_infections_both'] = results_df['new_infections_1'] + results_df['new_infections_2']
    #     max_generation = max(results_df.generation)
    #     max_sim_id = max(results_df.simulation_id)
    #     data_table = []
    #     for sim in range(max_sim_id+1):
    #         data_slice = list(results_df[results_df.simulation_id == sim]['new_infections_both'])
    #         if len(data_slice) < max_generation + 1:
    #             data_table.append(data_slice + [0]*(max_generation + 1 - len(data_slice)))
    #         else:
    #             data_table.append(data_slice)
    #     extinct_probs = []
    #     for gen in range(max_generation + 1):
    #         extinct_true_false_grouped_by_gen = [1 if data_table[sim][gen] == 0 else 0 for sim in range(max_sim_id+1)]
    #         extinct_probs.append(sum(extinct_true_false_grouped_by_gen)/len(extinct_true_false_grouped_by_gen)) # take the mean
    #     return pd.DataFrame({'gens': list(range(max_generation+1)),
    #                          'probs': extinct_probs})

    def extinction_probability(self, results_df, column_name):
        results_df['new_infections_both'] = results_df['new_infections_1'] + results_df['new_infections_2']
        max_generation = max(results_df.generation)
        max_sim_id = max(results_df.simulation_id)
        data_table = []
        for sim in range(max_sim_id+1):
            data_slice = list(results_df[results_df.simulation_id == sim][column_name])
            if len(data_slice) < max_generation + 1:
                data_table.append(data_slice + [0]*(max_generation + 1 - len(data_slice)))
            else:
                data_table.append(data_slice)
        extinct_probs = []
        for gen in range(max_generation + 1):
            extinct_true_false_grouped_by_gen = [1 if data_table[sim][gen] == 0 else 0 for sim in range(max_sim_id+1)]
            extinct_probs.append(sum(extinct_true_false_grouped_by_gen)/len(extinct_true_false_grouped_by_gen)) # take the mean
        return pd.DataFrame({'gens': list(range(max_generation+1)),
                             'probs': extinct_probs})

    def extinction_probability_mt(self, results_df):
        results_df['new_infections_both'] = results_df['new_infections_1'] + results_df['new_infections_2']
        column_names = ['new_infections_both', 'new_infections_1', 'new_infections_2']
        extinction_df = pd.DataFrame()
        for column_name in column_names:
            df = self.extinction_probability(results_df, column_name)
            extinction_df[column_name] = df.probs

        extinction_df['gens'] = df.gens
        return extinction_df

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
        generations_reinfection = []
        for simulation in range(max_simulation_id + 1):
            data_current = results_df[results_df.simulation_id == simulation]
            # print(np.count_nonzero(data_current.new_infections_1 == 0))
            # print('data', np.array(data_current.new_infections_1))
            # print('extint',self.zero_sequences(np.array(data_current.new_infections_1)))

            zero_sequence, gen_reinfection = self.calc_zero_sequences(np.array(data_current.new_infections_1))
            extinction_arr = extinction_arr + zero_sequence
            generations_reinfection = generations_reinfection + gen_reinfection

        return extinction_arr, generations_reinfection


    def calc_zero_sequences(self, arr):
        first_non_zero = None
        zero_sequences = []
        gen_reinfections = []
        current_sequence = 0
        for pos, num in enumerate(arr):
            if first_non_zero is None:
                if num != 0:
                    first_non_zero = num
            else:
                if num == 0:
                    current_sequence += 1
                else:
                    if current_sequence != 0:
                        zero_sequences.append(current_sequence)
                        gen_reinfections.append(pos)
                    current_sequence = 0
        # if current_sequence != 0:
        #     zero_sequences.append(current_sequence)
        return zero_sequences, gen_reinfections

    def get_lifetime_distribution(results):
        # to count the total frequencies for each values of total_infections
        max_generation_bp = results.generation.value_counts() / results.simulation.size
        max_generation_bp = max_generation_bp.sort_values()
        lifetime_distribution = pd.DataFrame({'gens': np.array(max_generation_bp.index),
                                              'probs': list(max_generation_bp.values)})
        return lifetime_distribution

    def continuous_extinction_prob(self, sim_results, column_name='new_infections_1'):
        max_generation = max(sim_results.generation)
        probs = []  # enforce value for the generation 0
        for gen in range(1, max_generation + 1):
            # take all the simulations which survived until that gen
            reinfection_rows = sim_results[
                (sim_results['generation'] == gen) & (sim_results[column_name].shift() == 0)]
            # survived until gen-1 and have zero offspring in gen-1
            count_zeros_generation_minus_1 = sim_results[
                (sim_results['generation'] == gen - 1) & (sim_results[column_name] == 0)]

            if len(count_zeros_generation_minus_1) > 0:
                prob = len(reinfection_rows) / len(count_zeros_generation_minus_1)
            else:
                prob = np.nan

            probs.append(prob)
        return pd.DataFrame({'gens': list(range(1, max_generation + 1)),
                             'probs': probs})

    def reinfection_prob(self, sim_results, column_name='new_infections_1'):
        continuous_ext = self.continuous_extinction_prob(sim_results, column_name)
        df = continuous_ext[['gens']].copy()
        df['probs'] = 1 - continuous_ext['probs']
        return df

    # def continuous_extinction_prob(self, sim_results, column_name = 'new_infections_1'):
    #     max_generation = max(sim_results.generation)
    #     probs = []  # enforce value for the generation 0
    #     for gen in range(1, max_generation + 1):
    #         # survived until that gen and had zero in gen-1
    #         sample1 = sim_results[
    #             (sim_results['generation'] == gen) & (sim_results[column_name].shift() == 0)]
    #         # survived until gen-1 and have zero offspring in gen-1
    #         count_zeros_generation_minus_1 = sim_results[
    #             (sim_results['generation'] == gen - 1) & (sim_results[column_name] == 0)]


    # def calculate_lifetime_distribution(df, offspring_type):
    #     # Filter the DataFrame to include only the desired offspring type
    #     offspring_df = df[df['offspring_type-' + str(offspring_type)] > 0]
    #
    #     # Group the data by the 'generation' and aggregate the sum of the offspring type
    #     grouped_df = offspring_df.groupby('generation')['offspring_type-' + str(offspring_type)].sum()
    #
    #     # Calculate the cumulative sum of the offspring type
    #     cumulative_sum = np.cumsum(grouped_df)
    #
    #     # Calculate the lifetime distribution
    #     lifetime_distribution = cumulative_sum / cumulative_sum.max()
    #
    #     return lifetime_distribution
