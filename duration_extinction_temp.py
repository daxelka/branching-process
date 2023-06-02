import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mtbp_analysis_class import MTBPAnalysis
from pgf_analysis_class import PGFAnalysis
from plotting_class import Plotting
import time

lin, lout = 8, 2
pin = 0.09

# from simulation
sim_results_full = pd.read_csv('sim_data/mtbp_pin_0.09.csv')
sim_results = sim_results_full.iloc[:, :5]

analysis = MTBPAnalysis()


def calculate_reinfection_prob(self, results_df, column_name):
    results_df['new_infections_both'] = results_df['new_infections_1'] + results_df['new_infections_2']
    max_generation = max(results_df.generation)
    max_sim_id = max(results_df.simulation_id)
    data_table = []
    for sim in range(max_sim_id + 1):
        data_slice = list(results_df[results_df.simulation_id == sim][column_name])
        if len(data_slice) < max_generation + 1:
            data_table.append(data_slice + [0] * (max_generation + 1 - len(data_slice)))
        else:
            data_table.append(data_slice)
    extinct_probs = []
    return data_table
    # for gen in range(max_generation + 1):
    #     extinct_true_false_grouped_by_gen = [1 if data_table[sim][gen] == 0 else 0 for sim in range(max_sim_id + 1)]
    #     extinct_probs.append(
    #         sum(extinct_true_false_grouped_by_gen) / len(extinct_true_false_grouped_by_gen))  # take the mean
    # return pd.DataFrame({'gens': list(range(max_generation + 1)),
    #                      'probs': extinct_probs})

def calculate_frequency(arr):
    series = pd.Series(arr)
    frequency = series.value_counts(normalize=True).to_dict()
    return frequency

# reinfection_prob = []
# for n in range(max(sim_results.generation)):
#     reinfection_prob.append(reinfection_probability(sim_results, n))

# plt.plot(list(range(max(sim_results.generation))), reinfection_prob)
# plt.show()

duration_extinction_arr, gen_reinfections_arr = analysis.duration_extinction(sim_results)

duration_extinction_freq = calculate_frequency(duration_extinction_arr)
gen_reinfections_freq = calculate_frequency(gen_reinfections_arr)

def plot_frequency(frequency):
    x = list(frequency.keys())
    y = list(frequency.values())

    plt.scatter(x, y)
    plt.xlabel('Numbers')
    plt.ylabel('Frequency')
    plt.title('Duration of extinction')
    plt.show()

plot_frequency(duration_extinction_freq)
plot_frequency(gen_reinfections_freq)

#
# def test_zero_sequences():
#     arr = [0, 0, 1, 0, 0, 0, 2, 3, 0, 0, 0, 0, 4, 5, 0]
#     expected_result = [3, 4]
#     result = analysis.calc_zero_sequences(arr)
#     print(result)
#     print(expected_result)
#
# test_zero_sequences()