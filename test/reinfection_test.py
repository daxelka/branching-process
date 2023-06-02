import numpy as np
import pandas as pd

sim_results = pd.DataFrame({
    'simulation_id': [0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3],
    'generation': [0, 1, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 3, 4],
    'new_infections_1': [1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 2, 0, 1, 0],
})


def calculate_reinfection_prob(sim_results, column_name='new_infections_1'):
    max_generation = max(sim_results.generation)
    reinfection_probs = [] # enforce value for the generation 0
    for gen in range(1, max_generation + 1):
        # take all the simulations which survived until that gen
        reinfection_rows = sim_results[
            (sim_results['generation'] == gen) & (sim_results[column_name].shift() == 0)]
        # survived until gen-1 and have zero offspring in gen-1
        count_zeros_generation_minus_1 = sim_results[
            (sim_results['generation'] == gen-1) & (sim_results[column_name] == 0)]

        if len(count_zeros_generation_minus_1) > 0:
            prob = len(reinfection_rows)/len(count_zeros_generation_minus_1)
        else:
            prob = np.nan

        reinfection_probs.append(prob)
    return pd.DataFrame({'gens': list(range(1, max_generation + 1)),
                         'probs': reinfection_probs})


reinfections = calculate_reinfection_prob(sim_results, 'new_infections_1')
print('Returned results: ', reinfections)
print('True result: ', pd.DataFrame({'gens': [1,2,3,4], 'probs': [np.nan, 0.5, 0.5, 0]}))

