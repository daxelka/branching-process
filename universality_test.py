import pandas as pd
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from bp_mt_class import BranchingProcessMultiType
from mtbp_analysis_class import MTBPAnalysis

# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 8, 2

# simulation parameters
n_simulations = 100000
probabilities_infection = [0.02, 0.04, 0.06, 0.08]
cascades_list = []

# analysis
analysis = MTBPAnalysis()

for probability in probabilities_infection:
    file_name = 'data/cascades/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_cascades.csv'
    cascades = pd.read_csv(file_name)
    cascades_list.append(cascades)

# plotting
for cascades in cascades_list:
    plt.scatter(cascades.index, cascades['cascades_both'])

plt.xlim([0,170])
plt.gca().yaxis.get_major_locator().set_params(integer=True)
plt.yscale('log')
plt.show()

# cleaning data from -inf and remove the early and late point, leave the middle
critical_cascades = cascades_list[-1]
critical_cascades['cascades_both_log'] = np.log(critical_cascades['cascades_both'])
critical_cascades['x_coord'] = list(critical_cascades.index)
critical_cascades_cleaned = critical_cascades[(critical_cascades['cascades_both_log'] != float('-inf'))
                                              & (critical_cascades['x_coord'] > 10)
                                              & (critical_cascades['x_coord'] < 170)]

# Fit a linear regression model (1st degree polynomial)
x_coord = list(critical_cascades_cleaned.index)
coefficients = np.polyfit(critical_cascades_cleaned['x_coord'], critical_cascades_cleaned['cascades_both_log'], 1)
slope = coefficients[0]
intercept = coefficients[1]

print(coefficients)

# Create a line in the log space
fit_line_log = np.polyval(coefficients, critical_cascades['x_coord'])
# Convert the fitted line back to the original scale
fit_line = np.exp(fit_line_log)

# plotting
for cascades in cascades_list:
    plt.scatter(cascades.index, cascades['cascades_both'])
# plt.scatter(critical_cascades['x_coord'], critical_cascades['cascades_both'])
plt.plot(critical_cascades['x_coord'], fit_line, color='red')
# plt.xlim([0,150])
plt.ylim([1,45107])
plt.gca().yaxis.get_major_locator().set_params(integer=True)
plt.yscale('log')
plt.xscale('log')
plt.show()


# plt.plot(x_coord, fit_line, color='red', label=f'Fitted Line: y = {slope:.2f}x + {intercept:.2f}')
# plt.yscale('log')
# plt.show()
