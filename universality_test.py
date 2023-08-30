import pandas as pd
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from bp_mt_class import BranchingProcessMultiType
from mtbp_analysis_class import MTBPAnalysis
from pgf_analysis_class import PGFAnalysis

# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 8, 2

# simulation parameters
n_simulations = int(5e05)
probabilities_infection = [0.06, 0.08, 0.09, 0.095]
cascades_list = []

# analysis
analysis = MTBPAnalysis()

for probability in probabilities_infection:
    file_name = 'data/cascades/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_cascades_cor.csv'
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
critical_cascades['x_coord'] = np.array(critical_cascades.index) + 1
critical_cascades['x_coord_log'] = np.log(critical_cascades['x_coord'])
critical_cascades_cleaned = critical_cascades[(critical_cascades['cascades_both_log'] != float('-inf'))
                                              & (critical_cascades['x_coord'] > 10)
                                              & (critical_cascades['x_coord'] < 100)]

# Fit a linear regression model (1st degree polynomial)
coefficients = np.polyfit(critical_cascades_cleaned['x_coord_log'], critical_cascades_cleaned['cascades_both_log'], 1)
slope = coefficients[0]
intercept = coefficients[1]

print(coefficients)

# Create a line in the log space
fit_line_log = np.polyval(coefficients, critical_cascades['x_coord_log'])
# Convert the fitted line back to the original scale
fit_line = np.exp(fit_line_log)

# plotting
for count, cascades in enumerate(cascades_list):
    plt.scatter(cascades.index, cascades['cascades_both'], s=2, label='p: '+str(probabilities_infection[count]))
# plt.scatter(critical_cascades['x_coord'], critical_cascades['cascades_both'])
plt.plot(critical_cascades['x_coord'], fit_line, color='red')
# plt.xlim([0,150])
# plt.ylim([1,45107])
# plt.gca().yaxis.get_major_locator().set_params(integer=True)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('cascade size')
plt.ylabel('probability')
plt.legend()
plt.show()

# # # pgf analysis
# cascades_large = pd.read_csv('data/cascades/p_0.08_n_sim_100000_cascades_cor.csv')
# cascades_small = pd.read_csv('data/cascades/p_0.08_n_sim_500000_cascades_cor.csv')
# cascades_tiny = pd.read_csv('data/cascades/p_0.08_n_sim_50000_cascades.csv')
#
# plt.scatter(cascades_small.index, cascades_small['cascades_both'], s=2, label=r'$n_{sim}:\; 5\cdot 10^5$')
# plt.scatter(cascades_large.index, cascades_large['cascades_both'], s=2, label=r'$n_{sim}:\; 10^5$')
# plt.scatter(cascades_tiny.index, cascades_tiny['cascades_both'], s=2, label=r'$n_{sim}:\; 5\cdot 10^4$')
# # plt.plot(critical_cascades['x_coord'], fit_line, color='red')
# # plt.xlim([0,150])
# # plt.ylim([1,45107])
# plt.gca().yaxis.get_major_locator().set_params(integer=True)
# plt.yscale('log')
# plt.xscale('log')
# plt.legend()
# plt.rc('text', usetex=True)
# plt.xlabel('cascade size')
# plt.ylabel('probability')
# # plt.xlim([0,50])
# plt.show()



