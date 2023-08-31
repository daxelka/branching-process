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
cascades_numerical_list = []
cascades_pgf_list = []

# Cascades from numerical simulations
analysis = MTBPAnalysis()

for probability in probabilities_infection:
    file_name = 'data/cascades/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_cascades_cor.csv'
    cascades = pd.read_csv(file_name)
    cascades_numerical_list.append(cascades)

# Cascades from PGF

for probability in probabilities_infection:
    file_name = 'data/cascades/p_' + str(probability) + '_cascades_pgf.csv'
    cascades = pd.read_csv(file_name)
    cascades_pgf_list.append(cascades)


# Fit a linear regression model (1st degree polynomial)
# def fit_line():

# # cleaning data from -inf and remove the early and late point, leave the middle
# critical_cascades = cascades_numerical_list[-1]
# critical_cascades['cascades_both_log'] = np.log(critical_cascades['cascades_both'])
# critical_cascades['x_coord'] = np.array(critical_cascades.index) + 1
# critical_cascades['x_coord_log'] = np.log(critical_cascades['x_coord'])
# critical_cascades_cleaned = critical_cascades[(critical_cascades['cascades_both_log'] != float('-inf'))
#                                               & (critical_cascades['x_coord'] > 10)
#                                               & (critical_cascades['x_coord'] < 100)]
#
# coefficients = np.polyfit(critical_cascades_cleaned['x_coord_log'], critical_cascades_cleaned['cascades_both_log'], 1)
# slope = coefficients[0]
# intercept = coefficients[1]
#
# print(coefficients)

# # Create a line in the log space
# fit_line_log = np.polyval(coefficients, critical_cascades['x_coord_log'])
# # Convert the fitted line back to the original scale
# fit_line = np.exp(fit_line_log)

# plotting
for count, prob in enumerate(probabilities_infection):
    # cascades_num = cascades_numerical_list[count]
    cascades_pgf = cascades_pgf_list[count]
    # plt.scatter(cascades_num.index + 1, cascades_num['cascades_both'], s=2, label='p: '+str(round(prob/0.1, 2))+' p*')
    plt.plot(cascades_pgf.cascades, cascades_pgf.prob)

# plt.plot(critical_cascades['x_coord'], fit_line, color='red')
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-6,1])
# plt.xlim([0,20])
plt.xlabel('cascade size')
plt.ylabel('probability')
plt.legend()
plt.show()





