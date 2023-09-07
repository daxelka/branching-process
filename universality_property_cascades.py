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
for probability in probabilities_infection:
    file_name = 'data/cascades/numerical/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_cascades.csv'
    cascades = pd.read_csv(file_name)
    cascades_numerical_list.append(cascades)

# Cascades from PGF
for probability in probabilities_infection:
    file_name = 'data/cascades/from_pgfs/p_' + str(probability) + '_cascades_pgf.csv'
    cascades = pd.read_csv(file_name)
    cascades_pgf_list.append(cascades)


def fit_line(x_array, y_array):
    x_array_log = np.log(x_array)
    y_array_log = np.log(y_array)
    sample = y_array_log != float('-inf')
    y_array_log_cleaned = y_array_log[sample]
    x_array_log_cleaned = x_array_log[sample]

    # Fit a linear regression model (1st degree polynomial)
    coefficients = np.polyfit(x_array_log_cleaned, y_array_log_cleaned, 1)

    # Create a line in the log space
    fitted_line_log = np.polyval(coefficients, x_array_log_cleaned)

    # Convert the fitted line back to the original scale
    fitted_line = np.exp(fitted_line_log)

    return coefficients, fitted_line


# Fitting line to near critical cascade
p_critical = 0.099
critical_cascade = pd.read_csv('data/cascades/from_pgfs/p_' + str(p_critical) + '_cascades_pgf.csv')
x_array = critical_cascade.cascades
y_array = critical_cascade.prob
coefficients, fitted_line = fit_line(x_array, y_array)
print(coefficients)


# plotting
for count, prob in enumerate(probabilities_infection):
    cascades_num = cascades_numerical_list[count]
    cascades_pgf = cascades_pgf_list[count]
    plt.scatter(cascades_num.index + 1, cascades_num['cascades_both'], s=2, label='p: '+str(round(prob/0.1, 2))+' p*')
    plt.plot(cascades_pgf.cascades, cascades_pgf.prob)

plt.plot(x_array, fitted_line, color='black', linewidth=3)
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-6,1])
plt.xlabel('cascade size')
plt.ylabel('probability')
plt.legend()
plt.show()





