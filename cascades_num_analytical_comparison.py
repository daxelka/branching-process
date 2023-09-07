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
probability = 0.08

# Numerical cascades
file_name = 'data/cascades/numerical/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_cascades.csv'
cascades_num = pd.read_csv(file_name)

# analytical cascades
file_name = 'data/cascades/from_pgfs/p_' + str(probability) + '_cascades_pgf.csv'
cascades_pgf = pd.read_csv(file_name)

# plotting
plt.scatter(cascades_num.index + 1, cascades_num['cascades_both'], s=2)
plt.plot(cascades_pgf.cascades, cascades_pgf.prob_both, label='both')
# plt.scatter(cascades_num.index + 1, cascades_num['cascades_1'], s=2)
# plt.plot(cascades_pgf.cascades, cascades_pgf.prob_1, label='type 1')
# plt.scatter(cascades_num.index + 1, cascades_num['cascades_2'], s=2)
# plt.plot(cascades_pgf.cascades, cascades_pgf.prob_2, label='type 2')

plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-6,1])
plt.xlim([1, 5e02])
plt.xlabel('cascade size')
plt.ylabel('probability')
plt.legend()
# plt.show()
plt.savefig('img/cascades_p0.08_presentation_1curveq.png', dpi=300)