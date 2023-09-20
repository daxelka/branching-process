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
probabilities = list(np.linspace(0.01, 0.09, 9))

analysis = PGFAnalysis()
r1_list = []
r2_list = []
r2_max_list = []

# for p in probabilities:
#     r1, r2 = analysis.reinfection_prob(list(range(10)), lambda_in, lambda_out, p)
#     r1_list.append(r1)
#     r2_list.append(r2)
#     r2_max_list.append(r2[1])
#
# print(r2_max_list)

# for count, p in enumerate(probabilities):
#     plt.plot(list(range(10)), r2_list[count])
# plt.show()
#
# plt.scatter(probabilities, r2_max_list)
# plt.show()

# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in = 8
lambda_out_list = [0, 1, 2, 3, 4, 5, 6, 7, 8]
probability = 0.06

r1_list = []
r2_list = []
r2_max_list = []

for l in lambda_out_list:
    r1, r2 = analysis.reinfection_prob(list(range(10)), lambda_in, l, probability)
    r1_list.append(r1)
    r2_list.append(r2)
    r2_max_list.append(r2[1])

print(r2_max_list)

plt.scatter(lambda_out_list, r2_max_list)
plt.show()

for count, l in enumerate(lambda_out_list):
    plt.plot(list(range(10)), r2_list[count])
plt.show()
