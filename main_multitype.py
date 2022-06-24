import pandas as pd
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from multitypy_branching_process import BranchingProcessMultiType

# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 6, 6
probability_in, probability_out = 0.05, 0.05

# simulation parameters
n_simulations = 10

# initiate branching process
bp = BranchingProcessMultiType(seed_1, seed_2,
                               lambda_in, lambda_out,
                               probability_in, probability_out)


# run branching process n_simulation types
t1 = time.time()
sim_results = bp.branching()
print('elapsed time:', time.time() - t1)
print(sim_results.iloc[:, 1:5])
