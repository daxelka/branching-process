import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt


random.seed(0)
n_trials = 1000
n_generations = 30

# create offspring distribution
n_offspring = np.array(range(5))
offspring_dist = np.array([0.35, 0.35, 0.1, 0.15, 0.05])

# average n_offspring
df = pd.DataFrame({'n_trial': list(range(1, n_trials)),
                    'extinct': 0,
                    'generations': 0,
                    'n_population': 0})
# create data to store simulation path
data_path = np.zeros((n_generations, n_trials))
data_path[0, :] = 1

# simulate branching process
for i in range(0, n_trials):
    population = 1
    generation = 0
    while (population > 0) and (generation < n_generations):
        # branching process
        population = np.sum(np.random.choice(n_offspring,
                                             size=population,
                                             replace=True,
                                             p=offspring_dist))
        generation += 1
        data_path[generation-1, i] = population

    df['generations'][i] = generation

    if population == 0:
        # if extinct, record extinction
        df['extinct'][i] = 1
    else:
        # if survived, record the population size
        df['n_population'][i] = population

print(df.head())

# analysis
# see which generation lives the longest
print(np.max(df['generations']))
# print(np.where(df['generations'] == np.max(df['generations'])))
# checkout all df columns sitisfying the condition
# print(df[df['generations']==30])

# estimate extinction probability
print(np.mean(df['extinct']))

# visualise
# Histogram
plt.hist(df['generations'])
plt.xlabel('generations')
plt.title('Branching process')
plt.show()
# # path
# plt.plot(data_path)
# plt.ylim(0, 100)
# plt.show()
