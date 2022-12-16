import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Davids_code/full_extin_dist.csv')

# print(df.head())
# print(list(set(df.lin)))

# my_vals = df[df['pin']==0.02]

# plt.plot(df.t[df['pin']==0.02], df.S_1[df['pin']==0.02])
plt.scatter(df.t[df['pin']==0.08], df.S_1[df['pin']==0.08])
# plt.plot(df.t[df['pin']==0.02], df.q_1[df['pin']==0.02])
# plt.plot(df.t[df['pin']==0.02], df.p_1[df['pin']==0.02])
plt.show()

