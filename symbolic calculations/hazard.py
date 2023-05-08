import sympy as sp
import numpy as np
import math


# B = sp.IndexedBase('B')
# B_values = {(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4}
#
# # Get the value of B[1,0]
# result = B[1, 0].subs(B_values).evalf()
#
# print(result)

# Define the symbols for the two random variables
x, y, i, j = sp.symbols('x y i j')
# x, y = sp.symbols('x y')
#
# # Define the probability mass function (PMF) of the bivariate distribution
P_values = {(0, 0): 1/4, (0, 1): 1/8, (1, 0): 1/8, (1, 1): 1/2}

P = sp.IndexedBase('P')
#
# # Set the values of B using subs method
P = P.subs(P_values)
# print(P[0,0].subs(P_values).evalf())


# # Define the bivariate probability generating function using the PMF
gxy = sp.Sum(sp.Sum(x**i * y**j * P[(i, j)], (i, 0, 1)), (j, 0, 1))

# Print the bivariate probability generating function
print("The bivariate probability generating function is:")
sp.pprint(gxy)
sp.pprint(gxy.subs({x: 1, y: 1}).expand())

