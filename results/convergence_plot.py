import numpy as np
from matplotlib import pyplot as plt

data = np.log(np.loadtxt('convergence.csv', delimiter=' ', skiprows=1))
coeffs = np.polyfit(data[:, 0], data[:, 1], 1)
print(coeffs)

plt.plot(data[:, 0], data[:, 1], 'ro')
plt.show()

