import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('rr.csv')
plt.figure(1)
plt.semilogx(data[:, 0], data[:, 1], 'ro', data[:, 0], data[:, 2], 'bo')

data = np.loadtxt('dr.csv')
plt.figure(2)
plt.semilogx(data[:, 0], data[:, 1], 'ro', data[:, 0], data[:, 2], 'bo')

data = np.loadtxt('total.csv')
plt.figure(3)
plt.semilogx(data[:, 0], data[:, 1], 'ro', data[:, 0], data[:, 2], 'bo')

plt.show()

