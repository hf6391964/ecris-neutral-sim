import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('data.csv', delimiter=',')

plt.semilogx(data[:, 0], data[:, 1])
plt.show()
