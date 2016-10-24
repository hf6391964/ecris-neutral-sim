import numpy as np
from matplotlib import pyplot as plt


data = np.loadtxt('solenoid_B_onaxis.csv', delimiter=';')

plt.plot(data[:, 0], data[:, 1])
plt.xlabel('z (m)')
plt.ylabel('B (T)')
plt.xlim([-0.16, 0.17])
plt.savefig('solenoid_B_onaxis.png')
