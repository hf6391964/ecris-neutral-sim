from matplotlib import pyplot as plt
import numpy as np

data_exc = np.loadtxt('exc.csv', delimiter=',')
data_ion = np.loadtxt('ion.csv', delimiter=',')

plt.figure(1)
pexc, = plt.semilogx(data_exc[:, 0], data_exc[:, 1])
pion, = plt.semilogx(data_ion[:, 0], data_ion[:, 1])
plt.xlabel('T (eV)')
plt.ylabel('Rate coefficient')
plt.legend([pexc, pion], ['Excitation', 'Ionization'])

plt.show()
