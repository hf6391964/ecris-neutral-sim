from math import sqrt
import numpy as np
from matplotlib import pyplot as plt

data = np.log(np.loadtxt('convergence.csv', delimiter=' ', skiprows=1))
coeffs, cov = np.polyfit(data[:, 0], data[:, 1], 1, full=False, cov=True)
print(coeffs)
print(cov)
xstart = data[0, 0] - 0.5
xend = data[-1, 0] + 0.5
x2 = np.array([xstart, xend])
y2 = np.polyval(coeffs, x2)
a = coeffs[0]
b = coeffs[1]
aerr = 0.03#sqrt(cov[0, 0])
berr = 0.4#sqrt(cov[1, 1])

fit_line, = plt.plot(x2, y2, 'b-')
data_points, = plt.plot(data[:, 0], data[:, 1], 'rx')
plt.xlabel('Logarithm of particle count, $\log( N )$')
plt.ylabel('Logarithm of global error, $\log (\delta)$')
plt.xlim(xstart, xend)
plt.legend([data_points, fit_line], ['Data points', 'Line fit to points'])
fittext = "Fit parameters:\n$f(N) = ax + b$\n$a = {0:.2f}\pm{1}$\n$b={2:.1f}\pm{3}$".format(a, aerr, b, berr)
plt.gca().text(0.05, 0.07, fittext, transform=plt.gca().transAxes)
plt.savefig('convergence_plot.png', bbox_inches='tight', dpi=150)
plt.show()

