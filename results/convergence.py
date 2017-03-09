import numpy as np
from math import sqrt


def load_data(fname):
    return np.loadtxt(fname, delimiter=';', skiprows=2, usecols=3)


DATA_FOLDER = 'ecr3/stationary'
counts = [ 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000 ]

data_files = [ 
    (n, '{0}/test_{1}_injection_stationary.csv'.format(DATA_FOLDER, n))
    for n in counts
]

reference_data = load_data(data_files[-1][1])

with open('convergence.csv', 'w') as outfile:
    for n, fname in data_files[0:-1]:
        current_data = load_data(fname)
        error = sqrt(np.sum(np.square(current_data - reference_data)))
        print('N = {0}, error = {1}'.format(n, error))
        outfile.write('{0} {1}\n'.format(n, error))

