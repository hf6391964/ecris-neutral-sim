from sys import argv, exit
from os import makedirs, path
from math import log10

import numpy as np
from matplotlib import pyplot as plt


def normalize(data, threshold):
    data_max = np.max(data)
    mask = data > data_max * threshold
    n_plasma_cells = mask.sum()
    plasma_sum = data[mask].sum()
    average_density = plasma_sum / n_plasma_cells
    normalization = 1 / average_density
    new_data = data * normalization
    print(np.max(new_data))


def readCsvLine(f, dtype=float, sep=';'):
    return [dtype(field) for field in f.readline().split(sep)]

if len(argv) < 2:
    exit('usage: python count.py filename.csv')

prefix = argv[1]
PLOTDIR = 'plots'

with open(argv[1], 'r') as f:
    if not path.isdir(PLOTDIR):
        makedirs(PLOTDIR)

    f.readline()
    nIntervals = readCsvLine(f, int)
    f.readline()
    xyzMin = readCsvLine(f)
    f.readline()
    xyzMax = readCsvLine(f)

    nx = nIntervals[0]
    ny = nIntervals[1]
    nz = nIntervals[2]
    nxyz = nx*ny*nz
    nDecimals = int(log10(nz) + 1)
    zvalues = [xyzMin[2] + i*(xyzMax[2] - xyzMin[2]) for i in range(nz)]
    nums = [str(iz).zfill(nDecimals) for iz in range(nz)]

    slices = np.empty(nxyz)

    for i in range(nxyz):
        slices[i] = readCsvLine(f)[0]
    slices.shape = nz, ny, nx

    cMax = np.max(slices[:, :, :])

    normalize(slices, 0.1)
    normalize(slices, 0.25)
    normalize(slices, 0.5)

    i = 200

    X = np.linspace(xyzMin[0], xyzMax[0], nIntervals[0] + 1)
    Y = np.linspace(xyzMin[1], xyzMax[1], nIntervals[1] + 1)
    XX, YY = np.meshgrid(X, Y)

    data = slices[i, :, :]

    heatmap = plt.pcolormesh(XX, YY, data, edgecolor='face', vmin=0,
                             vmax=cMax)

    plt.subplot(1)
    plt.xlabel('$X$ [m]')
    plt.ylabel('$Y$ [m]')
    plt.xlim(xyzMin[0], xyzMax[0])
    plt.ylim(xyzMin[1], xyzMax[1])
    plt.gca().set_aspect('equal')
    cb = plt.colorbar(heatmap, drawedges=False)
    cb.set_label('Total particle hit count')
    plt.show()

