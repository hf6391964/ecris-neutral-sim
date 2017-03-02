from sys import argv, exit
from os import makedirs, path
from math import log10

import numpy as np
from matplotlib import pyplot as plt


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

    X = np.linspace(xyzMin[0], xyzMax[0], nIntervals[0] + 1)
    Y = np.linspace(xyzMin[1], xyzMax[1], nIntervals[1] + 1)
    XX, YY = np.meshgrid(X, Y)

    for i, z in enumerate(zvalues):
        countPath = path.join(PLOTDIR,
                              'electron_count_{0}.png'.format(nums[i]))
        data = slices[i, :, :]

        heatmap = plt.pcolormesh(XX, YY, data, edgecolor='face', vmin=0,
                                 vmax=cMax)

        plt.xlabel('$X$ [m]')
        plt.ylabel('$Y$ [m]')
        plt.xlim(xyzMin[0], xyzMax[0])
        plt.ylim(xyzMin[1], xyzMax[1])
        plt.gca().set_aspect('equal')
        plt.title('Z = ' + str(z))
        cb = plt.colorbar(heatmap, drawedges=False)
        cb.set_label('Total particle hit count')
        plt.savefig(countPath, bbox_inches='tight')
        plt.clf()

