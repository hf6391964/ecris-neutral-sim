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

    i = 55
    slice_ax = 1
    horz_ax = 2
    vert_ax = 0

    X = np.linspace(xyzMin[horz_ax], xyzMax[horz_ax], nIntervals[horz_ax] + 1)
    Y = np.linspace(xyzMin[vert_ax], xyzMax[vert_ax], nIntervals[vert_ax] + 1)
    XX, YY = np.meshgrid(X, Y)

    data = np.transpose(slices, (2 - slice_ax, 2 - vert_ax, 2 - horz_ax))[i, :, :]

    heatmap = plt.pcolormesh(XX, YY, data, edgecolor='face', vmin=0,
                             vmax=cMax)

    axname = ['X', 'Y', 'Z']
    plt.xlabel('${0}$ [m]'.format(axname[horz_ax]))
    plt.ylabel('${0}$ [m]'.format(axname[vert_ax]))
    plt.xlim(xyzMin[horz_ax], xyzMax[horz_ax])
    plt.ylim(xyzMin[vert_ax], xyzMax[vert_ax])
    plt.gca().set_aspect('equal')
    plt.show()

