from sys import argv, exit
from os import makedirs, path

import numpy as np
from matplotlib import pyplot as plt

import myparser


PLOTDIR = 'plots'


def makePlots(prefix, dimensions, parsedData, sliceAx):
    nIntervals, xyzMin, xyzMax = dimensions
    xvalues, yvalues, zvalues, nums, slices = parsedData
    axValuesArr = [xvalues, yvalues, zvalues]

    if not path.isdir(PLOTDIR):
        makedirs(PLOTDIR)

    cMax = np.max(slices[:, :, :, 3])

    sliceAx = max(min(int(sliceAx), 2), 0)

    if sliceAx == 0:
        ax1 = 2
        ax2 = 1
    elif sliceAx == 1:
        ax1 = 2
        ax2 = 0
    else:
        ax1 = 0
        ax2 = 1

    axNames = ['X', 'Y', 'Z']

    if abs(xyzMax[ax1] - xyzMin[ax1]) < abs(xyzMax[ax2] - xyzMin[ax2]):
        ax1, ax2 = ax2, ax1

    dimidx = [2, 1, 0]
    slices = np.transpose(slices,
                          axes=[dimidx[sliceAx], dimidx[ax2], dimidx[ax1], 3])
    axvalues = axValuesArr[sliceAx]

    X = np.linspace(xyzMin[ax1], xyzMax[ax1], nIntervals[ax1] + 1)
    Y = np.linspace(xyzMin[ax2], xyzMax[ax2], nIntervals[ax2] + 1)
    XX, YY = np.meshgrid(X, Y)

    for i, z in enumerate(axvalues):
        countPath = path.join(PLOTDIR,
            '{0}_count_{1}_{2}.png'.format(prefix, axNames[sliceAx], nums[i])
        )
        data = slices[i, :, :, 3]

        heatmap = plt.pcolormesh(XX, YY, data, edgecolor='face', vmin=0,
                                 vmax=cMax)

        plt.xlabel('${0}$ [m]'.format(axNames[ax1]))
        plt.ylabel('${0}$ [m]'.format(axNames[ax2]))
        plt.xlim(xyzMin[ax1], xyzMax[ax1])
        plt.ylim(xyzMin[ax2], xyzMax[ax2])
        plt.gca().set_aspect('equal')
        plt.title('{0} = {1:.3f}'.format(axNames[sliceAx], z))
        cb = plt.colorbar(heatmap, drawedges=False)
        cb.set_label('Total particle hit count')
        plt.savefig(countPath, bbox_inches='tight')
        plt.clf()


if len(argv) < 2:
    exit('usage: python count.py <prefix> <slice axis x/y/z>')

prefix = argv[1]

sliceAxes = {'x': 0, 'y': 1, 'z': 2}
sliceAx = sliceAxes[argv[2].lower()]

dimensions = myparser.readDimensions(prefix)
parsedData = myparser.parseData(prefix, dimensions)

makePlots(prefix, dimensions, parsedData, sliceAx)

