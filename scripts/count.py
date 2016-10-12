from sys import argv, exit
from os import makedirs, path

import numpy as np
from matplotlib import pyplot as plt

import myparser


PLOTDIR = 'plots'


def makePlots(prefix, dimensions, parsedData):
    nIntervals, xyzMin, xyzMax = dimensions
    zvalues, nums, slices = parsedData

    if not path.isdir(PLOTDIR):
        makedirs(PLOTDIR)

    cMax = np.max(slices[:, :, :, 3])

    X = np.linspace(xyzMin[0], xyzMax[0], nIntervals[0] + 1)
    Y = np.linspace(xyzMin[1], xyzMax[1], nIntervals[1] + 1)
    XX, YY = np.meshgrid(X, Y)

    for i, z in enumerate(zvalues):
        countPath = path.join(PLOTDIR,
                              '{0}_count_{1}.png'.format(prefix, nums[i]))
        data = slices[i, :, :, 3]

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


if len(argv) < 2:
    exit('usage: python count.py <prefix>')

prefix = argv[1]

dimensions = myparser.readDimensions(prefix)
parsedData = myparser.parseData(prefix, dimensions)

makePlots(prefix, dimensions, parsedData)
