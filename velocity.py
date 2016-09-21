from sys import argv, exit

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np

import myparser


def showVelocityPlot(dimensions, parsedData):
    nIntervals, xyzMin, xyzMax = dimensions
    zvalues, nums, slices = parsedData

    dx = (xyzMax[0] - xyzMin[0]) / nIntervals[0]
    X = np.linspace(xyzMin[0], xyzMax[0], nIntervals[0]) + dx/2
    Y = np.linspace(xyzMin[1], xyzMax[1], nIntervals[1]) + dx/2
    Z = np.linspace(xyzMin[2], xyzMax[2], nIntervals[2]) + dx/2
    XX, YY, ZZ = np.meshgrid(X, Y, Z, indexing='ij')

    incr = 5

    # original data is Z, X, Y
    slices = slices.transpose((2, 1, 0, 3))
    u = slices[::incr, ::incr, ::incr, 0]
    v = slices[::incr, ::incr, ::incr, 1]
    w = slices[::incr, ::incr, ::incr, 2]
    XX = XX[::incr, ::incr, ::incr]
    YY = YY[::incr, ::incr, ::incr]
    ZZ = ZZ[::incr, ::incr, ::incr]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.quiver(XX, YY, ZZ, u, v, w, length=0.3)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

if len(argv) < 2:
    exit('usage: python velocity.py <prefix>')

prefix = argv[1]

dimensions = myparser.readDimensions(prefix)
parsedData = myparser.parseData(prefix, dimensions)

showVelocityPlot(dimensions, parsedData)
