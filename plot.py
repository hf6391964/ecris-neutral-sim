from sys import argv, exit
from os import path, makedirs

from math import log10
import numpy as np
from matplotlib import pyplot as plt


CSV_SEP = ';'
PLOTDIR = 'plots'


def readCsvLine(f, dtype=float, sep=CSV_SEP):
    return [dtype(field) for field in f.readline().split(CSV_SEP)]


def readDimensions(filename):
    with open(filename, 'r') as f:
        f.readline()
        nIntervals = readCsvLine(f, int)
        f.readline()
        xyzMin = readCsvLine(f)
        f.readline()
        xyzMax = readCsvLine(f)
        return nIntervals, xyzMin, xyzMax


def parseData(prefix, dimensions):
    nIntervals, xyzMin, xyzMax = dimensions
    nx = nIntervals[0]
    ny = nIntervals[1]
    nz = nIntervals[2]
    nxy = nx*ny
    nDecimals = int(log10(nz) + 1)
    zvalues = np.empty((nz))
    slices = np.empty((nz, nxy, 4))
    nums = []
    for iz in range(nz):
        num = str(iz).zfill(nDecimals)
        filename = '{0}_{1}'.format(prefix, num)
        nums.append(num)
        csvFilename = filename + '.csv'
        if not path.isfile(csvFilename):
            exit('File {0} not found, exiting'.format(filename))

        with open(csvFilename, 'r') as f:
            z = float(f.readline().split('=')[1])
            zvalues[iz] = z
            f.readline()

            for i in range(nxy):
                slices[iz, i] = readCsvLine(f)

    slices.shape = nz, ny, nx, 4
    return zvalues, nums, slices


def makePlots(prefix, dimensions, parsedData):
    nIntervals, xyzMin, xyzMax = dimensions
    zvalues, nums, slices = parsedData

    if not path.isdir(PLOTDIR):
        makedirs(PLOTDIR)

    cMax = np.max(slices[:, :, :, 3])

    for i, z in enumerate(zvalues):
        countPath = path.join(PLOTDIR,
                              '{0}_count_{1}.png'.format(prefix, nums[i]))
        velPath = path.join(PLOTDIR,
                            '{0}_vel_{1}.png'.format(prefix, nums[i]))

        X = np.linspace(xyzMin[0], xyzMax[0], nIntervals[0] + 1)
        Y = np.linspace(xyzMin[1], xyzMax[1], nIntervals[1] + 1)
        XX, YY = np.meshgrid(X, Y)
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
        # plt.show()
        plt.clf()


if len(argv) < 2:
    exit('usage: python plot.py <prefix>')

prefix = argv[1]
dimensionFileName = prefix + '_dimensions.csv'
if not path.isfile(dimensionFileName):
    exit('No dimension file found for given prefix, exiting')

dimensions = readDimensions(dimensionFileName)
parsedData = parseData(prefix, dimensions)

makePlots(prefix, dimensions, parsedData)
