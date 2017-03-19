from os import path
from math import log10
import numpy as np

CSV_SEP = ';'


def readCsvLine(f, dtype=float, sep=CSV_SEP):
    return [dtype(field) for field in f.readline().split(CSV_SEP)]


def readDimensions(prefix):
    filename = prefix + '_dimensions.csv'
    if not path.isfile(filename):
        exit('No dimension file found for given prefix, exiting')

    with open(filename, 'r') as f:
        f.readline()
        nIntervals = readCsvLine(f, int)
        f.readline()
        xyzMin = readCsvLine(f)
        f.readline()
        xyzMax = readCsvLine(f)
        return nIntervals, xyzMin, xyzMax


def parseData(csvFilename, dimensions):
    nIntervals, xyzMin, xyzMax = dimensions
    nx = nIntervals[0]
    ny = nIntervals[1]
    nz = nIntervals[2]
    nxyz = nx*ny*nz
    nDecimals = int(log10(nz) + 1)
    xvalues = [xyzMin[0] + i/nx*(xyzMax[0] - xyzMin[0]) for i in range(nx)]
    yvalues = [xyzMin[1] + i/ny*(xyzMax[1] - xyzMin[1]) for i in range(ny)]
    zvalues = [xyzMin[2] + i/nz*(xyzMax[2] - xyzMin[2]) for i in range(nz)]
    nums = [str(iz).zfill(nDecimals) for iz in range(nz)]
    slices = np.empty((nxyz, 4))

    if not path.isfile(csvFilename):
        exit('File {0} not found, exiting'.format(csvFilename))

    with open(csvFilename, 'r') as f:
        f.readline()
        t = f.readline()
        t = float(t[t.find('=') + 1 : -1])

        for i in range(nxyz):
            slices[i] = readCsvLine(f)

    slices.shape = nz, ny, nx, 4
    return xvalues, yvalues, zvalues, nums, t, slices
