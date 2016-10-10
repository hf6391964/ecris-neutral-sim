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


def parseData(prefix, dimensions):
    nIntervals, xyzMin, xyzMax = dimensions
    nx = nIntervals[0]
    ny = nIntervals[1]
    nz = nIntervals[2]
    nxyz = nx*ny*nz
    nDecimals = int(log10(nz) + 1)
    zvalues = [xyzMin[2] + i*(xyzMax[2] - xyzMin[2]) for i in range(nz)]
    nums = [str(iz).zfill(nDecimals) for iz in range(nz)]
    slices = np.empty((nxyz, 4))

    csvFilename = prefix + '_stationary.csv'
    if not path.isfile(csvFilename):
        exit('File {0} not found, exiting'.format(csvFilename))

    with open(csvFilename, 'r') as f:
        f.readline()

        for i in range(nxyz):
            slices[i] = readCsvLine(f)

    slices.shape = nz, ny, nx, 4
    return zvalues, nums, slices
