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
