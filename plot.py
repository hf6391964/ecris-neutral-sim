from sys import argv, exit
from os import path

from math import log10
import numpy as np
from matplotlib import pyplot as plt


CSV_SEP = ';'


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
    names = []
    for iz in range(nz):
        filename = '{0}_{1}'.format(prefix, str(iz).zfill(nDecimals))
        names.append(filename)
        csvFilename = filename + '.csv'
        if not path.isfile(csvFilename):
            exit('File {0} not found, exiting'.format(filename))

        with open(csvFilename, 'r') as f:
            z = float(f.readline().split('=')[1])
            zvalues[iz] = z
            f.readline()

            for i in range(nxy):
                slices[iz][i] = readCsvLine(f)

    return zvalues, names, slices


def makePlots(dimensions, parsedData):
    nIntervals, xyzMin, xyzMax = dimensions
    zvalues, names, slices = parsedData


if len(argv) < 2:
    exit('usage: python plot.py <prefix>')

prefix = argv[1]
dimensionFileName = prefix + '_dimensions.csv'
if not path.isfile(dimensionFileName):
    exit('No dimension file found for given prefix, exiting')

dimensions = readDimensions(dimensionFileName)
parsedData = parseData(prefix, dimensions)

makePlots(dimensions, parsedData)
