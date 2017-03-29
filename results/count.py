from sys import argv, exit
from os import makedirs, path, listdir

import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import myparser


PLOTDIR = 'plots'


def minmax(parsedData):
    _, _, _, _, _, slices = parsedData
    slices = ma.masked_values(slices, 0) * 1e-6
    return np.min(slices[:, :, :, 3]), np.max(slices[:, :, :, 3])


def makePlots(filename, dimensions, parsedData, sliceAx, cMin=None, cMax=None):
    nIntervals, xyzMin, xyzMax = dimensions
    xvalues, yvalues, zvalues, nums, t, slices = parsedData
    axValuesArr = [xvalues, yvalues, zvalues]

    idx_ = filename.rfind('_')
    prefix = filename[:idx_]
    plotdir = PLOTDIR + '_' + prefix

    if not path.isdir(plotdir):
        makedirs(plotdir)

    slices = ma.masked_values(slices, 0) * 1e-6

    if cMin == None:
        cMin = np.min(slices[:, :, :, 3])
    if cMax == None:
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

    X = np.linspace(xyzMin[ax1], xyzMax[ax1], nIntervals[ax1] + 1) * 1000
    Y = np.linspace(xyzMin[ax2], xyzMax[ax2], nIntervals[ax2] + 1) * 1000
    XX, YY = np.meshgrid(X, Y)

    idxpoint = filename.rfind('.')
    prefix1 = filename[:idxpoint]

    for i, z in list(enumerate(axvalues))[57:58]:#[1::2]:
        countPath = path.join(plotdir,
            '{0}_count_{1}_{2}.png'.format(prefix1, axNames[sliceAx], nums[i])
        )
        data = slices[i, :, :, 3]

        heatmap = plt.pcolormesh(XX, YY, data,
                                 norm=LogNorm(vmin=cMin, vmax=cMax),
                                 cmap=plt.get_cmap('inferno'))

        plt.xlabel('${0}$ (mm)'.format(axNames[ax1]))
        plt.ylabel('${0}$ (mm)'.format(axNames[ax2]))
        plt.xlim(xyzMin[ax1] * 1000, xyzMax[ax1] * 1000)
        plt.ylim(xyzMin[ax2] * 1000, xyzMax[ax2] * 1000)
        plt.gca().set_aspect('equal')
        plt.gca().set_facecolor((0, 0, 0))
        # plt.title('{0} = {1:.1f} mm'.format(axNames[sliceAx], z * 1000))
        plt.title('t = {0:.0f} ms'.format(t * 1000))
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = plt.colorbar(heatmap, drawedges=False, cax=cax)
        cb.set_label('Neutral density ($\\mathrm{cm}^{-3}$)')
        plt.savefig(countPath, bbox_inches='tight', dpi=150)
        plt.clf()

if __name__ == '__main__':
    if len(argv) < 2:
        exit('usage: python count.py <prefix> <slice axis x/y/z/all>')

    prefix = argv[1]

    sliceAxes = {'x': 0, 'y': 1, 'z': 2, 'all': 3}

    ax = argv[2].lower()
    if not ax in sliceAxes:
        exit('usage: python count.py <prefix> <slice axis x/y/z/all>')

    sliceAx = sliceAxes[ax]
    dimensions = myparser.readDimensions(prefix)
    parsedFiles = []
    cMin = float('inf')
    cMax = float('-inf')
    for filename in listdir('.'):
        if '_dimensions.csv' in filename or '.csv' not in filename or \
                not filename.startswith(prefix):
            continue

        print('Parsing file ' + filename)
        parsedData = myparser.parseData(filename, dimensions)
        curMin, curMax = minmax(parsedData)
        if curMin < cMin:
            cMin = curMin
        if curMax > cMax:
            cMax = curMax
        parsedFiles.append((filename, parsedData))
        print('min = {0}, max = {1}'.format(cMin, cMax))

    for filename, parsedData in parsedFiles:
        if sliceAx == 3:
            print('plotting all axes')
            for x in range(3):
                makePlots(filename, dimensions, parsedData, x, cMin, cMax)
        else:
            makePlots(filename, dimensions, parsedData, sliceAx, cMin, cMax)

