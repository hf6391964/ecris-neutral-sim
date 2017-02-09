import count
import myparser

sizes = [ 1000000, 2000000, 5000000, 10000000, 20000000, 50000000, 100000000, 200000000 ]

for size in sizes:
    prefix = 'test' + str(size)
    print('Plotting ' + prefix)
    for x in range(0, 3):
        dimensions = myparser.readDimensions(prefix)
        parsedData = myparser.parseData(prefix, dimensions)
        count.makePlots(prefix, dimensions, parsedData, x)

