import numpy as np

data = np.loadtxt('beam_intensities.txt')

normalized_data = [ relativeIntensity / (i + 1) for i, relativeIntensity in enumerate(data) ]
normalized_data /= sum(normalized_data)
print(normalized_data)
