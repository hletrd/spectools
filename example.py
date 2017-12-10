#!/usr/bin/env python

import spectools
from ccdproc import CCDData

import matplotlib.pyplot as plt
import numpy as np



#read images
image = spectools.SpecData('Sample_Himalia1.fits')
reference = spectools.SpecData('Sample_Neon.fits')
bias = spectools.SpecData('Sample_Bias.fits')


#bias calibration
#basic arithmetic operations are available for SpecData().
reference -= bias


#finding peaks
identified, peaks, max_intensity = reference.find_peaks(max_peaks=20)


#plot the result
x_plot = np.arange(0, len(identified))
plt.plot(x_plot, reference.identified)

for i in peaks:
	plt.plot((i, i), (identified[i]+0.01*max_intensity, identified[i]+0.05*max_intensity), color='r', ls='-', lw=1)

#plt.show()


#Gaussian fit
gauss = reference.fit(identified, peaks, spectools.SpecData.Fitter.Gaussian)

comparison = np.array(list(map(lambda x: (float(x[1]), float(peaks[x[0]])), enumerate(gauss))))
print(comparison)