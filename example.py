#!/usr/bin/env python

import spectools
from ccdproc import CCDData

import matplotlib.pyplot as plt
import numpy as np



#read images
himalia = spectools.SpecData('Sample_Himalia2.fits')
reference = spectools.SpecData('Sample_Hydrogen.fits')
bias = spectools.SpecData('Sample_Bias.fits')
star = spectools.SpecData('Sample_Star_76151.fits') #76446


#bias calibration
#basic arithmetic operations are available for SpecData().
reference -= bias
himalia -= bias
star -= bias

#finding peaks
identified, peaks, max_intensity = reference.find_peaks(max_peaks=4)


#plot the result
#x_plot = np.arange(0, len(identified))
#plt.plot(x_plot, reference.identified)

#for i in peaks:
#	plt.plot((i, i), (identified[i]+0.01*max_intensity, identified[i]+0.05*max_intensity), color='r', ls='-', lw=1)

#plt.show()


#Gaussian fit
gauss = reference.fit_gauss(identified, peaks)
comparison = np.array(list(map(lambda x: (float(x[1]), float(peaks[x[0]])), enumerate(gauss))))


#plot the result
#x_plot = np.arange(0, len(identified))
#plt.plot(x_plot, reference.identified)

#for i in peaks:
#	plt.plot((i, i), (identified[i]+0.01*max_intensity, identified[i]+0.05*max_intensity), color='r', ls='-', lw=1)
#plt.show()

#matching table of pixel value and wavelength
pixels = [149.83549377324192, 285.55542938303705, 534.32462237244499, 897.82309722927607]
wavelengths = np.array([4101.74, 4340.47, 4861.33, 6562.852])

#fit to chebyshev model
chebyshev = reference.fit_chebyshev(pixels, wavelengths)
x = np.arange(0, reference.get_length(1))
vals = np.polynomial.chebyshev.chebval(x, chebyshev['coeff'])


#plot the result
#x_plot = np.arange(0, reference.get_length(0))
#plt.plot(pixels, wavelengths)
#plt.plot(x, vals)
#plt.show()


#reduce axis to find peak of Himalia
#TODO: now assuming that dispersion axis is exactly perpendicular to slit axis. Rotational transform should be added for more accurate analysis.
reduced = himalia.reduce_axis(1)
identified_himalia, peaks_himalia, max_intensity_himalia = himalia.find_peaks(max_peaks=100, axis_target=1)

gauss = himalia.fit_gauss(identified_himalia, peaks_himalia)
himalia_yloc = gauss[8] #location on y axis of Himalia


#reduce axis to find peak of reference star
reduced = star.reduce_axis(1)
identified_star, peaks_star, max_intensity_star = star.find_peaks(max_peaks=100, axis_target=1)

gauss = star.fit_gauss(identified_star, peaks_star)
star_yloc = gauss[0] #location on y axis of reference star

#plot to see the result
#x_plot = np.arange(0, len(identified_star))
#plt.plot(x_plot, star.identified)

#for i in peaks_star:
#	plt.plot((i, i), (identified_star[i]+0.01*max_intensity_star, identified_star[i]+0.05*max_intensity_star), color='r', ls='-', lw=1)

#plt.show()


#mimic the behavior of IRAF APALL
star_spectra = star.linap(sum_axis=0, center=star_yloc, size=10)
himalia_spectra = himalia.linap(sum_axis=0, center=himalia_yloc, size=10)

#normalized the results with mean values
mean_star = np.mean(star_spectra)
mean_himalia = np.mean(himalia_spectra)
star_norm = star_spectra / mean_star
himalia_norm = himalia_spectra / mean_himalia

#plot to see the result

divided_spectra = himalia_norm / star_norm
#plt.plot(himalia_norm)
#plt.plot(star_norm)
plt.plot(vals[373:1000], divided_spectra[373:1000])
plt.show()
