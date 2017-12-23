#!/usr/bin/env python

import spectools
from ccdproc import CCDData

import matplotlib.pyplot as plt
import numpy as np



#read images
himalia1 = spectools.SpecData('Sample_Himalia1.fits')
himalia2 = spectools.SpecData('Sample_Himalia2.fits')

reference = spectools.SpecData('Sample_Hydrogen.fits')
neon = spectools.SpecData('Sample_Neon.fits')
bias = spectools.SpecData('Sample_Bias.fits')
star1 = spectools.SpecData('Sample_Star_76151.fits')
star2 = spectools.SpecData('Sample_Star_76446.fits')
star3 = spectools.SpecData('Sample_Star_77730.fits')



#bias calibration
#basic arithmetic operations are available for SpecData().
reference -= bias
himalia1 -= bias
himalia2 -= bias
star1 -= bias
star2 -= bias
star3 -= bias
neon -= bias

#finding peaks
identified, peaks, max_intensity = reference.find_peaks(max_peaks=4)

#plot the result
#x_plot = np.arange(0, len(identified))
#plt.plot(x_plot, identified)

#for i in peaks:
#	plt.plot((i, i), (identified[i]+0.01*max_intensity, identified[i]+0.05*max_intensity), color='r', ls='-', lw=1)

#plt.show()


identified_neon, peaks_neon, max_intensity_neon = neon.find_peaks(max_peaks=20)

#plot the result
#x_plot = np.arange(0, len(identified_neon))
#plt.plot(x_plot, identified_neon)

#for i in peaks_neon:
#	plt.plot((i, i), (identified_neon[i]+0.01*max_intensity_neon, identified_neon[i]+0.05*max_intensity_neon), color='r', ls='-', lw=1)
#plt.show()

#Gaussian fit
gauss = reference.fit_gauss(identified, peaks)
comparison = np.array(list(map(lambda x: (float(x[1]), float(peaks[x[0]])), enumerate(gauss))))

gauss_neon = neon.fit_gauss(identified_neon, peaks_neon)
comparison_neon = np.array(list(map(lambda x: (float(x[1]), float(peaks_neon[x[0]])), enumerate(gauss_neon))))


#plot the result
#x_plot = np.arange(0, len(identified))
#plt.plot(x_plot, identified)

#for i in peaks:
#	plt.plot((i, i), (identified[i]+0.01*max_intensity, identified[i]+0.05*max_intensity), color='r', ls='-', lw=1)
#plt.show()

#matching table of pixel value and wavelength
table = [
	[109.310905, 8634.648],
	[141.029251, 8495.36],
	[156.524582, 8410.4],
	[164.874344, 8377.607],
	[180.063324, 8300.326],
	[213.142975, 8136.4],
	[252.214828, 7943.3],
	[344.407532, 7488.872],
	[354.899384, 7438.899],
	[394.091492, 7245.167],
	[408.980316, 7173.939],
	[438.088776, 7032.4127],
	[459.093384, 6929.468],
	[503.110107, 6717.0428],
	[511.124329, 6678.2],
	[527.915466, 6598.9529],
	[541.708313, 6532.8824],
	[547.044922, 6506.5279],
	[568.953613, 6402.246],
	[583.027527, 6382.9914],
	[597.194397, 6334.4279],
	[607.705383, 6304.7892],
	[623.158020, 6266.495],
	[633.123962, 6217.2813],
	[637.707886, 6163.5939],
	[647.077393, 6029.9971],
	[665.066345, 5944.8342],
	[684.699707, 5852.4878]
]
pixels = []
wavelengths = []
for i in table:
	pixels.append(i[0])
	wavelengths.append(i[1])

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

reduced = himalia1.reduce_axis(1)
identified_himalia1, peaks_himalia1, max_intensity_himalia1 = himalia1.find_peaks(max_peaks=100, axis_target=1)
#TODO: full options.

#location on y axis of Himalia
gauss = himalia1.fit_gauss(identified_himalia1, peaks_himalia1)
himalia1_yloc = gauss[8] 

reduced = himalia2.reduce_axis(1)
identified_himalia2, peaks_himalia2, max_intensity_himalia2 = himalia2.find_peaks(max_peaks=100, axis_target=1)
#TODO: full options.

#location on y axis of Himalia
gauss = himalia2.fit_gauss(identified_himalia2, peaks_himalia2)
himalia2_yloc = gauss[11] 


#reduce axis to find peak of reference star
reduced_star1 = star1.reduce_axis(1)
identified_star1, peaks_star1, max_intensity_star1 = star1.find_peaks(max_peaks=1, axis_target=1)
reduced_star2 = star2.reduce_axis(1)
identified_star2, peaks_star2, max_intensity_star2 = star2.find_peaks(max_peaks=1, axis_target=1)
reduced_star3 = star3.reduce_axis(1)
identified_star3, peaks_star3, max_intensity_star3 = star3.find_peaks(max_peaks=2, axis_target=1)

#location on y axis of reference star
gauss = star1.fit_gauss(identified_star1, peaks_star1)
star1_yloc = gauss[0]
gauss = star2.fit_gauss(identified_star2, peaks_star2)
star2_yloc = gauss[0] 
gauss = star3.fit_gauss(identified_star3, peaks_star3)
star3_yloc = gauss[0] 

#plot to see the result
#x_plot = np.arange(0, len(identified_himalia2))
#plt.plot(x_plot, identified_himalia2)

#for i in peaks_himalia2:
#	plt.plot((i, i), (identified_himalia2[i]+0.01*max_intensity_himalia2, identified_himalia2[i]+0.05*max_intensity_himalia2), color='r', ls='-', lw=1)

#plt.show()


#mimic the behavior of IRAF APALL
star1_spectra = star1.linap(sum_axis=0, center=star1_yloc, size=10)
star2_spectra = star2.linap(sum_axis=0, center=star2_yloc, size=10)
star3_spectra = star3.linap(sum_axis=0, center=star3_yloc, size=10)

himalia1_spectra = himalia1.linap(sum_axis=0, center=himalia1_yloc, size=10)
himalia2_spectra = himalia2.linap(sum_axis=0, center=himalia2_yloc, size=10)

neon_test = neon.linap(sum_axis=0, center=125, size=10)
reference_spectra = reference.linap(sum_axis=0, center=130, size=10)

#normalized the results with mean values
mean_star1 = np.mean(star1_spectra)
mean_star2 = np.mean(star2_spectra)
mean_star3 = np.mean(star3_spectra)

mean_himalia1 = np.mean(himalia1_spectra)
mean_himalia2 = np.mean(himalia2_spectra)

star1_norm = star1_spectra / mean_star1
star2_norm = star2_spectra / mean_star2
star3_norm = star3_spectra / mean_star3

himalia1_norm = himalia1_spectra / mean_himalia1
himalia2_norm = himalia2_spectra / mean_himalia2

#plot to see the result

divided_spectra1 = himalia2_norm / star1_norm
divided_spectra2 = himalia2_norm / star2_norm
divided_spectra3 = himalia2_norm / star3_norm

#plt.plot(vals[139:864], star1_norm[139:864], label='HD 76151')
#plt.plot(vals[139:864], star2_norm[139:864], label='HD 76446')
#plt.plot(vals[139:864], star3_norm[139:864], label='HD 77730')

#plt.plot(vals[139:864], himalia1_norm[139:864], label='Image 1')
#plt.plot(vals[139:864], himalia2_norm[139:864], label='Image 2')

plt.plot(vals[139:864], divided_spectra1[139:864], label='HD 76151')
plt.plot(vals[139:864], divided_spectra2[139:864], label='HD 76446')
plt.plot(vals[139:864], divided_spectra3[139:864], label='HD 77730')



#plt.plot(vals[139:864], divided_spectra[139:864])
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Relative flux')
plt.title('Flux of Himalia on the second image compared to reference stars')
#plt.title('Normalized flux of reference stars')
plt.legend(loc=0)
plt.show()
