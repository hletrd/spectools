#!/usr/bin/env python

from ccdproc import CCDData
import numpy as np

from astropy.io import fits
from astropy.modeling.models import Linear1D
from astropy import stats

from skimage.feature import peak_local_max

if __name__ == "__main__":
	pass

class SpecData(CCDData):
	err_str = {
		'no_data': 'No data exists!',
		'unsupported_axis': 'invalid axis values!',
		'dimension': 'Dimensions of the array do not match!'
	}

	class Fitter:
		def __init__(self):
			self._const_Gaussian = 0
			pass

		@property
		def Gaussian(self):
			return self._const_Gaussian

		@Gaussian.setter
		def Gaussian(self, val):
			pass

		from astropy.modeling.fitting import LevMarLSQFitter
		from astropy.modeling.models import Gaussian1D, Chebyshev2D

		from numpy.polynomial.chebyshev import chebfit, chebval

		@staticmethod
		def fit_gauss(self, data, identified, peaks, fwhm=4):
			peak_fitted = []
			fitter = LevMarLSQFitter()
			for peak in peaks:
				g_init = Gaussian1D(
					amplitude=identified,
					mean=peak,
					stddev=fwhm * stats.gaussian_fwhm_to_sigma,
					bounds={
						'amplitude': (0, 2*identified[peak]),
						'mean':(peak-fwhm, peak+fwhm),
						'stddev':(0, fwhm)
						})
				identified_x = np.arange(0, len(identified))
				fitted = LINE_FITTER(g_init, identified_x, identified)
				peak_fitted.append(fitted.mean.value)

			return peak_fitted

		@staticmethod
		def fit_chebyshev(self, pixels, wavelengths, deg=4):
			if len(pixels) != len(wavelengths):
				SpecData.raise_error(self, SpecData.err_str['dimension'])
			cnt_data = len(pixels)
			coeff_id, fitted_full = chebfit(
				pixels, 
				wavelengths, 
				deg=deg,
				full=True
				)
			fitted_rms = np.sqrt(fitted_full[0][0]/cnt_data)
			rough_error = np.ptp(wavelengths) / np.ptp(pixels) / 2
			residual = (wavelengths - chebval(pixels, coeff_id))
			range_residual = np.max(np.abs(residual))

			return {
				'pixels': pixels,
				'residual': residual,
				'rough_error': rough_error,
				'range_residual': range_residual,
			}

	@staticmethod
	def raise_error(self, err):
		raise Exception('SpecData Error', err)

	def __init__(self, path='', ccddata=[], hdu=0, unit='adu'):
		self.ccddata = []
		if ccddata != []:
			self.ccddata = CCDData(ccddata, unit=unit)
		if path != '':
			self.load_fits(path, hdu=hdu, unit=unit)

	#overloaded operators
	def __add__(self, image):
		if self.ccddata == []:
			SpecData.raise_error(self.err_str['no_data'])
		self.ccddata = CCDData(self.ccddata.data + image.data, unit='adu')
		return self

	def __sub__(self, image):
		if self.ccddata == []:
			SpecData.raise_error(self.err_str['no_data'])
		self.ccddata = CCDData(self.ccddata.data - image.data, unit='adu')
		return self

	def __mul__(self, image):
		if self.ccddata == []:
			SpecData.raise_error(self.err_str['no_data'])
		self.ccddata = CCDData(self.ccddata.data * image.data, unit='adu')
		return self

	def __truediv__(self, image):
		if self.ccddata == []:
			SpecData.raise_error(self.err_str['no_data'])
		self.ccddata = CCDData(self.ccddata.data / image.data, unit='adu')
		return self

	#mimic extended CCDData
	@property
	def data(self):
		return self.ccddata.data

	#load_fits
	#hdu to handle FITS extensions.
	def load_fits(self, path, hdu=0, unit='adu'):
		image = fits.open(path)
		self.ccddata = CCDData(image[hdu].data, unit=unit)

	def load_2d(self, CCDData, unit='adu'):
		self.ccddata = CCDData(CCDData, unit=unit)

	#performs reduction of given axis.
	def reduce_axis(self, axis):
		if self.ccddata == []:
			SpecData.raise_error(self, SpecData.err_str['no_data'])
		self.ccddata = CCDData(self.ccddata.data.sum(axis=axis), unit='adu')

	def get_data(self):
		if self.ccddata == []:
			SpecData.raise_error(self, SpecData.err_str['no_data'])
		return self.ccddata

	def get_length(self, axis):
		if self.ccddata == []:
			SpecData.raise_error(self, SpecData.err_str['no_data'])
		return self.ccddata.data.shape[axis]

	def find_peaks(
		self,
		max_peaks=10, #max number of peaks to find
		axis_target=0, #axis to find peak
		size_slit=10, #size of the slit
		min_separation=5, #minimum separation between peaks
		sigma_clip=True, #if to perform sigma clip before finding peaks
		sigma_value=3, #sigma value to clip
		minamp_peak=0.01, #minimum amplitude to regard as peak
		):
		if self.ccddata == []:
			SpecData.raise_error(self, SpecData.err_str['no_data'])
		length_target = self.get_length(axis=axis_target)
		search_min = length_target//2 - size_slit//2 
		search_max = length_target//2 + size_slit//2

		if axis_target == 0:
			identified = np.median(self.ccddata[search_min:search_max, :], axis=axis_target)
		elif axis_target == 1:
			identified = np.median(self.ccddata[:, search_min:search_max], axis=axis_target)
		else:
			SpecData.raise_error(self, SpecData.err_str['unsupported_axis'])
			return

		identified_original = identified
		if sigma_clip == True:
			identified = stats.sigma_clip(identified, sigma=sigma_value, iters=5)
		max_intensity = np.max(identified)

		peaks = peak_local_max(
			identified,
			indices=True,
			num_peaks=max_peaks,
			min_distance=min_separation,
			threshold_abs=max_intensity * minamp_peak
			)

		return peaks

	#using CCDData of instance
	def fit(self, init, fitter=0):
		if self.ccddata == []:
			SpecData.raise_error(self, SpecData.err_str['no_data'])
		if fitter == self.Fitter.Gaussian:
			self.Fitter.fit_gauss(init)
