#!/usr/bin/env python

import spectools
from ccdproc import CCDData


#read images
image = spectools.SpecData('GG435_Himalia_Spectroscopy_OBJECT_NaCS0044108.fits')
reference = spectools.SpecData('GG435_Neon Lamp_Spectroscopy_COMPARISON_NaCS0043900.fits')
bias = spectools.SpecData('masterbias_-100.0_3.fits')

#bias calibration
#basic arithmetic operations are available for SpecData()!
reference -= bias

peaks = reference.find_peaks()