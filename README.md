## AO 2

* Under development.

### ```spectools``` package

* A simple, flexible Python package for basic spectral analysis.
* Dependencies: ```AstroPy```, ```ccdproc```, ```NumPy```, ```scikit-image```.
 * Developed on Python 3.4+

#### ```spectools.SpecData```

* Class to hold spectral data, and to analyze it.

##### ```SpecData.load_fits(path, hdu=0, unit='adu')```

* Load a fits file into the ```SpecData``` object. hdu=0 means primary FITS hdu. If FITS extension is used, hdu have to be set.

##### ```SpecData.load_2d(Data, unit='adu')```

* Load a 2-dimensional array into the ```SpecData``` object.

##### ```SpecData.reduce_axis(axis)```

* Reduce given axis of the ```SpecData``` object, and make an 1-dimensional array.

##### ```SpecData.get_length(axis)```

* Fetch the length of give axis of the ```SpecData``` object.

##### ```find_peaks(max_peaks=10, axis_target=0, size_slit=10, min_separation=5, sigma_clip=True, sigma_value=3, minamp_peak=0.01)```

* Find peaks from ```SpecData``` object.

#### ```spectools.RefLamps```

* Data of reference lamps.

##### ```RefLamps.neon```

* Neon reference lamp