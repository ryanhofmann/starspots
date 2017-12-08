# starspots
This repository contains tools for simulating the effects of unocculted starspots on transmission spectroscopy of exoplanets.
## Required Python packages
Before running any scripts, install [`pymacula`](https://github.com/timothydmorton/pymacula) (requires Fortran compiler) and [`ldtk`](https://github.com/hpparvi/ldtk), and make sure you have `numpy`, `matplotlib`, and `scipy` as well.
### Pymacula for Python 3
Pymacula is written for Python 2, and requires a slight modification to work with Python 3.
In the Pymacula install folder, open `pymacula/macula.py` and replace all occurences of `xrange` with `range`.
## PHOENIX spectra
All PHOENIX spectra grids can be downloaded [here](http://phoenix.astro.physik.uni-goettingen.de/?page_id=16). Spectra are downloaded as ZIP archives for each combination of grid resolution (angstroms or R, use R), metallicity, and alpha enhancement. Create a folder named `PHOENIX_spectra` in the working directory and extract the ZIP archive there, e.g. `PHOENIX_spectra/PHOENIX-ACES-AGSS-COND-2011_R10000FITS_Z-0.0`.
