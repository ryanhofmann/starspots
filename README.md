# starspots
This repository contains tools for simulating the effects of unocculted starspots on transmission spectroscopy of exoplanets.
## Required Python packages
Before running any scripts, install `pymacula` and `ldtk`, and make sure you have `numpy`, `matplotlib`, and `scipy` as well.
### Pymacula for Python 3
Pymacula is written for Python 2, and requires a slight modification to work with Python 3.
In the Pymacula install folder, open `pymacula/macula.py` and replace all occurences of `xrange` with `range`.
## PHOENIX spectra
All PHOENIX spectra grids can be downloaded [here](http://phoenix.astro.physik.uni-goettingen.de/?page_id=16).
