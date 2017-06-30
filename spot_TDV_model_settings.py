#!/usr/bin/env python3

import numpy as np

# Number of models to try
N = 100

# Set star properties
Tphot = 3270  # photosphere temperature
Tspot = Tphot - 800  # starspot temperature (same for all spots)
N_max = 10  # maximum number of starspots
R_spot = 5  # maximum starspot radius in degrees
logg = 5.0  # log surface gravity
z = 0.0  # metallicity
Peq = 30  # equatorial rotation period in days
Rs = 0.186  # stellar radius in solar units
Rp_Rs = 1.43*(1./109.)/Rs  # R_P/R_Star

# Set wavelength range and step value
xmin, xmax, dx = 300, 2500, 10  # nanometers

# Set time range and step value
tmin, tmax, dt = 0, 500, 0.05

# Set bandpass
band = "None"  # Can be "Kepler" or "MEarth", or "None" for no filtering

# Set amplitude/sigma of brightness modulations
filter_by = "sigma"  # can be "amp" or "sigma"
amp_min, amp_max = 0.0075, 0.0125
sigma_min, sigma_max = 0.0025, 0.0175

# Create new folder if none exists
foldername = "{:05d}_{:05d}_{:.3f}".format(Tphot, Tspot, Rp_Rs)  # default naming convention
# foldername = "TDV/LHS1140b_100K"

# Initialize RNG (comment out these lines for nonrepeatable results)
import numpy as np
np.random.seed(170303)

# Number of TDV curves for 3-panel plot
nlines = 20

