#!/usr/bin/env python3

"""
This script generates a plot comparing the transit depth variations
(TDV) caused by starspots to those caused by a transiting planet's
atmosphere. First, compute flux models with spot_TDV_model and use
the generated plots to select a model (the number in the plot's
filename). Locate the .pkl file with that number; this is your
flux_file. Next, use Exo-Transmit to compute spectra for your
planet's atmosphere; edit the appropriate lines below to identify
your spectra files. To change the compositions listed below, use
search-and-replace. Finally, set the remaining variables and
execute the script.
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle

# Set variables
planet = "LHS1140b"
dT = 100
flux_file = "/home/ryan/research/starspots/TDV/LHS1140b_100K/fluxes/014.pkl"
solar_file = "/home/ryan/research/starspots/Exo_Transmit/Spectra/LHS1140/300K_solar_LHS1140b.dat"
H2O_file = "/home/ryan/research/starspots/Exo_Transmit/Spectra/LHS1140/300K_H2O_LHS1140b.dat"
CO2_file = "/home/ryan/research/starspots/Exo_Transmit/Spectra/LHS1140/300K_CO2_LHS1140b.dat"
labels = ["solar", "H2O", "CO2"]
xmin, xmax, dx = 3000, 25000, 100

# Load flux array from file
params, flux = pickle.load(open(flux_file, "rb"))

# Load Exo-Transmit spectra
solar = np.genfromtxt(solar_file, skip_header=2)
H2O = np.genfromtxt(H2O_file, skip_header=2)
CO2 = np.genfromtxt(CO2_file, skip_header=2)
spectra = np.array([solar, H2O, CO2])

# Convert flux to transit depth
Rp_Rs = np.sqrt(0.01*np.min(spectra[:,:,1]))
TDV = flux**-1*Rp_Rs**2

# Plot TDV distribution
xs = np.arange(xmin,xmax,dx)
plt.figure(figsize=(16,9))
for i in range(int(len(TDV[0])/10)):
  plt.plot(xs, TDV[:, i*10]*100, '-k', alpha=0.03)

# Overplot transmission spectra
for spectrum, label in zip(spectra, labels):
  plt.plot(spectrum[:,0]*1e10, spectrum[:,1], label=label)

# Format plot
plt.axhline(Rp_Rs**2*100, color='k')
plt.xlim((xmin, xmax))
plt.ylim((0.49, 0.55))
plt.legend()
plt.xlabel("angstroms")
plt.ylabel("transit depth [%]")
plt.title("Starspot TDVs and atmospheric models for {}, DeltaT = {:d} K, ".format(planet, dT))
plt.savefig("comparison_{}_{:d}K.png".format(planet, dT))

