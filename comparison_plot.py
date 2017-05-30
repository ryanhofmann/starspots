#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pickle

# Set variables
flux_file = "/home/ryan/research/starspots/TDV/GJ1132b/fluxes/005.pkl"
#Rs = 1.
#Rp_Rs = (1./109.)/Rs
#Rp_Rs = 6.4/81.9
solar_file = "/home/ryan/research/starspots/Exo_Transmit/Spectra/GJ1132/500K_solar_GJ1132b.dat"
H2O_file = "/home/ryan/research/starspots/Exo_Transmit/Spectra/GJ1132/500K_H2O_GJ1132b.dat"
CO2_file = "/home/ryan/research/starspots/Exo_Transmit/Spectra/GJ1132/500K_CO2_GJ1132b.dat"
labels = ["solar", "H2O", "CO2"]
star = "GJ1132b"

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
xmin, xmax, dx = 3000, 25000, 100
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
plt.ylim((0.255, 0.295))
plt.legend()
plt.xlabel("angstroms")
plt.ylabel("transit depth [%]")
plt.title("Starspot TDVs and atmospheric models for {}".format(star))
plt.savefig("comparison_{}.png".format(star))

