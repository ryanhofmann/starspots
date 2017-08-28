#!/usr/bin/env python3

"""
This script creates a 2x2 grid of subplots illustrating the
wavelength-dependent effects of limb-darkening and spot contrast
on an observed spot transit.
"""

import numpy as np
import matplotlib.pyplot as plt
import pymacula
from ldtk import (BoxcarFilter, LDPSetCreator)
import line_ratio as lr

# Set star parameters
Tphot = 5800
logg = 4.5
z = 0.0

# Set spot parameters
Tspot1 = 5300
Tspot2 = 3800
lat = np.pi/4.
lon = 0
ingress = 0.2
egress = 0.2
lifetime = 0.4
tmax = 0.5
alpha_max = 5

# Set bandpass filters
blue, green, red = 350, 550, 750
filters = [BoxcarFilter("blue", blue-5, blue+5),
           BoxcarFilter("green", green-5, green+5),
           BoxcarFilter("red", red-5, red+5)]

# Get contrast ratios from spectra
ratio_1 = lr.specRatio("PHOENIX_spectra/", "R", Tphot, Tspot1, logg, z)
ratio_2 = lr.specRatio("PHOENIX_spectra/", "R", Tphot, Tspot2, logg, z)

# Determine contrast ratios in selected filters
xs = np.exp(1e-5*np.arange(len(ratio_1)))*300
contrast_b1 = np.mean(ratio_1[np.where(np.abs(xs-blue) <= 5)])
contrast_g1 = np.mean(ratio_1[np.where(np.abs(xs-green) <= 5)])
contrast_r1 = np.mean(ratio_1[np.where(np.abs(xs-red) <= 5)])
contrast_b2 = np.mean(ratio_2[np.where(np.abs(xs-blue) <= 5)])
contrast_g2 = np.mean(ratio_2[np.where(np.abs(xs-green) <= 5)])
contrast_r2 = np.mean(ratio_2[np.where(np.abs(xs-red) <= 5)])
contrasts = np.array([contrast_b1, contrast_g1, contrast_r1,
                      contrast_b2, contrast_g2, contrast_r2])

# Compute limb darkening coefficients
sc0 = LDPSetCreator(filters=filters, teff=(Tphot, 50), logg=(logg, 0.10), z=(z, 0.05))
sc1 = LDPSetCreator(filters=filters, teff=(Tspot1, 50), logg=(logg, 0.10), z=(z, 0.05))
sc2 = LDPSetCreator(filters=filters, teff=(Tspot2, 50), logg=(logg, 0.10), z=(z, 0.05))
ps0 = sc0.create_profiles(nsamples=500)
ps1 = sc1.create_profiles(nsamples=500)
ps2 = sc2.create_profiles(nsamples=500)
qc0,qe0 = ps0.coeffs_nl(do_mc=False)
qc1,qe1 = ps1.coeffs_nl(do_mc=False)
qc2,qe2 = ps2.coeffs_nl(do_mc=False)
qcs = np.array([qc0, qc1, qc2])

# Create starspot models
spot = [pymacula.Spot(lat=lat, lon=lon, alpha_max=alpha_max, ingress=ingress,
                     egress=egress, lifetime=lifetime, tmax=tmax)]
import copy
ref_b = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
ref_g = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
ref_r = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
ld_b = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
ld_g = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
ld_r = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
con_b = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
con_g = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
con_r = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
both_b = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
both_g = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
both_r = pymacula.MaculaModel(nspots=1, spots=copy.deepcopy(spot))
models = np.array([ref_b, ref_g, ref_r, ld_b, ld_g, ld_r,
                   con_b, con_g, con_r, both_b, both_g, both_r])

# Set contrast ratios
for i in range(6):
  models[i].spots[0].contrast = contrasts[i%3]
for i in range(6,12):
  models[i].spots[0].contrast = contrasts[i%3 + 3]

# Set limb darkening coefficients
for i in range(3,6):
  models[i].star.c1 = qcs[0, i%3, 0]
  models[i].star.c2 = qcs[0, i%3, 1]
  models[i].star.c3 = qcs[0, i%3, 2]
  models[i].star.c4 = qcs[0, i%3, 3]
  models[i].star.d1 = qcs[1, i%3, 0]
  models[i].star.d2 = qcs[1, i%3, 1]
  models[i].star.d3 = qcs[1, i%3, 2]
  models[i].star.d4 = qcs[1, i%3, 3]
for i in range(9,12):
  models[i].star.c1 = qcs[0, i%3, 0]
  models[i].star.c2 = qcs[0, i%3, 1]
  models[i].star.c3 = qcs[0, i%3, 2]
  models[i].star.c4 = qcs[0, i%3, 3]
  models[i].star.d1 = qcs[2, i%3, 0]
  models[i].star.d2 = qcs[2, i%3, 1]
  models[i].star.d3 = qcs[2, i%3, 2]
  models[i].star.d4 = qcs[2, i%3, 3]

# Plot reference light curves
ts = np.arange(0,500,0.05)
plt.subplot(221)
for i in range(3):
  plt.plot(ts, models[i](ts))
plt.xlim(235,265)
plt.ylim(0.994,1.0005)
plt.title("Low contrast, no LD")

# Plot limb darkening light curves
plt.subplot(222)
for i in range(3,6):
  plt.plot(ts, models[i](ts))
plt.xlim(235,265)
plt.ylim(0.994,1.0005)
plt.title("Low contrast, LD")

# Plot contrast light curves
plt.subplot(223)
for i in range(6,9):
  plt.plot(ts, models[i](ts))
plt.xlim(235,265)
plt.ylim(0.994,1.0005)
plt.title("High contrast, no LD")

# Plot combo light curves
plt.subplot(224)
for i in range(9,12):
  plt.plot(ts, models[i](ts))
plt.xlim(235,265)
plt.ylim(0.994,1.0005)
plt.title("High contrast, LD")

# Format plot
plt.suptitle("Tphot={:d}K, Tspot1={:d}K, Tspot2={:d}K,\nred={:d}nm, green={:d}nm, blue={:d}nm".format(Tphot, Tspot1, Tspot2, red, green, blue))
plt.tight_layout(rect=[0, 0.03, 1, 0.92])

# Save plot
plt.savefig("contrast_LD_grid_plot.png")
plt.savefig("contrast_LD_grid_plot.eps")
plt.clf()

