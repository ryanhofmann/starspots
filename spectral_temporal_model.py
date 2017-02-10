"""
This package is used to plot the relative flux of a spotted star,
resolved both spectrally (y-axis) and temporally (x-axis).

Before execution, define the spot parameters in this file.

When executed, a model is created with the specified spots.
LDTK is used to calculate the nonlinear LD coefficients, for
both the star and the spot, as a function of wavelength.
The appropriate PHOENIX spectra are used to calculate the
wavelength-dependent contrast ratio of the starspots.
The lightcurve for each wavelength bin is then calculated
using MACULA, with the corresponding LD coefficients and
contrast ratio.

The relative flux is then plotted as a 2D heatmap.
"""

import pymacula
import numpy as np
import ldtk
import line_ratio
import matplotlib.pyplot as plt

class model(object):
  def __init__(self, Tphot=5800, Tspot=3800, spots=None, nspots=1):
    """
    Class containing temperatures and spot parameters.
    Creates a MaculaModel with specified spot parameters and default star.
    """

    self.Tphot = Tphot
    self.Tspot = Tspot

    if spots is None:
      self.spots = [pymacula.Spot() for i in range(nspots)]
    elif type(spots) is dict:
      self.spots = [pymacula.Spot(**spots) for i in range(nspots)]
    else:
      self.spots = spots
    self.nspots = len(self.spots)

    self.macmod = pymacula.MaculaModel(spots=spots)


def ldParams(T=5800, logg=4.50, z=0.0, xmin=400, xmax=700, dx=10):
  """
  Computes nonlinear LD parameters for wavelength range (xmin, xmax, dx; in nm).
  Returns coefficients as 2D array.
  """

  # Define passbands
  filters = []
  for i in range(int((xmax-xmin)/dx)):
    filters.append(ldtk.BoxcarFilter('{:d}'.format(i), xmin + dx*i, xmin + dx*(i+1)))

  # Define star and download uncached spectra
  sc = ldtk.LDPSetCreator(teff=(T, 1), logg=(logg, 0.01), z=(z, 0.01), filters=filters)

  # Create LD profiles and estimate nonlinear coefficients
  ps = sc.create_profiles()
  cq, eq = ps.coeffs_nl(do_mc=False)

  # Return coefficients array
  return cq


if __name__=="__main__":

  # Set star properties
  Tphot = 5000
  logg = 4.50
  z = 0.0

  # Set spot properties
  Tspot = 4000
  nspots = 1
  lat = [50]
  lon = [0]
  ingress = [0.2]
  egress = [0.2]
  tmax = [0.5]
  lifetime = [0.4]
  spots = [pymacula.Spot(lat=lat[i]*np.pi/180., lon=lon[i]*np.pi/180., ingress=ingress[i], egress=egress[i], tmax=tmax[i], lifetime=lifetime[i]) for i in range(nspots)]

  # Set wavelength range and step value
  xmin, xmax, dx = 300, 2500, 10

  # Initialize model
  star = model(Tphot, Tspot, spots)

  # Compute limb darkening coefficients
  ld_star = ldParams(Tphot, logg, z, xmin, xmax, dx)
  ld_spot = ldParams(Tspot, logg, z, xmin, xmax, dx)

  # Get PHOENIX spectra
  preface = "PHOENIX_spectra/"
  spec1_full, spec2_full = line_ratio.getSpectra(preface, "R", Tphot, Tspot, logg, z)

  # Resample spectra with boxcar average
  x = np.exp(1e-5*np.arange(len(spec1_full)))*300
  nx = int((xmax-xmin)/dx)
  spec1 = np.zeros(nx)
  spec2 = np.zeros(nx)
  for i in range(nx):
    spec1[i] = np.mean(spec1_full[np.where(np.abs(x-(2*xmin+(2*i+1)*dx)/2) <= dx/2)])
    spec2[i] = np.mean(spec2_full[np.where(np.abs(x-(2*xmin+(2*i+1)*dx)/2) <= dx/2)])

  # Compute ratio of spectra
  ratio = spec2/spec1

  # Compute lightcurve for each spectral bin
  ts = np.arange(0, 500, 0.05)
  nt = int(500/0.05)
  flux = np.zeros([nx, nt])
  for i in range(nx):
    star.macmod.star.c1 = ld_star[i,0]
    star.macmod.star.c2 = ld_star[i,1]
    star.macmod.star.c3 = ld_star[i,2]
    star.macmod.star.c4 = ld_star[i,3]
    star.macmod.star.d1 = ld_spot[i,0]
    star.macmod.star.d2 = ld_spot[i,1]
    star.macmod.star.d3 = ld_spot[i,2]
    star.macmod.star.d4 = ld_spot[i,3]
    for j in range(nspots):
      star.macmod.spots[j].contrast = ratio[i]
    flux[i,:] = star.macmod(ts)

  # Plot results
  from matplotlib import cm
  from matplotlib import colorbar
  fig = plt.figure(figsize=(12,11))
  ax1 = plt.subplot2grid((10,9), (0,0), rowspan=2, colspan=8)
  plt.ylabel("relative flux")
  ax2 = plt.subplot2grid((10,9), (2,0), rowspan=4, colspan=8, sharex=ax1)
  plt.ylabel("wavelength [nm]")
  ax3 = plt.subplot2grid((10,9), (6,0), rowspan=4, colspan=8, sharex=ax1)
  plt.xlabel("time")
  plt.ylabel("wavelength [nm]")
  ax4 = plt.subplot2grid((10,9), (2,8), rowspan=4, colspan=1)
  ax5 = plt.subplot2grid((10,9), (6,8), rowspan=4, colspan=1)
  f1 = ax1.plot(ts, np.mean(flux, axis=0), '-k')
  f2 = ax2.imshow(flux, aspect='auto', extent=[0,500,xmin,xmax], origin='lower', cmap=cm.get_cmap('gray'))
  plt.colorbar(f2, cax=ax4)
  f3 = ax3.imshow(flux/np.mean(flux,axis=0), aspect='auto', extent=[0,500,xmin,xmax], origin='lower', cmap=cm.get_cmap('gray'))
  plt.colorbar(f3, cax=ax5)
  plt.tight_layout()
  plt.savefig("Tp_{:d}_Ts_{:d}.png".format(Tphot,Tspot))

