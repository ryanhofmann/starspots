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
import pickle
import os
import numpy as np
import ldtk
import line_ratio
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colorbar
from scipy.interpolate import interp1d

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

  # If previously computed profile exists, load from pickle
  fname = "ldParams/T{:05d}_g{:.2f}_z{:.1f}_xmin{:d}_xmax{:d}_dx{:d}.pkl".format(T, logg, z, xmin, xmax, dx)
  if os.path.isfile(fname):
    cq = pickle.load(open(fname,"rb"))
    return cq

  # Define passbands
  filters = []
  for i in range(int((xmax-xmin)/dx)):
    filters.append(ldtk.BoxcarFilter('{:d}'.format(i), xmin + dx*i, xmin + dx*(i+1)))

  # Define star and download uncached spectra
  sc = ldtk.LDPSetCreator(teff=(T, 50), logg=(logg, 0.20), z=(z, 0.05), filters=filters)

  # Create LD profiles and estimate nonlinear coefficients
  ps = sc.create_profiles()
  cq, eq = ps.coeffs_nl(do_mc=False)

  # Save profile to pickle for easy recall
  if not os.path.isdir("./ldParams"):
    os.mkdir("./ldParams")
  pickle.dump(cq, open(fname,"wb"))

  # Return coefficients array
  return cq


def computeFlux(Tphot,Tspot,logg,z,nspots,spots,xmin,xmax,dx,tmin,tmax,dt):
  """
  Computes a MACULA light curve for each wavelength in x.
  Returns 2D flux array.
  """

  # Initialize model
  star = model(Tphot, Tspot, spots)

  # Get PHOENIX spectra, interpolating if Tphot or Tspot not in PHOENIX grid points
  preface = "PHOENIX_spectra/"
  Temps = np.concatenate((np.arange(2300,7000,100), np.arange(7000,12200,200)))
  if Tphot in Temps and Tspot in Temps:
    spec1_full, spec2_full = line_ratio.getSpectra(preface, "R", Tphot, Tspot, logg, z)
  else:
    Tp1, Tp2 = np.floor(Tphot/100.)*100, np.ceil(Tphot/100.)*100
    Ts1, Ts2 = np.floor(Tspot/100.)*100, np.ceil(Tspot/100.)*100
    sp1, sp2 = line_ratio.getSpectra(preface, "R", int(Tp1), int(Tp2), logg, z)
    ss1, ss2 = line_ratio.getSpectra(preface, "R", int(Ts1), int(Ts2), logg, z)
    if Tp1 == Tp2:
      spec1_full = sp1
    else:
      Tpx = (Tphot - Tp1)/(Tp2 - Tp1)
      spec1_full = Tpx*sp2 + (1.-Tpx)*sp1
    if Ts1 == Ts2:
      spec2_full = ss1
    else:
      Tsx = (Tspot - Ts1)/(Ts2 - Ts1)
      spec2_full = Tsx*ss2 + (1.-Tsx)*ss1

  # Resample spectra with boxcar average
  x = np.exp(1e-5*np.arange(len(spec1_full)))*300
  nx = int((xmax-xmin)/dx)
  spec1 = np.zeros(nx)
  spec2 = np.zeros(nx)
  for i in range(nx):
    spec1[i] = np.mean(spec1_full[np.where(np.abs(x-xmin-(2*i+1)*dx/2) <= dx/2)])
    spec2[i] = np.mean(spec2_full[np.where(np.abs(x-xmin-(2*i+1)*dx/2) <= dx/2)])

  # Compute ratio of spectra
  ratio = spec2/spec1

  # Compute limb darkening coefficients
  ld_star = ldParams(Tphot, logg, z, xmin, xmax, dx)
  ld_spot = ldParams(Tspot, logg, z, xmin, xmax, dx)

  # Compute lightcurve for each spectral bin
  ts = np.arange(tmin, tmax, dt)
  nt = int((tmax-tmin)/dt)
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

  return flux


def plotFull(flux,xmin,xmax,dx,tmin,tmax,dt,Tphot,Tspot):
  """
  Plots the average light curve, full 2D flux array, and normalized flux array.
  Saves plot to png.
  """

  ts = np.arange(tmin,tmax,dt)
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
  f2 = ax2.imshow(flux, aspect='auto', extent=[tmin,tmax,xmin,xmax], origin='lower', cmap=cm.get_cmap('gray'))
  plt.colorbar(f2, cax=ax4)
  f3 = ax3.imshow(flux/np.mean(flux,axis=0), aspect='auto', extent=[tmin,tmax,xmin,xmax], origin='lower', cmap=cm.get_cmap('gray'))
  plt.colorbar(f3, cax=ax5)
  plt.tight_layout()
  plt.savefig("Tp_{:.0f}_Ts_{:.0f}.png".format(Tphot,Tspot))
  plt.clf()


def keplerFilter(dx):
  """
  Two-component linear fit to Kepler response function.
  Returns interpolated function.
  """

  m1, b1, m2, b2 = .002, -.41, -.0015, 1.62  # linear coefficients for Kepler function
  K_wl = np.arange(430, 880+dx, dx)
  K_response = np.concatenate((m1*K_wl[:15]+b1, m2*K_wl[15:]+b2))
  K_func = interp1d(K_wl, K_response, bounds_error=False, fill_value=0)

  return K_func


def plotFilters(flux,xmin,xmax,dx,tmin,tmax,dt,Tphot,Tspot):
  """
  Plots normalized flux array and light curves in different filters.
  Filters used: Kepler, mean, 1.1-1.7 microns.
  """

  mean_curve = np.mean(flux,axis=0)
  K_response = keplerFilter(dx)
  K_flux = np.zeros([46,int((tmax-tmin)/dt)])
  for i in range(46):
    K_flux[i] = K_response(430+dx*i)*flux[int((430-xmin)/dx)+i]
  K_curve = np.mean(K_flux,axis=0)/np.mean(K_flux[:,0])
  J_curve = np.mean(flux[int((1100-xmin)/dx):int((1700-xmin)/dx)+1], axis=0)
  fig = plt.figure(figsize=(12,9))
  ax1 = plt.subplot2grid((7,9), (0,0), rowspan=3, colspan=8)
  plt.ylabel("relative flux")
  ax2 = plt.subplot2grid((7,9), (3,0), rowspan=4, colspan=8, sharex=ax1)
  plt.xlabel("time")
  plt.ylabel("wavelength [nm]")
  ax3 = plt.subplot2grid((7,9), (3,8), rowspan=4, colspan=1)
  f1 = ax1.plot(ts, mean_curve, '-g', ts, K_curve, '-b', ts, J_curve, '-r')
  f2 = ax2.imshow(flux/mean_curve, aspect='auto', extent=[tmin,tmax,xmin,xmax], origin='lower', cmap=cm.get_cmap('gray'))
  plt.colorbar(f2, cax=ax3)
  plt.tight_layout()
  plt.savefig("Tp_{:.0f}_Ts_{:.0f}_filters.png".format(Tphot,Tspot))
  plt.clf()


def plotTDV(flux,xmin,xmax,dx,tmin,tmax,dt,Tphot,Tspot,Rp_Rs):
  """
  Plots transit depth variations over 1/4 rotation, from spot ingress to spot transit.
  """

  TdV = flux**-1*Rp_Rs**2
  wl = np.arange(xmin,xmax,dx)
  mean_curve = np.mean(flux, axis=0)
  mins = np.where(np.r_[True, mean_curve[1:] < mean_curve[:-1]] & np.r_[mean_curve[:-1] < mean_curve[1:], True] == True)[0]
  ind_min = mins[int(len(mins)/2.+1)]
  p = mins[int(len(mins)/2.)] - mins[int(len(mins)/2.-1)]
  times = np.linspace(ind_min-0.25*p, ind_min, num=9, endpoint=True)
  for time in times:
    plt.plot(wl, TdV[:,int(time)])
  plt.axhline(Rp_Rs**2, color='k')
  plt.xlabel("wavelength [nm]")
  plt.ylabel(r"$\delta_{obs}$ [normalized flux]")
  plt.tight_layout()
  plt.savefig("Tp_{:.0f}_Ts_{:.0f}_RpRs_{:.4f}_TDV.png".format(Tphot,Tspot,Rp_Rs))
  plt.clf()


if __name__=="__main__":

  # Set star properties
  Tphot = 5770
  logg = 4.50
  z = 0.0
  Rs = 1.0

  # Set spot properties
  Tspot = 4770
  nspots = 1
  lat = [50.]
  lon = [0.]
  alpha_max = [5.]  # angular radius in degrees
  ingress = [0.2]
  egress = [0.2]
  tmax = [0.5]
  lifetime = [0.4]
  spots = [pymacula.Spot(lat=lat[i]*np.pi/180., lon=lon[i]*np.pi/180., alpha_max=alpha_max[i], ingress=ingress[i], egress=egress[i], tmax=tmax[i], lifetime=lifetime[i]) for i in range(nspots)]

  # Set wavelength range and step value
  xmin, xmax, dx = 300, 2500, 10
  xs = np.arange(xmin,xmax,dx)

  # Set time range and step value
  tmin, tmax, dt = 0, 500, 0.05
  ts = np.arange(tmin,tmax,dt)

  # Compute fluxes
  flux = computeFlux(Tphot,Tspot,logg,z,nspots,spots,xmin,xmax,dx,tmin,tmax,dt)

  # Plot average light curve, 2D, normalized 2D
  plotFull(flux,xmin,xmax,dx,tmin,tmax,dt,Tphot,Tspot)

  # Plot average, Kepler, and 1.1-1.7 micron light curves above 2D
  # Approximate Kepler with piecewise-linear fit, 430-880 nm
  plotFilters(flux,xmin,xmax,dx,tmin,tmax,dt,Tphot,Tspot)

  # Plot transit depth variations for multiple phases
  Rp_Rs = (1./109.)/Rs
  plotTDV(flux,xmin,xmax,dx,tmin,tmax,dt,Tphot,Tspot,Rp_Rs)
