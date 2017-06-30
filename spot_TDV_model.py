#!/usr/bin/env python3

"""
This module is used to create a large number of starspot models
with a range of parameters, optionally constrained to a specified
level of brightness variation. The models are saved to disk, and
two plots are produced for each: a 3-panel plot, with the mean
lightcurve and several TDV curves adjacent to the full 2-D flux;
and a 2-panel plot, with the mean lightcurve and the full
distribution of TDV curves.

Before executing, define settings in "spot_TDV_model_settings.py".
"""

import numpy as np
import matplotlib.pyplot as plt
import spectral_temporal_model as stm
import pymacula
import os
import pickle
import time
from matplotlib import cm
from matplotlib import colorbar

if __name__=="__main__":

  # Get input values from settings file
  from spot_TDV_model_settings import *

  # Set wavelength range and step value
  xs = np.arange(xmin,xmax,dx)

  # Set time range and step value
  ts = np.arange(tmin,tmax,dt)

  # Create response function for specified bandpass
  if band == "Kepler":
    band_response = stm.keplerFilter(dx)
  elif band == "MEarth":
    band_response = stm.mearthFilter(fin="mearthtrans.txt")
  elif band == "None":
    pass
  else:
    print("Error: invalid bandpass specified")
    exit()


  # Create new folder if none exists
  if not os.path.isdir(foldername):
    os.mkdir(foldername)
    os.mkdir(foldername+'/fluxes')
    os.mkdir(foldername+'/plots')

  # For each simulation
  i = 0
  KMin, KMax = 0, 0
  TDVMin, TDVMax = 0, 0
  print("Computing fluxes")
  tstart = time.time()
  while i < N:
    try:
      # Set spot properties
      nspots = np.random.randint(1, N_max + 1)  # random integer between 1 and 10
      alpha_max = R_spot*np.ones(nspots)
      spots = [pymacula.Spot(alpha_max=alpha_max[j]) for j in range(nspots)]

      # Compute fluxes
      flux = stm.computeFlux(Tphot,Tspot,logg,z,nspots,spots,xmin,xmax,dx,tmin,tmax,dt,Peq)

      if band != "None":
        # Filter to selected bandpass
        if band == "Kepler":
          # Create Kepler lightcurve
          band_flux = np.zeros([46,int((tmax-tmin)/dt)])
          for j in range(46):
            band_flux[j] = band_response(430+dx*j)*flux[int((430-xmin)/dx)+j]
          band_curve = np.mean(band_flux,axis=0)/np.mean(band_flux[:,0])
        elif band == "MEarth":
          band_flux = np.zeros([42,int((tmax-tmin)/dt)])
          for j in range(42):
            band_flux[j] = band_response(650+dx*j)*flux[int((650-xmin)/dx)+j]
          band_curve = np.mean(band_flux,axis=0)/np.mean(band_flux[:,0])
      elif band == "None":
        band_curve = np.mean(flux,axis=0)/np.mean(flux[:,0])

      # Check band sigma
      sigma = np.var(band_curve)**0.5
      amp = (np.max(band_curve) - np.min(band_curve))/2.
      if filter_by == "amp":
        if amp < amp_min or amp > amp_max:
          continue
      if filter_by == "sigma":
        if sigma < sigma_min or sigma > sigma_max:
          continue

      # Update plot limits
      TDV = flux**-1*Rp_Rs**2
      TDV_lines = np.array([TDV[:, int((k+0.5)*tmax/dt/nlines)] for k in range(nlines)])
      Kmin = np.min(band_curve)
      Kmax = np.max(band_curve)
      TDVmin = np.min(TDV_lines)
      TDVmax = np.max(TDV_lines)
      if Kmax > KMax: KMax = Kmax
      if i == 0: KMin = KMax
      if Kmin < KMin: KMin = Kmin
      if TDVmax > TDVMax: TDVMax = TDVmax
      if i == 0: TDVMin = TDVMax
      if TDVmin < TDVMin: TDVMin = TDVmin

      # Save flux array
      fluxname = "{}/fluxes/{:03d}.pkl".format(foldername, i)
      params = (i, nspots, sigma, amp, KMin, KMax, TDVMin, TDVMax)
      pickle.dump((params, flux), open(fluxname, 'wb'))

      tnew = time.time()
      tleft = (tnew - tstart)/(i+1)*(N - (i+1))
      print("\r{:d}/{:d} fluxes computed, {:.0f} minutes remaining ".format(i+1, N, tleft/60.), end='')
      i += 1

    except KeyboardInterrupt:
      break

  print("\n{:d} fluxes computed".format(i))
  print("Creating plots")

  from tqdm import *
  fluxnames = os.listdir(foldername+'/fluxes')
  fluxnames.sort()
  final_params, flux = pickle.load(open(foldername+'/fluxes/'+fluxnames[-1], 'rb'))
  i, nspots, sigma, amp, KMin, KMax, TDVMin, TDVMax = final_params
  for j in tqdm(range(len(fluxnames))):

    # Load flux array
    params, flux = pickle.load(open(foldername+'/fluxes/'+fluxnames[j], 'rb'))
    i, nspots, sigma, amp, Kmin, Kmax, TDVmin, TDVmax = params

    # band_curve TDVdist plot
    if band=="Kepler":
      # Create Kepler lightcurve
      band_flux = np.zeros([46,int((tmax-tmin)/dt)])
      for k in range(46):
        band_flux[k] = band_response(430+dx*k)*flux[int((430-xmin)/dx)+k]
      band_curve = np.mean(band_flux,axis=0)/np.mean(band_flux[:,0])
    elif band=="MEarth":
      band_flux = np.zeros([42,int((tmax-tmin)/dt)])
      for k in range(42):
        band_flux[k] = band_response(650+dx*k)*flux[int((650-xmin)/dx)+k]
      band_curve = np.mean(band_flux,axis=0)/np.mean(band_flux[:,0])
    elif band=="None":
      band_curve = np.mean(flux,axis=0)/np.mean(flux[:,0])

    # Plot lightcurve
    plt.figure(figsize=(8,8))
    plt.subplot(211)
    plt.plot(ts, band_curve, '-k')
    plt.title("nspots={:d}, sigma={:.5f}, amp={:.5f}".format(nspots, sigma, amp))
    plt.ylim(KMin, KMax)
    plt.xlabel("time")
    plt.ylabel("{} flux".format(band))

    # Plot delta_obs distribution
    plt.subplot(212)
    TDV = flux**-1*Rp_Rs**2
    for k in range(int(len(TDV[0])/10)):
      plt.plot(xs, TDV[:, k*10], '-k', alpha=0.03)
    plt.axhline(Rp_Rs**2, color='k')
    plt.ylim(TDVMin, TDVMax)
    plt.xlabel("wavelength")
    plt.ylabel("transit depth")
    plt.tight_layout()

    # Save figure
    plt.savefig("{}/plots/{:d}spot_{:d}_{}curve_TDVdist.png".format(foldername, nspots, i, band))
    plt.clf()

    # 3 panel plot
    fig = plt.figure(figsize=(12, 8))
    ax1 = plt.subplot2grid((8,24), (0,6), rowspan=2, colspan=15)
    plt.ylabel("{} flux".format(band))
    ax2 = plt.subplot2grid((8,24), (2,0), rowspan=6, colspan=6)
    plt.ylabel("wavelength")
    plt.xlabel("transit depth")
    ax3 = plt.subplot2grid((8,24), (2,6), rowspan=6, colspan=15, sharex=ax1, sharey=ax2)
    plt.xlabel("time")
    ax4 = plt.subplot2grid((8,24), (2,21), rowspan=6, colspan=3)
    cmap = cm.get_cmap('rainbow')
    colors = [cmap((k+0.5)/nlines) for k in range(nlines)]
    f11 = ax1.plot(ts, band_curve, '-k')
    f12 = [ax1.axvline((k+0.5)*tmax/nlines, color=colors[k]) for k in range(nlines)]
    ax1.set_ylim([KMin, KMax])
    f21 = [ax2.plot(TDV[:, int((k+0.5)*tmax/dt/nlines)], xs, color=colors[k]) for k in range(nlines)]
    f22 = ax2.axvline(Rp_Rs**2, color='k')
    ax2.set_xlim([TDVMax, TDVMin])
    plt.setp(ax2.get_xticklabels(), rotation=40, horizontalalignment='right')
    f3 = ax3.imshow(flux/np.mean(flux,axis=0), aspect='auto', extent=[tmin,tmax,xmin,xmax], origin='lower', cmap=cm.get_cmap('plasma'))
    plt.colorbar(f3, cax=ax4)
    fig.suptitle("nspots={:d}, sigma={:.5f}, amp={:.5f}".format(nspots, sigma, amp))
    plt.tight_layout()
    fig.subplots_adjust(top=0.94)

    # Save figure
    plt.savefig("{}/plots/{:d}spot_{:d}_3panel.png".format(foldername, nspots, i))
    plt.clf()

