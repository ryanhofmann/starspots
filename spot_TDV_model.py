"""
Choose Tphot,Tspot.
For i in range(a lot):
  generate model(random spots, n=1-10)
  if sigma in Kepler bandpass = 1% +/- 0.1%:
    store range of delta_obs(lambda)
    3-panel plot:
      flux plot in center
      Kepler lightcurve on top
      delta_obs plot, rotated, at left, 10 evenly spaced times
        label times with vlines on flux plot
Plot all delta_obs(lambda) ranges
For 10 first models:
  plot flux vs time
  plot D vs lambda
Violin plot of D at set of lambda
  start with histograms of D
  start with every 100 nm (22 total histograms)
"""

import numpy as np
import matplotlib.pyplot as plt
import spectral_temporal_model as stm
import pymacula
import os
import pickle
from matplotlib import cm
from matplotlib import colorbar

if __name__=="__main__":

  # Number of models to try
  N = 100

  # Set star properties
  Tphot = 2560
  Tspot = 2300
  logg = 5.0
  z = 0.0
  Rs = 0.117
  Rp_Rs = (1./109.)/Rs

  # Set wavelength range and step value
  xmin, xmax, dx = 300, 2500, 10
  xs = np.arange(xmin,xmax,dx)

  # Set time range and step value
  tmin, tmax, dt = 0, 500, 0.05
  ts = np.arange(tmin,tmax,dt)

  # Create new folder if none exists
  foldername = "TDV/{:05d}_{:05d}_{:.3f}".format(Tphot, Tspot, Rp_Rs)
  if not os.path.isdir(foldername):
    os.mkdir(foldername)
    os.mkdir(foldername+'/fluxes')
    os.mkdir(foldername+'/plots')

  # Initialize RNG
  np.random.seed(170303)

  # For each simulation
  i = 0
  KMin, KMax = 0, 0
  TDVMin, TDVMax = 0, 0
  nlines = 20
  print("Computing fluxes")
  while i < N:
    try:
      # Set spot properties
      nspots = int(np.ceil(10*np.random.rand()))  # random integer between 1 and 10
      alpha_max = 10*np.ones(nspots)
      spots = [pymacula.Spot(alpha_max=alpha_max[j]) for j in range(nspots)]

      # Compute fluxes
      flux = stm.computeFlux(Tphot,Tspot,logg,z,nspots,spots,xmin,xmax,dx,tmin,tmax,dt)

      # Create Kepler lightcurve
      K_response = stm.keplerFilter(dx)
      K_flux = np.zeros([46,int((tmax-tmin)/dt)])
      for j in range(46):
        K_flux[j] = K_response(430+dx*j)*flux[int((430-xmin)/dx)+j]
      K_curve = np.mean(K_flux,axis=0)/np.mean(K_flux[:,0])

      # Check Kepler sigma
      sigma = np.var(K_curve)**0.5
      amp = (np.max(K_curve) - np.min(K_curve))/2.
#      if amp < 0.0075 or amp > 0.0125:
#        continue

      # Update plot limits
      TDV = flux**-1*Rp_Rs**2
      TDV_lines = np.array([TDV[:, int((k+0.5)*tmax/dt/nlines)] for k in range(nlines)])
      Kmin = np.min(K_curve)
      Kmax = np.max(K_curve)
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

      print("\r{:d}/{:d} fluxes computed".format(i+1, N), end='')
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

    # Kcurve TDVdist plot
    # Create Kepler lightcurve
    K_response = stm.keplerFilter(dx)
    K_flux = np.zeros([46,int((tmax-tmin)/dt)])
    for k in range(46):
      K_flux[k] = K_response(430+dx*k)*flux[int((430-xmin)/dx)+k]
    K_curve = np.mean(K_flux,axis=0)/np.mean(K_flux[:,0])

    # Plot lightcurve
    plt.figure(figsize=(8,8))
    plt.subplot(211)
    plt.plot(ts, K_curve, '-k')
    plt.title("nspots={:d}, sigma={:.5f}, amp={:.5f}".format(nspots, sigma, amp))
    plt.ylim(KMin, KMax)
    plt.xlabel("time")
    plt.ylabel("Kepler flux")

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
    plt.savefig("{}/plots/{:d}spot_{:d}_Kcurve_TDVdist.png".format(foldername, nspots, i))
    plt.clf()

    # 3 panel plot
    fig = plt.figure(figsize=(12, 8))
    ax1 = plt.subplot2grid((8,24), (0,6), rowspan=2, colspan=15)
    plt.ylabel("Kepler flux")
    ax2 = plt.subplot2grid((8,24), (2,0), rowspan=6, colspan=6)
    plt.ylabel("wavelength")
    plt.xlabel("transit depth")
    ax3 = plt.subplot2grid((8,24), (2,6), rowspan=6, colspan=15, sharex=ax1, sharey=ax2)
    plt.xlabel("time")
    ax4 = plt.subplot2grid((8,24), (2,21), rowspan=6, colspan=3)
    cmap = cm.get_cmap('rainbow')
    colors = [cmap((k+0.5)/nlines) for k in range(nlines)]
    f11 = ax1.plot(ts, K_curve, '-k')
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
