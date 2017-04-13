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

  # Initialize RNG
  np.random.seed(170303)

  # For each simulation
  from tqdm import *
  for i in tqdm(range(N)):

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
      if amp < 0.0075 or amp > 0.0125:
        continue

      # Plot lightcurve
      plt.figure(figsize=(8,10))
      plt.subplot(311)
      plt.plot(ts, K_curve, '-k')
      plt.title("nspots={:d}, sigma={:.5f}, amp={:.5f}".format(nspots, sigma, amp))
      plt.xlabel("time")
      plt.ylabel("Kepler flux")

      # Plot delta_obs mean and sigma
      plt.subplot(312)
      TDV = flux**-1*Rp_Rs**2
      TDV_mean = np.mean(TDV, axis=1)
      TDV_std = np.std(TDV, axis=1)
      TDV_K_mean = np.mean(K_curve)
      TDV_K_std = np.std(K_curve)
      plt.plot(xs, TDV_mean, '-k')
      plt.axhline(Rp_Rs**2, color='k')
      plt.fill_between(xs, TDV_mean-TDV_std, TDV_mean+TDV_std, color='gray')
      plt.xlabel("wavelength")
      plt.ylabel("transit depth")
      plt.subplot(313)
      plt.hist(TDV[0], bins=50, color='gray')
      plt.xlabel("transit depth")
      plt.ylabel("N at 300 nm")
      plt.tight_layout()

      # Save figure
      plt.savefig("spot_TDV_testing/{:d}_Kcurve_TDV.png".format(i))
      plt.clf()

    except KeyboardInterrupt:
      break

