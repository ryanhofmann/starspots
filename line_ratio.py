#!/usr/bin/python3

"""
Computes the ratio of two spectra at different temperatures.
Assumes med-res spectra extracted to working or child directory.
If run standalone, plots ratio vs. wavelength and exits.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def getSpectra(preface="", grid="R", Tphot=5800, Tspot=3800, g=4.5, Z=0.0, alpha=0.0):
  """
  Load spectra into arrays
  grid can be "A" for angstroms or "R" for resolving power
  See PHOENIX paper and site for valid T, g, Z, alpha
  """

  # Set folder from which to load spectra
  if alpha==0:
    folder_suffix = "/"
  else:
    folder_suffix = ".Alpha={:+1.2f}/".format(alpha)

  if Z<=0:
    Zstr = "-{:1.1f}".format(abs(Z))
  else:
    Zstr = "+{:1.1f}".format(Z)

  if grid=="A":
    folder = preface + "PHOENIX-ACES-AGSS-COND-2011_A1FITS_Z" + Zstr + folder_suffix
  elif grid=="R":
    folder = preface + "PHOENIX-ACES-AGSS-COND-2011_R10000FITS_Z" + Zstr + folder_suffix
  else:
    print("Error: invalid grid")
    exit(1)

  # Create file names
  if alpha==0:
    alpha_str = ""
  else:
    alpha_str = ".Alpha={:+1.2f}".format(alpha)

  f1 = folder + "lte{:05d}-{:1.2f}".format(Tphot, g) + Zstr + alpha_str + ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
  f2 = folder + "lte{:05d}-{:1.2f}".format(Tspot, g) + Zstr + alpha_str + ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"

  # Load spectra into arrays
  spec1 = fits.open(f1)[0].data
  spec2 = fits.open(f2)[0].data

  return spec1, spec2


def conv(spec, grid="R", res=10):
  """
  Smooth spectrum by convolving with gaussian
  grid : "A" for linear, "R" for logarithmic
  res is desired resolution in angstroms at 1 micron
  """

  from astropy.convolution import Gaussian1DKernel, convolve

  # Generate gaussian kernel for convolution
  sigma = 10*res
  gauss = Gaussian1DKernel(stddev=sigma)

  # Convolve spectrum with kernel
  spec_conv = convolve(spec, gauss, boundary='extend')

  return spec_conv


def specRatio(preface="", grid="R", Tphot=5800, Tspot=3800, g=4.5, Z=0.0, alpha=0.0, smooth=False, sigma=10):
  """
  Computes ratio of two spectra
  """

  # Read spectra into arrays
  spec1, spec2 = getSpectra(preface, grid, Tphot, Tspot, g, Z, alpha)

  # Apply smoothing if specified
  if smooth:
    spec1 = conv(spec1, grid=grid, res=sigma)
    spec2 = conv(spec2, grid=grid, res=sigma)

  # Compute ratio
  ratio = spec2/spec1

  return ratio


if __name__=="__main__":

  # Get inputs
  preface = input("Parent directory [./]: ") or ""
  grid = input("grid [R]: ") or "R"
  Tphot = input("Photosphere temperature [5800]: ") or '5800'
  Tphot = int(Tphot)
  Tspot = input("Spot temperature [3800]: ") or '3800'
  Tspot = int(Tspot)
  g = input("Log surface gravity [4.5]: ") or '4.5'
  g = float(g)
  Z = input("Metallicity [0.0]: ") or '0.0'
  Z = float(Z)
  alpha = input("Alpha enhancement [0.0]: ") or '0.0'
  alpha = float(alpha)
  smooth = input("Apply gaussian smoothing? [no]: ") or "no"
  if smooth=="yes":
    sigma = input("Resolution (angstroms) [10]: ") or '10'
    sigma = int(sigma)
    smooth = True
  else:
    sigma = 0
    smooth = False

  # Calculate ratio
  ratio = specRatio(preface, grid, Tphot, Tspot, g, Z, alpha, smooth, sigma)

  # Plot ratio
  x = np.arange(len(ratio))
  if grid=="R":
    ang = np.exp(1e-5*x)*3000
  else:
    ang = 0.1*x + 3000
  plt.plot(ang, ratio)
  plt.xlabel(r"wavelength [$\AA$]")
  plt.ylabel("flux")
  plt.show()

