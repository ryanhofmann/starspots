#!/usr/bin/env python3

import numpy as np
import pymacula

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

# Generate starspots
spots = [pymacula.Spot(lat=lat[i]*np.pi/180., lon=lon[i]*np.pi/180., alpha_max=alpha_max[i], ingress=ingress[i], egress=egress[i], tmax=tmax[i], lifetime=lifetime[i]) for i in range(nspots)]

# Set wavelength range and step value
xmin, xmax, dx = 300, 2500, 10

# Set time range and step value
tmin, tmax, dt = 0, 500, 0.05

