#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:46:14 2022

Code to explore aspects of JWST NIRCam 

@author: christian
"""

import numpy as np
from astropy.io import fits
import astropy.units as uni
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.table import Table, join, vstack, Column, MaskedColumn
import matplotlib.colors as colors
from astropy.coordinates import Angle
import pdb
import matplotlib.pyplot as plt


arcsec_to_radian = 4.8481e-6

short_res = 0.031 # arcsec
long_res = 0.063 # arcsec
short_psf_fwhm = short_res*2
long_psf_fwhm = long_res*2

dist_array = np.arange(0.1, 50, 0.5)*10**6 # pc


short_pc_res = (short_res*arcsec_to_radian)*dist_array
long_pc_res = (long_res*arcsec_to_radian)*dist_array
short_psf_fwhm_pc = (short_psf_fwhm*arcsec_to_radian)*dist_array
long_psf_fwhm_pc = (long_psf_fwhm*arcsec_to_radian)*dist_array

plt.figure(dpi=500)
plt.plot(dist_array/10**6, short_pc_res, 'b', label='0.6-2.3 micron, pixel-size')
plt.plot(dist_array/10**6, short_psf_fwhm_pc, '--b', label='0.6-2.3 micron, PSF FWHM')
plt.plot(dist_array/10**6, long_pc_res, 'r', label='2.4-5.0 micron, pixel-size')
plt.plot(dist_array/10**6, long_psf_fwhm_pc, '--r', label='2.4-5.0 micron, PSF FWHM')
plt.xlabel('Distance [Mpc]')
plt.ylabel('NIRCam Spatial Resolution [pc]')
plt.legend()