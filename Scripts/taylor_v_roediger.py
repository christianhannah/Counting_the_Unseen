#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 18:48:16 2022

@author: christian
"""

import matplotlib.pyplot as plt
import numpy as np
import mge1d_util as u
from astropy.io import fits
from astropy.io import ascii
import csv
import pdb
from scipy.optimize import curve_fit
from linmix import linmix
import multiprocessing
from matplotlib.colors import LogNorm
import random
import re



# =============================================================================
# # determine the range of g-i colors in our data
# slope_ext = '2x_pixel_scale'
# phys_ext = '_or_10pc_extinction_corr_nsa_ml_w_vi_color'
# gal_file = '../Result_Tables/all_gal_data_'+slope_ext+phys_ext+'.fits'
# hdul = fits.open(gal_file)
# head = hdul[0].data
# dat = hdul[1].data
# hdul.close()
# =============================================================================

# read in fits table with model galaxy parameters
model_gal_filename = '../Result_Tables/final_gal_data.fits'
hdul = fits.open(model_gal_filename)
gal_data = hdul[1].data 
hdul.close()

names_c = gal_data['name']
gicolors_c = gal_data['g_i']   


# =============================================================================
# plt.figure(dpi=500)
# plt.hist(gicolors_c[np.where(gicolors_c != -999)], bins=20)
# plt.xlabel('g-i')
# plt.ylabel('# of galaxies')
# =============================================================================

gi_min = np.min(gicolors_c[np.where(np.isnan(gicolors_c) == False)])
gi_max = np.max(gicolors_c[np.where(np.isnan(gicolors_c) == False)])

g_i = np.arange(gi_min-0.3,2.0,0.05)
log_ML_taylor = -0.68 + 0.70*g_i
log_ML_roed = -0.831 + 0.979*g_i

# =============================================================================
# plt.figure(dpi=500)
# plt.plot(g_i, log_ML_taylor, color='c', label='Taylor+11')
# plt.plot(g_i, log_ML_roed, color='m', label='Roediger+15')
# plt.plot([gi_min, gi_min], [-1.75,2], linestyle='--', color='k')
# plt.plot([gi_max, gi_max], [-1.75,2], linestyle='--', color='k')
# plt.xlabel('g-i')
# plt.ylabel('log(M/L$_i$ [M$_\odot$/L$_{i,\odot}$])')
# plt.legend()
# =============================================================================


fig, ax = plt.subplots(dpi=500)
ax.plot(g_i, 10**log_ML_taylor, color='c', label='Taylor+11')
ax.plot(g_i, 10**log_ML_roed, color='m', label='Roediger+15')
#ax.plot([gi_min, gi_min], [-1.75,2.25], linestyle='--', color='k')
#ax.plot([gi_max, gi_max], [-1.75,2.25], linestyle='--', color='k')
ax.set_xlabel('g-i')
#ax.set_ylabel('log(M/L$_i$ [M$_\odot$/L$_{i,\odot}$])')
ax.set_ylabel('M/L$_i$ [M$_\odot$/L$_{i,\odot}$]')
ax.legend(loc='upper center')
ax.set_xlim(0.25,1.5)
ax.set_ylim(0,5)
ax2 = ax.twinx()
ax2.hist(gicolors_c[np.where(gicolors_c != -999)], bins=40, alpha=0.3, color='0.4')
ax2.set_ylabel('# of galaxies')
ax2.plot([gi_min, gi_min], [0,12.5], linestyle='--', color='k')
ax2.plot([gi_max, gi_max], [0,12.5], linestyle='--', color='k')
plt.show()

