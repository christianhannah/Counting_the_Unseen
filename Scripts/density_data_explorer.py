#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 13:42:01 2022

@author: christian
"""

import matplotlib.pyplot as plt
import numpy as np
import mge1d_util as u
from astropy.io import fits
import warnings
warnings.filterwarnings("ignore")

#%%

slope_ext = '2x_pixel_scale'
phys_ext = '_or_10pc_extinction_corr_nsa_ml_w_vi_color'

gal_file = '../Result_Tables/all_gal_data_'+slope_ext+phys_ext+'.fits'
hdul = fits.open(gal_file)
head = hdul[0].data
dat = hdul[1].data
hdul.close()

names = dat['name']
vmags = dat['vmag']
dists = dat['dist'] # Mpc
MLs = dat['ml']
ML_types = dat['ml_type']
slopes = dat['slope']
cen_dens = dat['dens_at_5pc'] # M_sol/pc^3
lograds = dat['lograd'] # log(pc)
logdens = dat['logdens'] # log(M_sol/pc^3)
stone_slopes = dat['stone_slope']
stone_dens_at_5pc = dat['stone_dens_at_5pc']
NSC_comp = dat['NSC_comp']
all_MLs = dat['all_mls']
all_ML_types = dat['all_ml_types']
dist_flags = dat['dist_flag']
SBs = dat['SBs']
filts = dat['filt']


# add in some code here to read in Pechetti+20 density data once I have it
#%%
# read in fits table with model galaxy parameters
model_gal_filename = '../Result_Tables/final_gal_data.fits'
hdul = fits.open(model_gal_filename)
gal_data = hdul[1].data 
hdul.close()



#%%


# let's isolate the SB profiles of the nucleated galaxies in F814W
i_inds = np.where(filts == 'F814W')
nuc_inds = np.where(NSC_comp == 1)
final_inds = np.intersect1d(i_inds,nuc_inds)

names_nuc = names[final_inds]
SBs_nuc = SBs[final_inds]
lograds_nuc = lograds[final_inds]
dists_nuc = dists[final_inds]


# first, let's convolve the SB profile with the PSF of NIRCam 
import astropy.convolution as apcon
sigma2fwhm=2*np.sqrt(2*np.log(2))
rad_spacing = lograds_nuc[:,1]-lograds_nuc[:,0]
psf_sig_short = (2*0.031*4.8481*dists_nuc/rad_spacing)/sigma2fwhm
pix_sig_short = 0.031/sigma2fwhm
psf_sig_long = 2*0.063/sigma2fwhm
pix_sig_long = 0.063/sigma2fwhm
SBs_conv = np.zeros_like(SBs_nuc)
for i in range(len(names_nuc)):
    SBs_temp = np.concatenate((np.flip(SBs_nuc[i,:]),SBs_nuc[i,1:]))
    SBs_c = apcon.convolve(SBs_temp, apcon.Gaussian1DKernel(psf_sig_short[i]))
    SBs_conv[i,:] = SBs_c[99:]

cen_SBs = SBs_conv[:,0]

# convert the SBs back to mag/arcsec^2
M_sol = 4.52 # ABmag F814W
cen_SBs_1 = -2.5*np.log10(cen_SBs)+4.52+21.572



