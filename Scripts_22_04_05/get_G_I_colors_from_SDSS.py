#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 18:22:44 2022

@author: christian
"""

from astropy.io import fits
import numpy as np
import pdb
import mge1d_util as u
#from astroML.datasets import fetch_sdss_galaxy_colors



#%%
#dat = fetch_sdss_galaxy_colors()

#%%

# =============================================================================
# read in NASA Sloan Atlas

nsa_file = '/Users/christian/OneDrive/Desktop/TDE Code/nsa_v1_0_1.fits'
hdul = fits.open(nsa_file)
head_nsa = hdul[0].data
data_nsa = hdul[1].data
hdul.close()

# read in our data from our sample

gal_file = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data.fits'
hdul = fits.open(gal_file)
head_gal = hdul[0].data
data_gal = hdul[1].data
hdul.close()

# =============================================================================
# =============================================================================



# =============================================================================
# Perform cross matching via RA and DEC to get colors for our galaxies
RAs = data_gal['RA']
DECs = data_gal['DEC']

G_mags = np.zeros(len(RAs))
I_mags = np.zeros(len(RAs))

no_match_count = 0

for i in range(len(RAs)):
    # do a search for the galaxy based on RA and DEC
    gal_ind_ra = u.find_nearest(data_nsa['RA'],RAs[i])
    gal_ind_dec = u.find_nearest(data_nsa['DEC'],DECs[i])
    
    if gal_ind_ra == gal_ind_dec: # check that the same object was found in both searches
        # convert nanomaggies to magnitudes
        # [u,g,r,i,z,,fd,nd] filters for sersic_nmgy array
        g_nmgy = data_nsa['SERSIC_NMGY'][gal_ind_ra,1]
        G_mags[i] = -2.5*np.log10(g_nmgy) + 22.5
        i_nmgy = data_nsa['SERSIC_NMGY'][gal_ind_ra,3]
        I_mags[i] = -2.5*np.log10(i_nmgy) + 22.5
    else:
        no_match_count += 1
        print('No Match')




# =============================================================================
# =============================================================================

