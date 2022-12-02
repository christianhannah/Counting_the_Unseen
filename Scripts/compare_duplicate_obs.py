#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 09:47:15 2022

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

import warnings
warnings.filterwarnings("ignore")

#%%


# Read in table ffor Pechetti+20 galaxies with RA and Dec
all_names = []
all_RAs = []
all_DECs = []
pechetti_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Data_Sets/Pechetti_20_tables/pechetti_20_ra_dec_NED.csv'
# fields: name	instrument	channel	filter	pa	l_pa	u_pa	ell	l_ell	u_ell	n	l_n	u_n	reff	l_reff	
#         u_reff	ieff	l_ieff	u_ieff	m	l_m	u_m	v	u_v	i	u_i	logm	u_logm
with open(pechetti_filename, 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        #pdb.set_trace()
        if row[0] != 'No.': # marks the start of data row
            all_names.append(re.split('>|<',row[1])[2])
            all_RAs.append(float(row[2]))
            all_DECs.append(float(row[3]))
       
all_names = np.array(all_names)
all_RAs = np.array(all_RAs)
all_DECs = np.array(all_DECs)

# Adjust some of the names to match the names in the Pechetti+20 paper
all_names[np.where(all_names == 'NGC 3115 DW01')] = 'NGC 3115B'
all_names[np.where(all_names == 'ESO 274- G 001')] = 'ESO 274-1'
all_names[np.where(all_names == 'MESSIER 051a')] = 'NGC 5194'
all_names[np.where(all_names == 'MESSIER 051b')] = 'NGC 5195'
all_names[np.where(all_names == 'MESSIER 101')] = 'NGC 5457'
all_names[np.where(all_names == 'MESSIER 063')] = 'NGC 5055'
all_names[np.where(all_names == 'MESSIER 083')] = 'NGC 5236'
all_names[np.where(all_names == 'Circinus Galaxy')] = 'Circinus'

# Remove the spaces in the names
for i in range(len(all_names)):
    s = all_names[i].split()
    if len(s) > 1:
        all_names[i] = s[0]+s[1]


# from Pechetti+20
gal_name = ['Circinus','ESO274-1','IC5052','IC5332','NGC2784','NGC2787',
             'NGC2903','NGC3115','NGC3115B','NGC3184','NGC3274','NGC3344',
             'NGC3593','NGC4242','NGC4460','NGC4517','NGC4592','NGC4600',
             'NGC4605','NGC4941','NGC5055','NGC5068','NGC5194','NGC5195',
             'NGC5236','NGC5238','NGC5457','NGC6503','NGC7713']


# now let's construct RA and DEC arrays for the galaxies in the order they appear above
RAs = np.zeros(len(gal_name))
DECs = np.zeros(len(gal_name))
temp_names = np.zeros(len(gal_name)).astype(str)
for i in range(len(gal_name)):
    ind = np.where(all_names == gal_name[i])[0]
    RAs[i] = all_RAs[ind][0]
    DECs[i] = all_DECs[ind][0]
    #temp_names[i] = all_names[ind][0]



#%%
# put in density/slope values manually as well as errors
log_dens_5pc = np.array([3.84,2.78,3.06,-99,-99,3.93,4.03,-99,3.21,-99,2.35,3.83,-99,
                         2.03,-99,2.61,2.75,3.22,3.65,3.75,3.60,2.68,-99,3.93,3.66,
                         2.04,3.47,3.42,1.83])

log_dens_5pc_err = np.array([.10,.14,.11,-99,-99,.11,.10,-99,.10,-99,.11,.11,-99,.10,
                             -99,.11,.12,.10,.12,.11,.12,.11,-99,.10,.12,.10,.10,.16,
                             .10])
gamma = -1*np.array([.92,2.16,2.06,-99,-99,1.42,1.56,-99,1.86,-99,2.63,1.43,-99,1.80,
                     -99,2.44,2.87,2.09,2.75,1.80,1.87,1.68,-99,1.23,2.48,2.76,1.82,
                     2.75,3.22])
gamma_err = np.array([.12,.07,.04,-99,-99,.03,.08,-99,.04,-99,.04,.11,-99,.12,-99,
                     .03,.02,.04,.02,.12,.06,.26,-99,.24,.04,.03,.02,.15,.01])


#%%
# =============================================================================
# Read in our data with extinction correction nand new M/L method
slope_ext = '2x_pixel_scale'
phys_ext = '_or_10pc_extinction_corr_nsa_ml'
names_c, slopes_c, cen_dens_c, vmags_c, stone_slopes_c, \
    stone_dens_at_5pc_c, NSC_comp_c, lograds_c, logdens_c, dists_c = u.get_our_data(slope_ext,phys_ext,True)

names_c[np.where(names_c == 'PGC 02888')] = 'PGC 028887'
names_c[np.where(names_c == 'PGC 05039')] = 'PGC 050395'

#%%

# get the indices for duplicate obs galaxies
# Lauer+05 and Pechetti+17
i_2778 = np.where(names_c == 'NGC 2778')[0]
i_4458 = np.where(names_c == 'NGC 4458')[0]
i_4660 = np.where(names_c == 'NGC 4660')[0]
i_7457 = np.where(names_c == 'NGC 7457')[0]
# Lauer+05 and Pechetti+17
i_4387 = np.where(names_c == 'NGC 4387')[0]
i_4434 = np.where(names_c == 'NGC 4434')[0]
i_4551 = np.where(names_c == 'NGC 4551')[0]
# Lauer+05 and Pechetti+20
i_3115 = np.where(names_c == 'NGC 3115')[0]

# let's plot the duplicate observations between the samples
# NGC 2778
plt.figure(dpi=500)
plt.plot(lograds_c[i_2778[0],:], logdens_c[i_2778[0],:], color='c', label='Lauer+05')
plt.plot(lograds_c[i_2778[1],:], logdens_c[i_2778[1],:], color='m', label='Pechetti+17')

stone_data = u.get_stone_data('NGC2778')
if stone_data != 0:
    plt.plot(np.log10(stone_data[1]), np.log10(stone_data[2]), color='k', linestyle='--')

plt.legend()
plt.xlim(np.min(lograds_c[i_2778[1],:])-0.2,np.max(lograds_c[i_2778[1],:])+0.2)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.text(0.76,0.58, names_c[i_2778[0]], transform=plt.gcf().transFigure, size=10)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Duplicate_OBS/NGC2778_density_comparison.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# NGC 4458
plt.figure(dpi=500)
plt.plot(lograds_c[i_4458[0],:], logdens_c[i_4458[0],:], color='c', label='Lauer+05')
plt.plot(lograds_c[i_4458[1],:], logdens_c[i_4458[1],:], color='m', label='Pechetti+17')

stone_data = u.get_stone_data('NGC4458')
if stone_data != 0:
    plt.plot(np.log10(stone_data[1]), np.log10(stone_data[2]), color='k', linestyle='--')

plt.legend()
plt.xlim(np.min(lograds_c[i_4458[1],:])-0.2,np.max(lograds_c[i_4458[1],:])+0.2)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.text(0.76,0.58, names_c[i_4458[0]], transform=plt.gcf().transFigure, size=10)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Duplicate_OBS/NGC4458_density_comparison.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# NGC 4660
plt.figure(dpi=500)
plt.plot(lograds_c[i_4660[0],:], logdens_c[i_4660[0],:], color='c', label='Lauer+05')
plt.plot(lograds_c[i_4660[1],:], logdens_c[i_4660[1],:], color='m', label='Pechetti+17')

stone_data = u.get_stone_data('NGC4660')
if stone_data != 0:
    plt.plot(np.log10(stone_data[1]), np.log10(stone_data[2]), color='k', linestyle='--')

plt.legend()
plt.xlim(np.min(lograds_c[i_4660[1],:])-0.2,np.max(lograds_c[i_4660[1],:])+0.2)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.text(0.76,0.58, names_c[i_4660[0]], transform=plt.gcf().transFigure, size=10)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Duplicate_OBS/NGC4660_density_comparison.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# NGC 7457
plt.figure(dpi=500)
plt.plot(lograds_c[i_7457[0],:], logdens_c[i_7457[0],:], color='c', label='Lauer+05')
plt.plot(lograds_c[i_7457[1],:], logdens_c[i_7457[1],:], color='m', label='Pechetti+17')

stone_data = u.get_stone_data('NGC7457')
if stone_data != 0:
    plt.plot(np.log10(stone_data[1]), np.log10(stone_data[2]), color='k', linestyle='--')

plt.legend()
plt.xlim(np.min(lograds_c[i_7457[1],:])-0.2,np.max(lograds_c[i_7457[1],:])+0.2)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.text(0.76,0.58, names_c[i_7457[0]], transform=plt.gcf().transFigure, size=10)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Duplicate_OBS/NGC7457_density_comparison.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)




# NGC 4387
plt.figure(dpi=500)
plt.plot(lograds_c[i_4387[0],:], logdens_c[i_4387[0],:], color='c', label='Lauer+05')
plt.plot(lograds_c[i_4387[1],:], logdens_c[i_4387[1],:], color='m', label='Pechetti+17')

stone_data = u.get_stone_data('NGC4387')
if stone_data != 0:
    plt.plot(np.log10(stone_data[1]), np.log10(stone_data[2]), color='k', linestyle='--')

plt.legend()
plt.xlim(np.min(lograds_c[i_4387[1],:])-0.2,np.max(lograds_c[i_4387[1],:])+0.2)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.text(0.76,0.58, names_c[i_4387[0]], transform=plt.gcf().transFigure, size=10)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Duplicate_OBS/NGC4387_density_comparison.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# NGC 4434
plt.figure(dpi=500)
plt.plot(lograds_c[i_4434[0],:], logdens_c[i_4434[0],:], color='c', label='Lauer+05')
plt.plot(lograds_c[i_4434[1],:], logdens_c[i_4434[1],:], color='m', label='Pechetti+17')

stone_data = u.get_stone_data('NGC4434')
if stone_data != 0:
    plt.plot(np.log10(stone_data[1]), np.log10(stone_data[2]), color='k', linestyle='--')

plt.legend()
plt.xlim(np.min(lograds_c[i_4434[1],:])-0.2,np.max(lograds_c[i_4434[1],:])+0.2)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.text(0.76,0.58, names_c[i_4434[0]], transform=plt.gcf().transFigure, size=10)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Duplicate_OBS/NGC4434_density_comparison.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# NGC 4551
plt.figure(dpi=500)
plt.plot(lograds_c[i_4551[0],:], logdens_c[i_4551[0],:], color='c', label='Lauer+05')
plt.plot(lograds_c[i_4551[1],:], logdens_c[i_4551[1],:], color='m', label='Pechetti+17')

stone_data = u.get_stone_data('NGC4551')
if stone_data != 0:
    plt.plot(np.log10(stone_data[1]), np.log10(stone_data[2]), color='k', linestyle='--')

plt.legend()
plt.xlim(np.min(lograds_c[i_4551[1],:])-0.2,np.max(lograds_c[i_4551[1],:])+0.2)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.text(0.76,0.58, names_c[i_4551[0]], transform=plt.gcf().transFigure, size=10)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Duplicate_OBS/NGC4551_density_comparison.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

#%%

# let's compare M/L's before and after new color relation
gal_file = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data_2x_pixel_scale_or_10pc_extinction_corr_nsa_ml.fits'
hdul = fits.open(gal_file)
head = hdul[0].data
dat = hdul[1].data
hdul.close()

# inds for galaxies included after radius cut of 5 pc
cut_inds = np.where(dat['lograd'][:,0] <= np.log10(5))[0]

mls_before = dat['ml'][cut_inds]
ml_types_before = dat['ml_type'][cut_inds]

gal_file = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data_2x_pixel_scale_or_10pc_extinction_corr_nsa_ml_w_vi_color.fits'
hdul = fits.open(gal_file)
head = hdul[0].data
dat = hdul[1].data
hdul.close()

mls_after = dat['ml'][cut_inds]
ml_types_after = dat['ml_type'][cut_inds]
dists_after = dat['dist'][cut_inds]

inds = np.where(ml_types_after != -99)[0]

#mls_after[np.where(ml_types_after == 5)] *= 4.58

plt.figure(dpi=500)
plt.scatter(mls_before[inds], mls_after[inds], c=ml_types_after[inds], cmap='viridis_r', marker='.')
plt.colorbar()
plt.ylabel('M/L$_{New}$')
plt.xlabel('M/L$_{Original}$')


gal_file = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data_2x_pixel_scale_or_10pc.fits'
hdul = fits.open(gal_file)
head = hdul[0].data
dat = hdul[1].data
hdul.close()
#%%
dists_before = dat['dist'][cut_inds]

plt.figure(dpi=500)
plt.plot(dists_before, dists_after, linestyle='', marker='.')
plt.xlabel('Original Distances')
plt.ylabel('David Distances')





