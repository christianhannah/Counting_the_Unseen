#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 13:46:57 2021

@author: christian
"""

from mgefit.mge_fit_1d import mge_fit_1d
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import mge1d_util as u
from scipy import stats

# read in fits table with general galaxy properties
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/gal_properties_Lauer_2005.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMAG','sigma','N','Filt','PS96','Prof','Simbad','NED','_RA','_DE','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data = hdul[1].data 
hdul.close()

names = data['name']
# remove the galaxies whose centers were too obscurred for analysis in Lauer05
names = names[np.where(names != 'NGC 2768')] # no data for this galaxy
names = names[np.where(names != 'NGC 3557')] # no data for this galaxy
names = names[np.where(names != 'NGC 4125')] # no data for this galaxy
names = names[np.where(names != 'NGC 4786')] # no data for this galaxy
names = names[np.where(names != 'NGC 4936')] # no data for this galaxy
names = names[np.where(names != 'NGC 5018')] # no data for this galaxy
names = names[np.where(names != 'NGC 5322')] # no data for this galaxy
names = names[np.where(names != 'NGC 5845')] # no data for this galaxy
names = names[np.where(names != 'NGC 6776')] # no data for this galaxy
names = names[np.where(names != 'NGC 7626')] # no data for this galaxy
names = names[np.where(names != 'IC 3370')] # no data for this galaxy
names = names[np.where(names != 'IC 4296')] # no data for this galaxy

# read in fits table with general galaxy properties
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/gal_properties_Lauer_2005.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMAG','sigma','N','Filt','PS96','Prof','Simbad','NED','_RA','_DE','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data2 = hdul[1].data 
hdul.close()

# read in fits table with general galaxy properties including R_eff from Lauer+2007
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2007_data/gal_properties_Lauer_2007.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMag','logR','r1','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data3 = hdul[1].data 
hdul.close()

david_logmasses = []
nick_logmasses = []
david_dists = []
nick_dists = []
ML_types = []
MLs = []
L_vs = []
M_vs = []

for i in range(len(names)):
    n = names[i].split(' ')
    david_logmasses.append(u.get_galaxy_logmass(n[0]+n[1]))
    david_dists.append(u.get_galaxy_distance(n[0]+n[1]))

    if len(n[1]) == 3:
        gal_name_mod = 'NGC 0'+n[1]
    else:
        gal_name_mod = names[i]
    
    # store the data for specific galaxy
    gal_data = data2[np.where(data2['name']==names[i])]
    gal_data_w_reff = data3[np.where(data3['name']==gal_name_mod)]
    
    # compute M/L ratio for the galaxy
    #M_L = 11 # temporary
    L_v = 10**(0.4*(-gal_data_w_reff['VMAG'][0] + 4.83)) # L_sol_v
    L_vs.append(L_v)
    M_vs.append(gal_data_w_reff['VMAG'][0])
    if len(gal_data_w_reff) > 0:
        logRe = gal_data_w_reff['logR'][0] # log(pc)
    else:
        logRe = 0.0
    disp = gal_data['sigma'][0] # km/s
    G = 4.301e-3 # pc (km/s)^2 M_sol^-1
    if logRe != 0.0:
        Re = 10**logRe
        M_L = (2*disp**2*Re)/(G*L_v)
        ML_type = 0
    else:
        M_L = 4.9*(L_v/10**10)**0.18
        ML_type = 1
        
    MLs.append(M_L)
    ML_types.append(ML_type)
    nick_logmasses.append(np.log10(L_v*M_L))
    nick_dists.append(gal_data_w_reff['Dist'][0])
    


# adjust the log mass arrays to exclude the galaxies that aren't in David's table
d_mass_1 = np.array(david_logmasses)
n_mass_1 = np.array(nick_logmasses)
ML_types_1 = np.array(ML_types)
MLs_1 = np.array(MLs)
L_vs_1 = np.array(L_vs)
M_vs_1 = np.array(M_vs)

keep_inds = np.where(d_mass_1 != 0)
names_new = names[keep_inds]
d_mass = d_mass_1[keep_inds]
n_mass = n_mass_1[keep_inds]
ML_types = ML_types_1[keep_inds]
MLs = MLs_1[keep_inds]
L_vs = L_vs_1[keep_inds]
M_vs = M_vs_1[keep_inds]

#%%
plt.figure(dpi=500)
plt.scatter(n_mass, d_mass, c=ML_types, cmap='viridis')
#plt.plot(n_mass, d_mass, linestyle='', marker='.')
plt.plot([9.4,11.6],[9.4,11.6], linestyle='--', color='r')
plt.xlabel('log(M$_\star$ [M$_\odot$]) from Stone')
plt.ylabel('log(M$_*$ [M$_\odot$]) from David')
plt.colorbar(label='M/L type; 0 = R$_{eff}$, 1 = L$_v$')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_tests/gal_masses_Lauer05_david_v_stone.pdf',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_tests/gal_masses_Lauer05_david_v_stone.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


#%%
# investigate residuals between the two masses
resids = d_mass - n_mass

std = np.std(resids)

inds = []
for i in range(len(resids)):
    if np.abs(resids[i]) > 2*std:
        inds.append(i)
inds = np.array(inds)

outliers = names_new[inds]


mad = stats.median_absolute_deviation(resids)

plt.figure(dpi=500)
plt.hist(resids, bins=60)
plt.plot([-mad, -mad], [0,7], linestyle='--', color='r', label='MAD')
plt.plot([mad, mad], [0,7], linestyle='--', color='r')
plt.plot([-std, -std], [0,7], linestyle='--', color='g', label='STD')
plt.plot([std, std], [0,7], linestyle='--', color='g')
plt.plot([0,0],[0,7], linestyle='-.', color='k')
plt.xlabel('M$_{David}$ - M$_{Nick}$')
plt.ylabel('Counts')
plt.legend()


#clean up distance arrays to remove the gals where david doesn't have distances
nick_dists_1 = np.array(nick_dists)
david_dists_1 = np.array(david_dists)
keeps = np.where(david_dists_1 != 0)
nick_dists = nick_dists_1[keeps]
david_dists = david_dists_1[keeps]

# single out data points that are outliers (NGC 1023, 3945, 4406)
n1023 = nick_dists[np.where(names == 'NGC 1023')[0]]
n3945 = nick_dists[np.where(names == 'NGC 3945')[0]]
n4406 = nick_dists[np.where(names == 'NGC 4406')[0]]
d1023 = david_dists[np.where(names == 'NGC 1023')[0]]
d3945 = david_dists[np.where(names == 'NGC 3945')[0]]
d4406 = david_dists[np.where(names == 'NGC 4406')[0]]


plt.figure(dpi=500)
plt.plot(nick_dists, david_dists, linestyle='', marker='.', color='b')
plt.plot([10,54],[10,54], linestyle='-.', alpha=0.7, color='k')
plt.plot(n1023,d1023,color='r',marker='*')
plt.plot(n3945,d3945,color='r',marker='*')
plt.plot(n4406,d4406,color='r',marker='*')
plt.xlabel('Lauer+07 Distances [Mpc]')
plt.ylabel('David Distances [Mpc]')

#%%
# read in fits table with galaxy masses from David
filename = '/Users/christian/OneDrive/Desktop/TDE Code/Galaxy Masses/50mpcsample.fits'
# fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
hdul = fits.open(filename)  # open a FITS file
data = hdul[1].data 
hdul.close()

