#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:11:38 2021

@author: christian
"""

import glob
import csv
import mge_fit_Lauer95_data as m95
import mge_fit_Lauer05_data as m05
import mge_fit_Laine03_data as m03
import mge_fit_Lauer05_data_minus_inner_radius as mr
import density_profiles_pechetti_17 as mp17
import mge1d_util as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import pdb
from tqdm import tqdm
import sys 
import os
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import ascii
from astropy.stats import median_absolute_deviation as mad

# specify largest radii to fit in pc
max_rad = 100000 # pc

# specify physical and minimum angular scale multiplier for the slope fits
phys_scale = 10
phys_ext = '_or_{}pc'.format(phys_scale)
#phys_ext = '10pc'
scale_mult = 2
slope_ext = '{}x_pixel_scale'.format(scale_mult)
#slope_ext = 'pixel_scale'
#slope_ext = ''
#%%
# =============================================================================
# ======================= Lauer et al. (2005) =================================
# =============================================================================

# read in fits table with general galaxy properties
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/gal_properties_Lauer_2005.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMAG','sigma','N','Filt','PS96','Prof','Simbad','NED','_RA','_DE','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data_05 = hdul[1].data 
hdul.close()
#%%
old_stdout = sys.stdout
#%%
names_05 = data_05['name']
# remove the galaxies whose centers were too obscurred for analysis in Lauer05
names_05 = names_05[np.where(names_05 != 'NGC 2768')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 3557')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 4125')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 4786')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 4936')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 5018')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 5322')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 5845')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 6776')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 7626')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'IC 3370')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'IC 4296')] # no data for this galaxy
names_05 = names_05[np.where(names_05 != 'NGC 3706')] # weird sampling from Lauer+02
names_05 = names_05[np.where(names_05 != 'NGC 4406')] # weird sampling from Lauer+02
names_05 = names_05[np.where(names_05 != 'NGC 6876')] # weird sampling from Lauer+02

#%%
dens_profiles_05 = np.zeros((100,len(names_05)))
radii_05 = np.zeros((100,len(names_05)))
lums_05 = np.zeros(len(names_05))
ML_types_05 = np.zeros(len(names_05))
MLs_05 = np.zeros(len(names_05))
dists_05 = np.zeros(len(names_05))
ras_05 = np.zeros(len(names_05))
decs_05 = np.zeros(len(names_05))
SBs_true_05 = []
rad_true_05 = []
SBs_interp_05 = np.zeros((100,len(names_05)))
rad_interp_05 = np.zeros((100,len(names_05)))
g_heights_05 = []
g_sigs_05 = []
SBs_interp_converted_05 = np.zeros((100,len(names_05)))
SBs_true_converted_05 = []
for i in tqdm(range(len(names_05)), position=0, leave=True):
    s = names_05[i].split(' ')
    sys.stdout = open(os.devnull, "w")
    dens_profiles_05[:,i], radii_05[:,i], lums_05[i], MLs_05[i], \
        ML_types_05[i], dists_05[i], ras_05[i], decs_05[i],\
        true_SB, true_rad, SBs_interp_05[:,i], rad_interp_05[:,i],\
        heights, sigs, SBs_interp_converted_05[:,i],\
        true_SB_converted = m05.get_density_profile(s[0], s[1],max_rad)
    SBs_true_05.append(true_SB)
    rad_true_05.append(true_rad)
    g_heights_05.append(heights)
    g_sigs_05.append(sigs)
    SBs_true_converted_05.append(true_SB_converted)
    #dens_profiles[:,i], radii[:,i], lums[i] = mr.get_density_profile(s[0], s[1])
    sys.stdout = old_stdout
    #pdb.set_trace()

# =============================================================================
# =============================================================================
# =============================================================================
#%%
# =============================================================================
# ======================= Lauer et al. (1995)==================================
# =============================================================================

SB_files = glob.glob('/Users/christian/OneDrive/Desktop/TDE Code/Lauer1995_data/*.comp.txt')
names_95 = []
for i in range(len(SB_files)):
    names_95.append(SB_files[i].split('/')[-1].split('.')[0][0:5]) 
    
# remove the duplicates that are in laine03 (NGC 2832, NGC 4889, A2052-M1)
names_95 = np.array(names_95)
#names_95 = names_95[np.where(names_95 != 'n2832')] # same coverage as Laine+03
#names_95 = names_95[np.where(names_95 != 'n4889')] # same coverage as Laine+03
names_95 = names_95[np.where(names_95 != 'a1831')] # not deep coverage 
#names_95 = names_95[np.where(names_95 != 'a2052')] # deeper coverage than Laine+03


# NOTE: The file for ngc 4486 (M87) conatined 2 data sets so it has been removed
#%%
#pdb.set_trace()
dens_profiles_95 = np.zeros((100,len(names_95)))
radii_95 = np.zeros((100,len(names_95)))
lums_95 = np.zeros(len(names_95))
ML_types_95 = np.zeros(len(names_95))
MLs_95 = np.zeros(len(names_95))
dists_95 = np.zeros(len(names_95))
ras_95 = np.zeros(len(names_95))
decs_95 = np.zeros(len(names_95))
SBs_true_95 = []
rad_true_95 = []
SBs_interp_95 = np.zeros((100,len(names_95)))
rad_interp_95 = np.zeros((100,len(names_95)))
g_heights_95 = []
g_sigs_95 = []
SBs_interp_converted_95 = np.zeros((100,len(names_95)))
SBs_true_converted_95 = []
for i in tqdm(range(len(names_95)), position=0, leave=True):
    sys.stdout = open(os.devnull, "w")
    dens_profiles_95[:,i], radii_95[:,i], lums_95[i], \
        MLs_95[i], ML_types_95[i], dists_95[i], ras_95[i], \
        decs_95[i], true_SB, true_rad, SBs_interp_95[:,i], rad_interp_95[:,i],\
        heights, sigs, SBs_interp_converted_95[:,i],\
        true_SB_converted = m95.get_density_profile(names_95[i],max_rad)
    SBs_true_95.append(true_SB)
    rad_true_95.append(true_rad)
    g_heights_95.append(heights)
    g_sigs_95.append(sigs)
    SBs_true_converted_95.append(true_SB_converted)
    sys.stdout = old_stdout
    #pdb.set_trace()


# =============================================================================
# =============================================================================
# =============================================================================


#%%
# =============================================================================
# ====================== Pechetti et al. (2017) ===============================
# =============================================================================

# =============================================================================
# # read names for pechetti galaxies
# mge_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Pechetti17_data/pechetti_2017_data.fit'
# hdul = fits.open(mge_filename)
# data1 = hdul[1].data  
# names_17 = np.unique(data1['Name'])
# 
# dens_profiles_17 = np.zeros((100,len(names_17)))
# radii_17 = np.zeros((100,len(names_17)))
# lums_17 = np.zeros(len(names_17))
# ML_types_17 = np.zeros(len(names_17))
# MLs_17 = np.zeros(len(names_17))
# dists_17 = np.zeros(len(names_17))
# ras_17 = np.zeros(len(names_17))
# decs_17 = np.zeros(len(names_17))
# filts_17 = np.zeros(len(names_17)).astype(str)
# for i in tqdm(range(len(names_17)), position=0, leave=True):
#     sys.stdout = open(os.devnull, "w")
#     dens_profiles_17[:,i], radii_17[:,i], lums_17[i], \
#         MLs_17[i], ML_types_17[i], dists_17[i], ras_17[i], \
#         decs_17[i], filts_17[i] = mp17.get_density_profile(names_17[i])
#     sys.stdout = old_stdout
#     #pdb.set_trace()
# =============================================================================

# =============================================================================
# =============================================================================
# =============================================================================


#%%
# =============================================================================
# ====================== Laine et al. (2003) ==================================
# =============================================================================

# =============================================================================
# #SB_files = glob.glob('/Users/christian/OneDrive/Desktop/TDE Code/Laine2003_BCG_data/bcg_profs/*.txt')
# names_03 = np.array(['A3744-M1', 'A3395-M1', 'A3747-M1', 'IC 5353', 'A3554-M1',
#                      'NGC 7014', 'A3144-M1', 'NGC 1500', 'A0261-M1', 'A1631-M1',
#                      'A3556-M1', 'IC 1565', 'IC 0664', 'IC 1733', 'NGC 4889',
#                      'A0548-M1', 'NGC 6086', 'A3532-M1', 'NGC 4696', 'NGC 7578B',
#                      'NGC 7647', 'A1308-M1', 'A0189-M1', 'A0160-M1', 'A0376-M1',
#                      'A3736-M1', 'IC 4931', 'IC 1695', 'NGC 6173', 'NGC 0910',
#                      'A2593-M1', 'A0147-M1', 'A0999-M1', 'A0634-M1', 'IC 0115',
#                      'A0397-M1', 'A2147-M1', 'IC 0712', 'NGC 2832', 'NGC 0545',
#                      'A0419-M1', 'A3716-M1', 'A3528-M1', 'A3677-M1', 'A0168-M1',
#                      'IC 1633', 'A0912-M1', 'A3571-M1', 'A3559-M1', 'A3558-M1',
#                      'A3570-M1', 'A3564-M1', 'A0496-M1', 'A2040-M1',
#                      'A0119-M1', 'A3376-M1', 'A1983-M1', 'A0295-M1',
#                      'A2247-M1', 'A3367-M1', 'IC 2738', 'IC 0613', 'A3562-M1',
#                      'NGC 3551', 'A0533-M1', 'A1836-M1'])
#                      #, 'A2052-M1'])
# 
# # duplicate obs in Lauer05 removed -> NGC 3842, IC 4329 (same coverage, just newer observation)
# =============================================================================
#%%
# read in fits table with general galaxy properties including R_eff from Lauer+2007
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2007_data/gal_properties_Lauer_2007.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMag','logR','r1','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data_07 = hdul[1].data 
hdul.close()
# =============================================================================
# #%%
# # remove galaxies from Laine+03 sample that aren't in Lauer+07
# in_lauer = []
# for i in range(len(names_03)):
#     if len(data_07[np.where(data_07['name']==names_03[i])]) != 0:
#         in_lauer.append(names_03[i])
# names_03 = np.array(in_lauer)
# #%%
# 
# dens_profiles_03 = np.zeros((100,len(names_03)))
# radii_03 = np.zeros((100,len(names_03)))
# lums_03 = np.zeros(len(names_03))
# ML_types_03 = np.zeros(len(names_03))
# MLs_03 = np.zeros(len(names_03))
# dists_03 = np.zeros(len(names_03))
# ras_03 = np.zeros(len(names_03))
# decs_03 = np.zeros(len(names_03))
# SBs_true_03 = []
# rad_true_03 = []
# SBs_interp_03 = np.zeros((100,len(names_03)))
# rad_interp_03 = np.zeros((100,len(names_03)))
# g_heights_03 = []
# g_sigs_03 = []
# SBs_interp_converted_03 = np.zeros((100,len(names_03)))
# for i in tqdm(range(len(names_03)), position=0, leave=True):
#     sys.stdout = open(os.devnull, "w")
#     dens_profiles_03[:,i], radii_03[:,i], lums_03[i], MLs_03[i], \
#         ML_types_03[i], dists_03[i], ras_03[i], decs_03[i], true_SB, true_rad,\
#         SBs_interp_03[:,i], rad_interp_03[:,i],heights,\
#         sigs, SBs_interp_converted_03[:,i] = m03.get_density_profile(names_03[i],max_rad)
#     SBs_true_03.append(true_SB)
#     rad_true_03.append(true_rad)
#     g_heights_03.append(heights)
#     g_sigs_03.append(sigs)
#     sys.stdout = old_stdout
#     #pdb.set_trace()
# =============================================================================

# =============================================================================
# =============================================================================
# =============================================================================

#%%
# =============================================================================
# plot all density profiles together by galaxy v-band magnitude

ls = []
all_lograds = []
all_logdens = []
# add lauer05 profiles to list
for i in range(len(names_05)):
    a = np.zeros((len(radii_05[:,2]),2))
    a[:,0] = np.log10(radii_05[:,i])
    a[:,1] = np.log10(dens_profiles_05[:,i])
    all_lograds.append(np.log10(radii_05[:,i]))
    all_logdens.append(np.log10(dens_profiles_05[:,i]))
    ls.append(a)
# add lauer95 profiles to list
for i in range(len(names_95)):
    a = np.zeros((len(radii_95[:,2]),2))
    a[:,0] = np.log10(radii_95[:,i])
    a[:,1] = np.log10(dens_profiles_95[:,i])
    all_lograds.append(np.log10(radii_95[:,i]))
    all_logdens.append(np.log10(dens_profiles_95[:,i]))
    ls.append(a)
# =============================================================================
# # add pechetti17 profiles to list
# for i in range(len(names_17)):
#     a = np.zeros((len(radii_17[:,2]),2))
#     a[:,0] = np.log10(radii_17[:,i])
#     a[:,1] = np.log10(dens_profiles_17[:,i])
#     all_lograds.append(np.log10(radii_17[:,i]))
#     all_logdens.append(np.log10(dens_profiles_17[:,i]))
#     ls.append(a)
# =============================================================================
# =============================================================================
# # add laine03 profiles to list
# for i in range(len(names_03)):
#     a = np.zeros((len(radii_03[:,2]),2))
#     a[:,0] = np.log10(radii_03[:,i])
#     a[:,1] = np.log10(dens_profiles_03[:,i])
#     all_lograds.append(np.log10(radii_03[:,i]))
#     all_logdens.append(np.log10(dens_profiles_03[:,i]))
#     ls.append(a)
# =============================================================================

# concatenate all the v-band magnitudes for colorbar
lumys = np.concatenate((lums_05,lums_95))#,lums_17))#,lums_03))

fig, ax = plt.subplots(dpi=500)
lines = LineCollection(ls, array=lumys, cmap='viridis', alpha=0.6)
ax.add_collection(lines)
fig.colorbar(lines, label='M$_{v}$ [mag]')
ax.set_xlabel('log(Radius [pc])')
ax.set_ylabel('log(Density [M$_{\odot}$/pc$^3$])')
ax.autoscale()
plt.show()
plt.close()
# =============================================================================

# =============================================================================
#%%

# get the inner density and powerlaw slope of each density profile

# measure slopes within the HST pixel scale so need to compute this
# conversion factor between radians and arcsec
arcsec_to_radian = 4.8481e-6
hst_ps_05 = (0.0456*arcsec_to_radian)*dists_05
hst_ps_95 = (0.0456*arcsec_to_radian)*dists_95
#hst_ps_17 = (0.05*arcsec_to_radian)*dists_17
all_hst_ps = np.concatenate((hst_ps_05,hst_ps_95))#,hst_ps_17))


central_dens = []
slopes = []
for i in range(len(ls)):
    s, cd, blah, blah = u.get_density_and_slope(all_lograds[i], all_logdens[i], all_hst_ps[i],
                                    phys_scale, scale_mult)
    slopes.append(s)
    central_dens.append(cd)
slopes = np.array(slopes)
central_dens = np.array(central_dens)

# break up these arrays by samle for pdf plotting
slopes_05 = slopes[0:len(names_05)]
slopes_95 = slopes[len(names_05):len(names_05)+len(names_95)]
slopes_17 = slopes[len(names_05)+len(names_95):]
#slopes_03 = slopes[len(names_95):]
central_dens_05 = central_dens[0:len(names_05)]
central_dens_95 = central_dens[len(names_05):len(names_05)+len(names_95)]
central_dens_17 = central_dens[len(names_05)+len(names_95):]
#central_dens_03 = central_dens[len(names_95):]

# =============================================================================
#%% 
# =============================================================================
# loop through 05 sample galaxies where density profiles were just computed to 
# check for the corresponding profile from Stone&Metzger2016
resids = []
gal_logmasses = []
stone_radii_05 = []
stone_densities_05 = []
for i in range(len(names_05)):
    n = names_05[i].split()
    if u.check_for_galaxy_in_stone(n[0]+n[1]) == False:
        print(n[0]+' '+n[1]+' not in Stone&Metzger(2016).')
        
        stone_radii_05.append(np.array([0]))
        stone_densities_05.append(np.array([0]))
        
    else:
        # read data from Stone & Metzger (2016) who deprojected densities using
        # nuker law SB profiles (essentially double power-laws) and assuming 
        # spherical symmetry
        stone_data =  u.get_stone_data(n[0]+n[1])   
        stone_rad = stone_data[1]
        stone_dens = stone_data[2]
        
        stone_radii_05.append(stone_rad)
        stone_densities_05.append(stone_dens)
        
        new_stone_dens = np.interp(radii_05[:,i], stone_rad, stone_dens)
        b = np.zeros((len(radii_05[:,i]),2))
        b[:,0] = np.log10(radii_05[:,i])
        b[:,1] = np.log10(dens_profiles_05[:,i]) - np.log10(new_stone_dens)
        resids.append(b)
        
        # store log mass for each galaxy
        gal_logmasses.append(u.get_galaxy_logmass(n[0]+n[1]))
      

#%%
# plot residuals between stone profiles and ours for 05 sample colored by ML_type
fig, ax = plt.subplots(dpi=500)
lines = LineCollection(resids, array=ML_types_05, cmap='viridis', alpha=0.6)
ax.add_collection(lines)
fig.colorbar(lines, label='M/L type; 0 = R$_{eff}$, 1 = L$_v$')
ax.set_xlabel('log(Radius [pc])')
ax.set_ylabel('log(Density Residuals [M$_{\odot}$/pc$^3$])')
ax.autoscale()
plt.show()
plt.close()
# =============================================================================


#%%
# =============================================================================
# plot all density profiles with stone profiles where available

# =============================================================================
# for i in range(len(names_05)):
#     plt.figure(dpi=500)
#     plt.plot(all_lograds[i],all_logdens[i],label=names_05[i])
#     if len(stone_radii_05[i]) != 1:
#         plt.plot(np.log10(stone_radii_05[i]),np.log10(stone_densities_05[i]),label='Stone&Metzger')
#     plt.xlabel('log(Radius [pc]')
#     plt.ylabel('log(Density Residuals [M$_{\odot}$/pc$^3$])')
#     plt.legend()
# 
# =============================================================================
# =============================================================================
#%%


# =============================================================================
# plot residuals between stone profiles and ours for 05 sample colored by galaxy
# mass where available

gal_logmasses = np.array(gal_logmasses)
mass_inds = np.where(gal_logmasses != 0)[0]
gal_logmasses_1 = gal_logmasses[mass_inds]
resids_1 = []
for i in range(len(resids)):
    if i in mass_inds:
        resids_1.append(resids[i])
    
fig, ax = plt.subplots(dpi=500)
lines = LineCollection(resids_1, array=gal_logmasses_1, cmap='viridis', alpha=0.6)
ax.add_collection(lines)
fig.colorbar(lines, label='log(Galaxy Mass [M$_\odot$])')
ax.set_xlabel('log(Radius [pc])')
ax.set_ylabel('log(Density Residuals [M$_{\odot}$/pc$^3$])')
ax.autoscale()
plt.show()
plt.close()
# =============================================================================

# =============================================================================
# let's shift residuals to zero at the norm radius
norm_resids = []
norm_rad = np.log10(500) # log(pc)
for i in range(len(resids)):
    x = np.zeros(len(resids[i]))
    y = np.zeros(len(resids[i]))
    for j in range(len(resids[i])):
        x[j] = resids[i][j][0]
        y[j] = resids[i][j][1]        
    norm_ind = u.find_nearest(x, norm_rad)
    y_norm = y-y[norm_ind]
    b = np.zeros((len(x),2))
    b[:,0] = x
    b[:,1] = y_norm
    norm_resids.append(b)

fig, ax = plt.subplots(dpi=500)
lines = LineCollection(norm_resids, array=ML_types_05, cmap='viridis', alpha=0.6)
ax.add_collection(lines)
fig.colorbar(lines, label='M/L type; 0 = R$_{eff}$, 1 = L$_v$')
ax.set_xlabel('log(Radius [pc])')
ax.set_ylabel('log(Density Residuals [M$_{\odot}$/pc$^3$])')
ax.autoscale()
plt.show()
plt.close()
# =============================================================================

#%%
# =============================================================================
# let's shift residuals with gal_logmasses to zero at the norm radius
norm_resids_1 = []
norm_rad = np.log10(300) # log(pc)
for i in range(len(resids_1)):
    x = np.zeros(len(resids_1[i]))
    y = np.zeros(len(resids_1[i]))
    for j in range(len(resids_1[i])):
        x[j] = resids_1[i][j][0]
        y[j] = resids_1[i][j][1]        
    norm_ind = u.find_nearest(x, norm_rad)
    y_norm = y-y[norm_ind]
    b = np.zeros((len(x),2))
    b[:,0] = x
    b[:,1] = y_norm
    norm_resids_1.append(b)


fig, ax = plt.subplots(dpi=500)
lines = LineCollection(norm_resids_1, array=gal_logmasses_1, cmap='viridis', alpha=0.6)
ax.add_collection(lines)
fig.colorbar(lines, label='log(Galaxy Mass [M$_\odot$])')
ax.set_xlabel('log(Radius [pc])')
ax.set_ylabel('log(Density Residuals [M$_{\odot}$/pc$^3$])')
ax.autoscale()
plt.show()
plt.close()
# =============================================================================

#%%
# =============================================================================
# compare lists of galaxy names to see what we have and what we need

# convert names_95 to the Lauer2007 naming convention
names_95_c = []
for i in range(len(names_95)):
    if names_95[i][0] == 'n':
        names_95_c.append('NGC '+names_95[i][1:])
    elif names_95[i][0] == 'v':
        names_95_c.append('VCC '+names_95[i][1:])
    elif names_95[i][0] == 'a':
        names_95_c.append('A'+names_95[i][1:]+'-M1')
names_95_c = np.array(names_95_c)
  #%%  
# correct the first few names from lauer05 for the same reason
names_05_c = names_05
for i in range(5):
    s = names_05_c[i].split()
    names_05_c[i] = 'NGC 0'+s[1]
    
#%%    
all_names = np.concatenate((names_05_c,names_95_c))#,names_17))#,names_03))

# remove any weird spaces at the end of the string names
for i in range(len(all_names)):
     a = all_names[i][4:]
     all_names[i] = all_names[i][0:4]+a.replace(" ","")

# =============================================================================
# # find duplicates between lauer95 and laine03
# m1 = np.zeros_like(all_names, dtype=bool)
# m1[np.unique(all_names, return_index=True)[1]] = True
# print(all_names[~m1])
# 
# #%%
# 
# all_lauer07_names = data_07['name']
# all_names_c = all_names#[m1]
# missing_gals = []
# gals_we_have = []
# 
# # =============================================================================
# # # remove any weird spaces at the end of the string names
# # for i in range(len(all_names_c)):
# #      a = all_names_c[i][5:]
# #      all_names_c[i] = all_names_c[i][0:5]+a.replace(" ","")
# # =============================================================================
#      
# m2 = np.zeros_like(all_names_c, dtype=bool)
# m2[np.unique(all_names_c, return_index=True)[1]] = True
# print(all_names_c[~m2])    
# 
# #%%
# 
# for i in range(len(all_lauer07_names)):
#     have = False
#     for j in range(len(all_names_c)):
#         if all_lauer07_names[i] == all_names_c[j]:
#             have = True
#             gals_we_have.append(all_lauer07_names[i])
#             break
#     if have == False:
#         #pdb.set_trace()
#         missing_gals.append(all_lauer07_names[i])
# 
# missing_gals = np.array(missing_gals)
# gals_we_have = np.array(gals_we_have)
# 
# 
# extra_gals = []
# for i in range(len(all_names_c)):
#     h = False
#     for j in range(len(gals_we_have)):
#         if all_names_c[i] == gals_we_have[j]:
#             h = True
#             break
#     if h == False:
#         extra_gals.append(all_names_c[i])
# extra_gals = np.array(extra_gals)
# 
# =============================================================================
#%%

# let's find the density profiles that are also in Stone&Metzger 2016 and compute the slopes
all_names_nospace = []
stone_slopes = np.zeros(len(all_names))
stone_dens_at_5pc = np.zeros(len(all_names))
stone_radii = []
stone_densities = []
for i in range(len(all_names)):
    all_names_nospace.append(u.format_gal_name(all_names[i]))
    stone_data =  u.get_stone_data(all_names_nospace[i])   
    if stone_data != 0:
        stone_rad = np.log10(stone_data[1])
        stone_dens = np.log10(stone_data[2])
        stone_slopes[i], stone_dens_at_5pc[i] = u.get_density_and_slope_simple(stone_rad, stone_dens)
        stone_radii.append(stone_data[1])
        stone_densities.append(stone_data[2])
    else:
        stone_slopes[i] = -99
        stone_dens_at_5pc[i] = -99
        stone_radii.append(np.array([0]))
        stone_densities.append(np.array([0]))
#%%
# let's mark the galaxies with noted NSC components
data05_3 = ascii.read(table='/Users/christian/OneDrive/Desktop/TDE Code/Lauer05_Table3.txt', format='aastex') 
NSC_comp = np.zeros_like(all_names)
for i in range(len(all_names)):
    found = False
    for j in range(len(data05_3['Gal Name'])):
        if data05_3['Gal Name'][j] == all_names[i]:
            found = True
            if data05_3['Nuc'][j] == '$\\times$':
                NSC_comp[i] = 1
            else:
                NSC_comp[i] = 0
    if not found:
        NSC_comp[i] = -1


#%%
# =============================================================================
# store the relevant data to a fits table
c1 = fits.Column(name='name', array=all_names, format='10A')
c2 = fits.Column(name='vmag', array=lumys, format='D', unit='mag')

# concatenate distances, M/L ratios, and M/L types
all_dists = np.concatenate((dists_05,dists_95))#,dists_17))#,dists_03))
all_MLs = np.concatenate((MLs_05,MLs_95))#,MLs_17))#,MLs_03))
all_ML_types = np.concatenate((ML_types_05,ML_types_95))#, ML_types_17))#,ML_types_03))
all_ras = np.concatenate((ras_05,ras_95))#, ras_17))#,ras_03))
all_decs = np.concatenate((decs_05,decs_95))#,decs_17))#,decs_03))


c3 = fits.Column(name='dist', array=all_dists/10**6, format='D', unit='Mpc')
c4 = fits.Column(name='ml', array=all_MLs, format='D')
c5 = fits.Column(name='ml_type', array=all_ML_types, format='I')
c6 = fits.Column(name='slope', array=slopes, format='D')
c7 = fits.Column(name='dens_at_5pc', array=central_dens, format='D', unit='M_sol/pc^3')
c8 = fits.Column(name='lograd', array=all_lograds, format='100D', unit='log(pc)')
c9 = fits.Column(name='logdens', array=all_logdens, format='100D', unit='log(M_sol/pc^3)')
c10 = fits.Column(name='RA', array=all_ras, format='D', unit='deg')
c11 = fits.Column(name='DEC', array=all_decs, format='D', unit='deg')
c12 = fits.Column(name='stone_slope', array=stone_slopes, format='D')
c13 = fits.Column(name='stone_dens_at_5pc', array=stone_dens_at_5pc, format='D', unit='M_sol/pc^3')
c14 = fits.Column(name='NSC_comp', array=NSC_comp, format='I')


t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14])
t.writeto('/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data_'+slope_ext+phys_ext+'.fits',clobber=True)

# =============================================================================

#%%
# =============================================================================
# generate a plot of 3rd original data point vs 0.05" at each galaxy distance

# conversion factor between radians and arcsec
arcsec_to_radian = 4.8481e-6

hst_scales = []
third_data_pt = []
second_data_pt = []
all_names_5pc = []
names_05_5pc = []
names_95_5pc = []
for i in range(len(names_05)):
    # apply the 5pc data cut
    if all_lograds[i][0] <= np.log10(5):
        third_data_pt.append((rad_true_05[i][2]*arcsec_to_radian)*dists_05[i])
        second_data_pt.append((rad_true_05[i][1]*arcsec_to_radian)*dists_05[i])
        hst_scales.append((0.0456*arcsec_to_radian)*dists_05[i])
        all_names_5pc.append(names_05[i])
        names_05_5pc.append(names_05[i])
for i in range(len(names_95)):
    # apply the 5pc data cut
    if all_lograds[len(names_05)+i][0] <= np.log10(5):
        third_data_pt.append((rad_true_95[i][2]*arcsec_to_radian)*dists_95[i])
        second_data_pt.append((rad_true_95[i][1]*arcsec_to_radian)*dists_95[i])
        hst_scales.append((0.0456*arcsec_to_radian)*dists_95[i])
        all_names_5pc.append(names_95[i])
        names_95_5pc.append(names_95[i])

hst_scales = np.array(hst_scales)
third_data_pt = np.array(third_data_pt)
second_data_pt = np.array(second_data_pt)
all_names_5pc = np.array(all_names_5pc)

plt.figure(dpi=500)
plt.plot(hst_scales[0:len(names_05_5pc)], third_data_pt[0:len(names_05_5pc)], linestyle='', 
         marker='.', color='c',label='Lauer+05')
plt.plot(hst_scales[len(names_05_5pc):], third_data_pt[len(names_05_5pc):], linestyle='', 
         marker='.', color='b',label='Lauer+95')
plt.plot([0,14],[0,14],linestyle='--',color='r')
plt.xlabel('HST Pixel Scale [pc]')
plt.ylabel('Radial Coverage of 3rd Innermost Data Point [pc]')
plt.legend()

plt.figure(dpi=500)
plt.plot(hst_scales[0:len(names_05_5pc)], second_data_pt[0:len(names_05_5pc)], linestyle='',
         marker='.', color='c',label='Lauer+05')
plt.plot(hst_scales[len(names_05_5pc):], second_data_pt[len(names_05_5pc):], linestyle='', 
         marker='.', color='b',label='Lauer+95')
plt.plot([0,14],[0,14],linestyle='--',color='r')
plt.xlabel('HST Pixel Scale [pc]')
plt.ylabel('Radial Coverage of 2nd Innermost Data Point [pc]')
plt.legend()

# =============================================================================

#%%
# =============================================================================
# plot distribution of innermost radii
min_rads = np.zeros(len(names_05)+len(names_95))
min_2_rads = np.zeros(len(names_05)+len(names_95))
min_3_rads = np.zeros(len(names_05)+len(names_95))
for i in range(len(names_05)):
    min_rads[i] = rad_true_05[i][0]
    min_2_rads[i] = rad_true_05[i][1]
    min_3_rads[i] = rad_true_05[i][2]
for i in range(len(names_95)):
    min_rads[len(names_05)+i] = rad_true_95[i][0]
    min_2_rads[len(names_05)+i] = rad_true_95[i][1]
    min_3_rads[len(names_05)+i] = rad_true_95[i][2]


a = np.unique(min_rads)
b = np.unique(min_2_rads)
c = np.unique(min_3_rads)
nums_1 = len(np.where(min_rads == a[0])[0])
nums_2 = len(np.where(min_rads == a[1])[0])
nums_3 = len(np.where(min_rads == a[2])[0])
#nums_4 = len(np.where(min_rads == a[3])[0])
#nums_5 = len(np.where(min_rads == a[4])[0])


dbl_samp_names = np.array(['NGC 0596', 'NGC 1023','NGC 1399','NGC 1426',
                           'NGC 1700','NGC 2974','NGC 3585',
                           'NGC 3610','NGC 3640','NGC 3706','NGC 4026',
                           'NGC 4478',
                           'NGC 5061','NGC 5557','NGC 5576','NGC 7213',
                           'NGC 7785'])
#,'NGC 5018','NGC 4786','NGC 4125','NGC 2768','IC 3370']) # no data for excluded names


# =============================================================================
# plt.figure(dpi=500)
# #plt.plot(min_rads, min_2_rads, linestyle='', marker='.')
# plt.hist(min_rads, bins=30, alpha=0.7)
# plt.hist(min_2_rads, bins=200, alpha=0.7)
# #plt.hist(c, alpha=0.7)
# #plt.xlim(0,0.07)
# 
# plt.figure(dpi=500)
# plt.plot(min_rads, min_2_rads, linestyle='', marker='.',alpha=0.7,color='c')
# plt.xlabel('Innermost Radius ["]')
# plt.ylabel('2nd Innermost Radius["]')
# =============================================================================

# get the 5 different data samplings from our data
samp_1 = np.array([a[0], min_2_rads[np.where(min_rads == a[0])][0],
                  min_3_rads[np.where(min_rads == a[0])][0]])
y_1 = np.array([0,0,0])
samp_2 = np.array([a[1], min_2_rads[np.where(min_rads == a[1])][0],
                  min_3_rads[np.where(min_rads == a[1])][0]])
y_2 = np.array([1,1,1])
samp_3 = np.array([a[2], min_2_rads[np.where(min_rads == a[2])][0],
                  min_3_rads[np.where(min_rads == a[2])][0]])
y_3 = np.array([2,2,2])
#samp_4 = np.array([a[3], min_2_rads[np.where(min_rads == a[3])][0],
#                  min_3_rads[np.where(min_rads == a[3])][0]])
#y_4 = np.array([3,3,3])
#samp_5 = np.array([a[4], min_2_rads[np.where(min_rads == a[4])][0],
#                  min_3_rads[np.where(min_rads == a[4])][0]])
#y_5 = np.array([4,4,4])

plt.figure(dpi=500)
plt.title('Radial Sampling of our data (* indicates Lauer+95)')
plt.plot(samp_1,y_1,linestyle='',marker='.',label='{}'.format(nums_1))
plt.plot(samp_2,y_2,linestyle='',marker='.',label='{}*'.format(nums_2))
plt.plot(samp_3,y_3,linestyle='',marker='.',label='{}'.format(nums_3))
#plt.plot(samp_4,y_4,linestyle='',marker='.',label='{}*'.format(nums_4))
#plt.plot(samp_5,y_5,linestyle='',marker='.',label='{}'.format(nums_5))#
plt.plot([0.0456, 0.0456],[0,4],linestyle='--',color='k',alpha=0.3)
plt.text(0.047,1.5,'HST Pixel Scale',rotation=270,size=8)
plt.plot([0.0456, 0.0456],[0,4],linestyle='--',color='k',alpha=0.3)
plt.text(0.047,1.5,'HST Pixel Scale',rotation=270,size=8)
plt.legend()
plt.xlabel('Radius ["]')
plt.ylabel('Arbitrary Offsets')


# =============================================================================


#%%

# =============================================================================
# define gauss function to plot MGEs in pdf
def gauss(height, sig, r):
        return height*np.exp((-0.5)*r**2/sig**2)


# generate a nice PDF containing SB plots, MGE, Density profiles
fname_results = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_results_'+slope_ext+phys_ext+'.pdf'
with PdfPages(fname_results) as pdf:

    for i in range(len(names_05)):
        
        if all_lograds[i][0] <= np.log10(5):
        
            # let's plot the orginal SB data with the interpolated
            plt.figure(dpi=500)
            # add relevant text
            plt.text(0.05, 0.9, names_05[i], transform=plt.gcf().transFigure, size=10)
            plt.text(0.23, 0.9, 'M_v = {:.2f}'.format(lums_05[i]), transform=plt.gcf().transFigure, size=10)
            plt.text(0.43, 0.9, 'Dist = {:.2f} Mpc'.format(dists_05[i]/10**6), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.7,0.9, 'Inner Rad = {:.2f} pc'.format(radii_05[0,i]), 
                     transform=plt.gcf().transFigure, size=10)
            nsc_text = 'No NSC'
            if NSC_comp[i] == '1':
                nsc_text = 'NSC'
            elif NSC_comp[i] == '-1':
                nsc_text = 'Unknown'
            plt.text(0.9,0.02, nsc_text, 
                     transform=plt.gcf().transFigure, size=10)
            
            plt.plot(np.log10(rad_true_05[i]), SBs_true_05[i],linestyle='',marker='.',
                     label='Original')
            plt.plot(np.log10(rad_interp_05[:,i]), SBs_interp_05[:,i],label='Interpolated')
        
            plt.ylim(np.max(SBs_true_05[i])+0.2, np.min(SBs_true_05[i])-0.2)
            plt.ylabel('$\mu$ [mag/arcsec$^2$]')
            plt.xlabel('log(Radius [arcsec])')
            plt.legend()
        
            pdf.savefig()
            plt.close()
        
        
            # let's plot the MGE results 
            plt.figure(dpi=500)
            # add relevant text
            plt.text(0.05, 0.9, names_05[i], transform=plt.gcf().transFigure, size=10)
            plt.text(0.23, 0.9, 'M_v = {:.2f}'.format(lums_05[i]), transform=plt.gcf().transFigure, size=10)
            plt.text(0.43, 0.9, 'Dist = {:.2f} Mpc'.format(dists_05[i]/10**6), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.7,0.9, 'Inner Rad = {:.2f} pc'.format(radii_05[0,i]), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.9,0.02, nsc_text, 
                     transform=plt.gcf().transFigure, size=10)
            
            gausses = [] 
            for k in range(len(g_heights_05[i])):
                gausses.append(gauss(g_heights_05[i][k],g_sigs_05[i][k],rad_interp_05[:,i]))

            summed_gauss = np.zeros_like(rad_interp_05[:,i])
            for k in range(len(rad_interp_05[:,i])):
                for j in range(len(g_heights_05[i])):
                    summed_gauss[k] += gausses[j][k]

            plt.plot(rad_interp_05[:,i], SBs_interp_converted_05[:,i], color='b', linestyle='', marker='o', label='Data')
            for k in range(len(g_heights_05[i])):
                plt.plot(rad_interp_05[:,i],gausses[k],linestyle='--')

            plt.plot(rad_interp_05[:,i], summed_gauss, color='orange',linewidth=3.0, label='MGE fit')#, alpha=0.7)    
            
            plt.yscale('log')
            plt.xscale('log')
            plt.ylim(min(SBs_interp_converted_05[:,i])-
                     min(SBs_interp_converted_05[:,i])*0.5,
                     max(SBs_interp_converted_05[:,i])+
                     max(SBs_interp_converted_05[:,i])*0.5)
            plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
            plt.xlabel('Radius ["]')
            plt.legend()
            
            mge_resids = 100*(np.abs(SBs_interp_converted_05[:,i] - summed_gauss)/SBs_interp_converted_05[:,i])
            mge_mad = mad(mge_resids)
            plt.text(0.65,0.70, 'MAD = {:.4f}%'.format(mge_mad), 
                     transform=plt.gcf().transFigure, size=8)
            
            pdf.savefig()
            plt.close()



            # let's plot the MGE results with orginal data 
            plt.figure(dpi=500)
            # add relevant text
            plt.text(0.05, 0.9, names_05[i], transform=plt.gcf().transFigure, size=10)
            plt.text(0.23, 0.9, 'M_v = {:.2f}'.format(lums_05[i]), transform=plt.gcf().transFigure, size=10)
            plt.text(0.43, 0.9, 'Dist = {:.2f} Mpc'.format(dists_05[i]/10**6), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.7,0.9, 'Inner Rad = {:.2f} pc'.format(radii_05[0,i]), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.9,0.02, nsc_text, 
                     transform=plt.gcf().transFigure, size=10)

            plt.plot(rad_true_05[i], SBs_true_converted_05[i], color='b', linestyle='', marker='o')
            for k in range(len(g_heights_05[i])):
                plt.plot(rad_interp_05[:,i],gausses[k],linestyle='--')

            plt.plot(rad_interp_05[:,i], summed_gauss, color='orange', linewidth=2.0)#, alpha=0.7)    
            
            plt.yscale('log')
            plt.xscale('log')
            plt.ylim(min(SBs_interp_converted_05[:,i])-
                     min(SBs_interp_converted_05[:,i])*0.5,
                     max(SBs_interp_converted_05[:,i])+
                     max(SBs_interp_converted_05[:,i])*0.5)
            plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
            plt.xlabel('Radius ["]')
            
            mge_resamp = np.interp(rad_true_05[i], rad_interp_05[:,i], summed_gauss)
            
            mge_resids = 100*(np.abs(SBs_true_converted_05[i] - mge_resamp)/SBs_true_converted_05[i])
            mge_mad = mad(mge_resids)
            plt.text(0.65,0.75, 'MAD = {:.4f}%'.format(mge_mad), 
                     transform=plt.gcf().transFigure, size=8)
            
            pdf.savefig()
            plt.close()



            # let's plot the density profiles with stone's if available
            plt.figure(dpi=500)
            # add relevant text
            plt.text(0.05, 0.9, names_05[i], transform=plt.gcf().transFigure, size=10)
            plt.text(0.23, 0.9, 'M_v = {:.2f}'.format(lums_05[i]), transform=plt.gcf().transFigure, size=10)
            plt.text(0.43, 0.9, 'Dist = {:.2f} Mpc'.format(dists_05[i]/10**6), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.7,0.9, 'Inner Rad = {:.2f} pc'.format(radii_05[0,i]), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.9,0.02, nsc_text, 
                     transform=plt.gcf().transFigure, size=10)
            
            plt.plot(all_lograds[i],all_logdens[i],label='This Work')
            if len(stone_radii[i]) != 1:
                plt.plot(np.log10(stone_radii[i]),np.log10(stone_densities[i]),label='Stone&Metzger2016')
            plt.xlabel('log(Radius [pc])')
            plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
            plt.legend()
            plt.xlim(np.min(all_lograds[i])-0.2,np.max(all_lograds[i])+0.2)    
            
            if scale_mult*all_hst_ps[i] < phys_scale:
                rad_lim = phys_scale
            else:
                rad_lim = scale_mult*all_hst_ps[i]
            
            plt.fill_between([np.min(all_lograds[i]),np.log10(rad_lim)],
                             np.min(all_logdens[i]),np.max(all_logdens[i])+
                             0.2*np.max(all_logdens[i]), color='c',alpha=0.4)
            
            
            plt.text(0.63,0.66, 'Slope = {:.2f}'.format(slopes_05[i]), 
                     transform=plt.gcf().transFigure, size=8)
            plt.text(0.63,0.72, 'log(Dens@5pc) = {:.2f}'.format(np.log10(central_dens_05[i])), 
                     transform=plt.gcf().transFigure, size=8)
        
            if stone_dens_at_5pc[i] != -99:
                plt.text(0.63,0.63, 'Stone Slope = {:.2f}'.format(stone_slopes[i]), 
                         transform=plt.gcf().transFigure, size=8)
                plt.text(0.63,0.69, 'Stone log(Dens@5pc) = {:.2f}'.format(np.log10(stone_dens_at_5pc[i])), 
                         transform=plt.gcf().transFigure, size=8)
            
            pdf.savefig()
            plt.close()
        
        

    for i in range(len(names_95)):
        
        if all_lograds[len(names_05)+i][0] <= np.log10(5):
        
            # let's plot the orginal SB data with the interpolated
            plt.figure(dpi=500)
            # add relevant text
            plt.text(0.05, 0.9, names_95[i], transform=plt.gcf().transFigure, size=10)
            plt.text(0.23, 0.9, 'M_v = {:.2f}'.format(lums_95[i]), transform=plt.gcf().transFigure, size=10)
            plt.text(0.43, 0.9, 'Dist = {:.2f} Mpc'.format(dists_95[i]/10**6), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.7,0.9, 'Inner Rad = {:.2f} pc'.format(radii_95[0,i]), 
                     transform=plt.gcf().transFigure, size=10)
            nsc_text = 'No NSC'
            if NSC_comp[len(names_05)+i] == '1':
                nsc_text = 'NSC'
            elif NSC_comp[len(names_05)+i] == '-1':
                nsc_text = 'Unknown'
            plt.text(0.9,0.02, nsc_text, 
                     transform=plt.gcf().transFigure, size=10)
            
            plt.plot(np.log10(rad_true_95[i]), SBs_true_95[i],linestyle='',marker='.',
                     label='Original')
            plt.plot(np.log10(rad_interp_95[:,i]), SBs_interp_95[:,i],label='Interpolated')
            
            plt.ylim(np.max(SBs_true_95[i])+0.2, np.min(SBs_true_95[i])-0.2)
            plt.ylabel('$\mu$ [mag/arcsec$^2$]')
            plt.xlabel('log(Radius [arcsec])')
            plt.legend()
            
            pdf.savefig()
            plt.close()
            
            
            # let's plot the MGE results 
            plt.figure(dpi=500)
            # add relevant text
            plt.text(0.05, 0.9, names_95[i], transform=plt.gcf().transFigure, size=10)
            plt.text(0.23, 0.9, 'M_v = {:.2f}'.format(lums_95[i]), transform=plt.gcf().transFigure, size=10)
            plt.text(0.43, 0.9, 'Dist = {:.2f} Mpc'.format(dists_95[i]/10**6), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.7,0.9, 'Inner Rad = {:.2f} pc'.format(radii_95[0,i]), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.9,0.02, nsc_text, 
                     transform=plt.gcf().transFigure, size=10)
            
            gausses = [] 
            for k in range(len(g_heights_95[i])):
                gausses.append(gauss(g_heights_95[i][k],g_sigs_95[i][k],rad_interp_95[:,i]))
                
            summed_gauss = np.zeros_like(rad_interp_95[:,i])
            for k in range(len(rad_interp_95[:,i])):
                for j in range(len(g_heights_95[i])):
                    summed_gauss[k] += gausses[j][k]

            plt.plot(rad_interp_95[:,i], SBs_interp_converted_95[:,i], color='b', linestyle='', marker='o')
            for k in range(len(g_heights_95[i])):
                plt.plot(rad_interp_95[:,i],gausses[k])

            plt.plot(rad_interp_95[:,i], summed_gauss, color='orange',)#, alpha=0.7)    
            
            plt.yscale('log')
            plt.xscale('log')
            plt.ylim(min(SBs_interp_converted_95[:,i])-
                     min(SBs_interp_converted_95[:,i])*0.5,
                     max(SBs_interp_converted_95[:,i])+
                     max(SBs_interp_converted_95[:,i])*0.5)
            plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
            plt.xlabel('Radius ["]')
            
            mge_resids = 100*(np.abs(SBs_interp_converted_95[:,i] - summed_gauss)/SBs_interp_converted_95[:,i])
            mge_mad = mad(mge_resids)
            plt.text(0.65,0.75, 'MAD = {:.4f}%'.format(mge_mad), 
                     transform=plt.gcf().transFigure, size=8)
        
            pdf.savefig()
            plt.close()



            # let's plot the MGE results with original data 
            plt.figure(dpi=500)
            # add relevant text
            plt.text(0.05, 0.9, names_95[i], transform=plt.gcf().transFigure, size=10)
            plt.text(0.23, 0.9, 'M_v = {:.2f}'.format(lums_95[i]), transform=plt.gcf().transFigure, size=10)
            plt.text(0.43, 0.9, 'Dist = {:.2f} Mpc'.format(dists_95[i]/10**6), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.7,0.9, 'Inner Rad = {:.2f} pc'.format(radii_95[0,i]), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.9,0.02, nsc_text, 
                     transform=plt.gcf().transFigure, size=10)
            
            plt.plot(rad_true_95[i], SBs_true_converted_95[i], color='b', linestyle='', marker='o')
            for k in range(len(g_heights_95[i])):
                plt.plot(rad_interp_95[:,i],gausses[k])

            plt.plot(rad_interp_95[:,i], summed_gauss, color='orange',)#, alpha=0.7)    
            
            plt.yscale('log')
            plt.xscale('log')
            plt.ylim(min(SBs_interp_converted_95[:,i])-
                     min(SBs_interp_converted_95[:,i])*0.5,
                     max(SBs_interp_converted_95[:,i])+
                     max(SBs_interp_converted_95[:,i])*0.5)
            plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
            plt.xlabel('Radius ["]')
            
            mge_resamp = np.interp(rad_true_95[i], rad_interp_95[:,i], summed_gauss)
            
            mge_resids = 100*(np.abs(SBs_true_converted_95[i] - mge_resamp)/SBs_true_converted_95[i])
            mge_mad = mad(mge_resids)
            plt.text(0.65,0.75, 'MAD = {:.4f}%'.format(mge_mad), 
                     transform=plt.gcf().transFigure, size=8)
            
            pdf.savefig()
            plt.close()



            # let's plot the density profiles with stone's if available
            plt.figure(dpi=500)
            # add relevant text
            plt.text(0.05, 0.9, names_95[i], transform=plt.gcf().transFigure, size=10)
            plt.text(0.23, 0.9, 'M_v = {:.2f}'.format(lums_95[i]), transform=plt.gcf().transFigure, size=10)
            plt.text(0.43, 0.9, 'Dist = {:.2f} Mpc'.format(dists_95[i]/10**6), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.7,0.9, 'Inner Rad = {:.2f} pc'.format(radii_95[0,i]), 
                     transform=plt.gcf().transFigure, size=10)
            plt.text(0.9,0.02, nsc_text, 
                     transform=plt.gcf().transFigure, size=10)
            
            plt.plot(np.log10(radii_95[:,i]),np.log10(dens_profiles_95[:,i]),label='This Work')
            if len(stone_radii[len(names_05)+i]) != 1:
                plt.plot(np.log10(stone_radii[len(names_05)+i]),
                         np.log10(stone_densities[len(names_05)+i]),label='Stone&Metzger2016')
            plt.xlabel('log(Radius [pc])')
            plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
            plt.legend()
            plt.xlim(np.min(np.log10(radii_95[:,i]))-0.2,np.max(np.log10(radii_95[:,i]))+0.2)
        
            if scale_mult*all_hst_ps[len(names_05)+i] < phys_scale:
                rad_lim = phys_scale
            else:
                rad_lim = scale_mult*all_hst_ps[len(names_05)+i]
        
            plt.fill_between([np.min(np.log10(radii_95[:,i])),np.log10(rad_lim)],
                             np.min(np.log10(dens_profiles_95[:,i])),np.max(np.log10(dens_profiles_95[:,i]))+
                             0.2*np.max(np.log10(dens_profiles_95[:,i])), color='c',alpha=0.4)
        
            plt.text(0.63,0.66, 'Slope = {:.2f}'.format(slopes_95[i]), 
                     transform=plt.gcf().transFigure, size=8)
            plt.text(0.63,0.72, 'log(Dens@5pc) = {:.2f}'.format(np.log10(central_dens_95[i])), 
                     transform=plt.gcf().transFigure, size=8)
            if stone_dens_at_5pc[len(names_05)+i] != -99:
                plt.text(0.63,0.63, 'Stone Slope = {:.2f}'.format(stone_slopes[len(names_05)+i]), 
                         transform=plt.gcf().transFigure, size=8)
                plt.text(0.63,0.69, 'Stone log(Dens@5pc) = {:.2f}'.format(np.log10(stone_dens_at_5pc[len(names_05)+i])), 
                         transform=plt.gcf().transFigure, size=8)
        
            pdf.savefig()
            plt.close()
        
    
# =============================================================================



#%%

# =============================================================================
# gal_file = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data.fits'
# hdul = fits.open(gal_file)
# head = hdul[0].data
# dat = hdul[1].data
# hdul.close()
# 
# # =============================================================================
# #%%
# # =============================================================================
# # plot the galaxy distributions we have
# # distribution of luminosities
# plt.figure(dpi=500)
# plt.hist(lumys,bins=50)
# plt.xlabel('M$_v$ [mag]')
# plt.ylabel('# of Galaxies')
# 
# # distribution of masses
# plt.figure(dpi=500)
# plt.hist(gal_logmasses_1,bins=50)
# plt.xlabel('log(M$_*$ [M$_\odot$])')
# plt.ylabel('# of Galaxies')
# # =============================================================================
# #%%
# # =============================================================================
# # plot the distribution of innermost radii
# all_inner_radii = np.concatenate((radii_05[0,:],radii_95[0,:],radii_03[0,:]))
# plt.figure(dpi=500)
# plt.hist(all_inner_radii,bins=100)
# plt.xlabel('Innermost Radius [pc]')
# plt.ylabel('# of Galaxies')
# 
# # plot the distribution of innermost radii without the largest outlier
# all_inner_radii = np.concatenate((radii_05[0,:],radii_95[0,:],radii_03[0,:]))
# plt.figure(dpi=500)
# plt.hist(all_inner_radii[0:-1],bins=100)
# plt.xlabel('Innermost Radius [pc]')
# plt.ylabel('# of Galaxies')
# 
# arcsec_to_radian = 4.8481e-6
# rad_0_A2052 = (0.227*arcsec_to_radian)*150.9*10**6#120.21*10**6
# rad_0_A0119 = (0.023*arcsec_to_radian)*177.9*10**6#120.21*10**6
# # =============================================================================
# #%%
# # =============================================================================
# # plot all central_densities and slopes as a function of magnitude and galaxy_mass
# plt.figure(dpi=500)
# plt.plot(lumys, np.log10(central_dens), linestyle='', marker='.')
# plt.xlabel('M$_v$ [mag]')
# plt.ylabel('log(Innermost 3D Density)')
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Central_Density_and_Slope_Plots/central_dens_vs_Mv.pdf',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Central_Density_and_Slope_Plots/central_dens_vs_Mv.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# central_dens_w_mass = central_dens[mass_inds]
# plt.figure(dpi=500)
# plt.plot(gal_logmasses_1, np.log10(central_dens_w_mass), linestyle='', marker='.')
# plt.xlabel('log(M$_*$ [M$_\odot$])')
# plt.ylabel('log(Innermost 3D Density)')
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Central_Density_and_Slope_Plots/central_dens_vs_logmass.pdf',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Central_Density_and_Slope_Plots/central_dens_vs_logmass.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# plt.figure(dpi=500)
# plt.plot(lumys, slopes, linestyle='', marker='.')
# plt.xlabel('M$_v$ [mag]')
# plt.ylabel('Power-law ($\gamma$)')
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Central_Density_and_Slope_Plots/slope_vs_Mv.pdf',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Central_Density_and_Slope_Plots/slope_vs_Mv.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# slopes_w_mass = slopes[mass_inds]
# plt.figure(dpi=500)
# plt.plot(gal_logmasses_1, slopes_w_mass, linestyle='', marker='.')
# plt.xlabel('log(M$_*$ [M$_\odot$])')
# plt.ylabel('Power-law ($\gamma$)')
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Central_Density_and_Slope_Plots/slope_vs_logmass.pdf',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Central_Density_and_Slope_Plots/slope_vs_logmass.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# # =============================================================================
# #%%
# # =============================================================================
# # let's bin the density profiles together based on luminosity
# # define the bins
# bin_1 = np.array([-17,-15])
# bin_2 = np.array([-19,-17.00001])
# bin_3 = np.array([-21,-19.00001])
# bin_4 = np.array([-23,-21.00001])
# bin_5 = np.array([-25.4,-23.00001])
# 
# # get the indices of glaxies in each bin
# inds_1 = []
# inds_2 = []
# inds_3 = []
# inds_4 = []
# inds_5 = []
# for i in range(len(all_names_c)):
#     if lumys[i] > bin_1[0] and lumys[i] <= bin_1[1]:
#         inds_1.append(i)
#     elif lumys[i] > bin_2[0] and lumys[i] <= bin_2[1]:
#         inds_2.append(i)
#     elif lumys[i] > bin_3[0] and lumys[i] <= bin_3[1]:
#         inds_3.append(i)
#     elif lumys[i] > bin_4[0] and lumys[i] <= bin_4[1]:
#         inds_4.append(i)
#     elif lumys[i] > bin_5[0] and lumys[i] <= bin_5[1]:
#         inds_5.append(i)
# inds_1 = np.array(inds_1)
# inds_2 = np.array(inds_2)
# inds_3 = np.array(inds_3)
# inds_4 = np.array(inds_4)
# inds_5 = np.array(inds_5)
#     
# all_lograds = np.array(all_lograds)
# all_logdens = np.array(all_logdens)
# 
# rads_1 = all_lograds[inds_1]
# dens_1 = all_logdens[inds_1]
# rads_2 = all_lograds[inds_2]
# dens_2 = all_logdens[inds_2]
# rads_3 = all_lograds[inds_3]
# dens_3 = all_logdens[inds_3]
# rads_4 = all_lograds[inds_4]
# dens_4 = all_logdens[inds_4]
# rads_5 = all_lograds[inds_5]
# dens_5 = all_logdens[inds_5]
# 
# density_1 = []
# density_2 = []
# density_3 = []
# density_4 = []
# density_5 = []
# 
# num_rad = 200
# 
# r_min = np.min(rads_1)
# r_max = np.max(rads_1)
# radii_1 = np.linspace(r_min, r_max, num_rad)
# dens_1_i = np.zeros((len(inds_1),num_rad))
# mask = np.zeros((len(inds_1),len(radii_1)), dtype=bool)
# for i in range(len(inds_1)):
#     for j in range(len(rads_1[i])):
#         if rads_1[i][j] >= r_min:
#             break
#     mask[i,j:] = True
#     dens_1_i[i,:] = np.interp(radii_1, rads_1[i],dens_1[i])
# for i in range(len(radii_1)):
#     density_1[i] = np.mean(dens_1_i)
#     
#     
#     
#     
#     
# # =============================================================================
# 
# =============================================================================












