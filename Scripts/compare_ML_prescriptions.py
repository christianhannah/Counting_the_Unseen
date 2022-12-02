#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 15:26:38 2022

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


# Read in our data
slope_ext = '2x_pixel_scale'
phys_ext = '_or_10pc_extinction_corr_nsa_ml_w_vi_color'
names_c, slopes_c, cen_dens_c, vmags_c, stone_slopes_c, \
    stone_dens_at_5pc_c, NSC_comp_c, lograds_c, logdens_c, dists_c, \
    all_MLs, all_ML_types = u.get_our_data(slope_ext,phys_ext,True)

names_c[np.where(names_c == 'PGC 02888')] = 'PGC 028887'
names_c[np.where(names_c == 'PGC 05039')] = 'PGC 050395'


names_0 = []
MLs_0 = []
names_1 = []
MLs_1 = []
names_2 = []
MLs_2 = []
names_3 = []
MLs_3 = []
names_4 = []
MLs_4 = []
names_5 = []
MLs_5 = []
for i in range(len(names_c)):
    for j in range(len(all_ML_types[i,:])):
        if all_ML_types[i,j] == 0:
            names_0.append(names_c[i])
            MLs_0.append(all_MLs[i,j])
        if all_ML_types[i,j] == 1:
            names_1.append(names_c[i])
            MLs_1.append(all_MLs[i,j])
        if all_ML_types[i,j] == 2:
            names_2.append(names_c[i])
            MLs_2.append(all_MLs[i,j])
        if all_ML_types[i,j] == 3:
            names_3.append(names_c[i])
            MLs_3.append(all_MLs[i,j])
        if all_ML_types[i,j] == 4:
            names_4.append(names_c[i])
            MLs_4.append(all_MLs[i,j])
        if all_ML_types[i,j] == 5:
            names_5.append(names_c[i])
            MLs_5.append(all_MLs[i,j])
        
names_0 = np.array(names_0)
MLs_0 = np.array(MLs_0)
names_1 = np.array(names_1)
MLs_1 = np.array(MLs_1)
names_2 = np.array(names_2)
MLs_2 = np.array(MLs_2)
names_3 = np.array(names_3)
MLs_3 = np.array(MLs_3)
names_4 = np.array(names_4)
MLs_4 = np.array(MLs_4)
names_5 = np.array(names_5)
MLs_5 = np.array(MLs_5)      
        
            
#%%
# find galaxies with multiple M/Ls for plotting
idx_0_1 = []
idx_1_0 = []
names_0_1 = []
for i in range(len(names_0)):
    if names_0[i] in names_1:
        names_0_1.append(names_0[i])
        idx_0_1.append(i)
        idx_1_0.append(np.where(names_1 == names_0[i])[0])
idx_0_1 = np.array(idx_0_1)
idx_1_0 = np.array(idx_1_0)


idx_0_2 = []
idx_2_0 = []
names_0_2 = []
for i in range(len(names_0)):
    if names_0[i] in names_2:
        names_0_2.append(names_0[i])
        idx_0_2.append(i)
        idx_2_0.append(np.where(names_2 == names_0[i])[0])
idx_0_2 = np.array(idx_0_2)
idx_2_0 = np.array(idx_2_0)


idx_0_3 = []
idx_3_0 = []
names_0_3 = []
for i in range(len(names_0)):
    if names_0[i] in names_3:
        names_0_3.append(names_0[i])
        idx_0_3.append(i)
        idx_3_0.append(np.where(names_3 == names_0[i])[0])
idx_0_3 = np.array(idx_0_3)
idx_3_0 = np.array(idx_3_0)


idx_0_4 = []
idx_4_0 = []
names_0_4 = []
for i in range(len(names_0)):
    if names_0[i] in names_4:
        names_0_4.append(names_0[i])
        idx_0_4.append(i)
        idx_4_0.append(np.where(names_4 == names_0[i])[0])
idx_0_4 = np.array(idx_0_4)
idx_4_0 = np.array(idx_4_0)


idx_0_5 = []
idx_5_0 = []
names_0_5 = []
for i in range(len(names_0)):
    if names_0[i] in names_5:
        names_0_5.append(names_0[i])
        idx_0_5.append(i)
        idx_5_0.append(np.where(names_5 == names_0[i])[0])
idx_0_5 = np.array(idx_0_5)
idx_5_0 = np.array(idx_5_0)


idx_1_2 = []
idx_2_1 = []
names_1_2 = []
for i in range(len(names_1)):
    if names_1[i] in names_2:
        names_1_2.append(names_1[i])
        idx_1_2.append(i)
        idx_2_1.append(np.where(names_2 == names_1[i])[0])
idx_1_2 = np.array(idx_1_2)
idx_2_1 = np.array(idx_2_1)


        
idx_1_3 = []
idx_3_1 = []
names_1_3 = []
for i in range(len(names_1)):
    if names_1[i] in names_3:
        names_1_3.append(names_1[i])
        idx_1_3.append(i)
        idx_3_1.append(np.where(names_3 == names_1[i])[0])
idx_1_3 = np.array(idx_1_3)
idx_3_1 = np.array(idx_3_1)

      
idx_1_4 = []
idx_4_1 = []
names_1_4 = []
for i in range(len(names_1)):
    if names_1[i] in names_4:
        names_1_4.append(names_1[i])
        idx_1_4.append(i)
        idx_4_1.append(np.where(names_4 == names_1[i])[0])
idx_1_4 = np.array(idx_1_4)
idx_4_1 = np.array(idx_4_1)

  
idx_1_5 = []
idx_5_1 = []
names_1_5 = []
for i in range(len(names_1)):
    if names_1[i] in names_5:
        names_1_5.append(names_1[i])
        idx_1_5.append(i)
        idx_5_1.append(np.where(names_5 == names_1[i])[0])
idx_1_5 = np.array(idx_1_5)
idx_5_1 = np.array(idx_5_1)

      
idx_2_3 = []
idx_3_2 = []
names_2_3 = []
for i in range(len(names_2)):
    if names_2[i] in names_3:
        names_2_3.append(names_2[i])
        idx_2_3.append(i)
        idx_3_2.append(np.where(names_3 == names_2[i])[0])
idx_2_3 = np.array(idx_2_3)
idx_3_2 = np.array(idx_3_2)


idx_2_4 = []
idx_4_2 = []
names_2_4 = []
for i in range(len(names_2)):
    if names_2[i] in names_4:
        names_2_4.append(names_2[i])
        idx_2_4.append(i)
        idx_4_2.append(np.where(names_4 == names_2[i])[0])
idx_2_4 = np.array(idx_2_4)
idx_4_2 = np.array(idx_4_2)

        
idx_2_5 = []
idx_5_2 = []
names_2_5 = []
for i in range(len(names_2)):
    if names_2[i] in names_5:
        names_2_5.append(names_2[i])
        idx_2_5.append(i)
        idx_5_2.append(np.where(names_5 == names_2[i])[0])
idx_2_5 = np.array(idx_2_5)
idx_5_2 = np.array(idx_5_2)


idx_3_4 = []
idx_4_3 = []
names_3_4 = []
for i in range(len(names_3)):
    if names_3[i] in names_4:
        names_3_4.append(names_3[i])
        idx_3_4.append(i)
        idx_4_3.append(np.where(names_4 == names_3[i])[0])
idx_3_4 = np.array(idx_3_4)
idx_4_3 = np.array(idx_4_3)


idx_3_5 = []
idx_5_3 = []
names_3_5 = []
for i in range(len(names_3)):
    if names_3[i] in names_5:
        names_3_5.append(names_3[i])
        idx_3_5.append(i)
        idx_5_3.append(np.where(names_5 == names_3[i])[0])
idx_3_5 = np.array(idx_3_5)
idx_5_3 = np.array(idx_5_3)

       
idx_4_5 = []
idx_5_4 = []
names_4_5 = []
for i in range(len(names_4)):
    if names_4[i] in names_5:
        names_4_5.append(names_4[i])
        idx_4_5.append(i)
        idx_5_4.append(np.where(names_5 == names_4[i])[0])
idx_4_5 = np.array(idx_4_5)
idx_5_4 = np.array(idx_5_4)


# let's se how the galaxies fall with rankings
names_0_only = []
MLs_0_only = []
names_1_only = []
MLs_1_only = []
names_2_only = []
MLs_2_only = []
names_3_only = []
MLs_3_only = []
names_4_only = []
MLs_4_only = []
names_5_only = []
MLs_5_only = []
for i in range(len(names_c)):
    if 2 in all_ML_types[i,:]:
        names_2_only.append(names_c[i])
    elif 3 in all_ML_types[i,:]:
        names_3_only.append(names_c[i])
    elif 5 in all_ML_types[i,:]:
        names_5_only.append(names_c[i])
    elif 4 in all_ML_types[i,:]:
        names_4_only.append(names_c[i])
    elif 0 in all_ML_types[i,:]:
        names_0_only.append(names_c[i])
    elif 1 in all_ML_types[i,:]:
        names_1_only.append(names_c[i])
        
names_0_only = np.array(names_0_only)
MLs_0_only = np.array(MLs_0_only)
names_1_only = np.array(names_1_only)
MLs_1_only = np.array(MLs_1_only)
names_2_only = np.array(names_2_only)
MLs_2_only = np.array(MLs_2_only)
names_3_only = np.array(names_3_only)
MLs_3_only = np.array(MLs_3_only)
names_4_only = np.array(names_4_only)
MLs_4_only = np.array(MLs_4_only)
names_5_only = np.array(names_5_only)
MLs_5_only = np.array(MLs_5_only)      






plt.figure(dpi=500)
plt.title('0 vs 1')
plt.plot(MLs_0[idx_0_1], MLs_1[idx_1_0], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.ylabel('M/L$_v$ from Stone&Metzger2016 (virial)')
plt.xlabel('M/L$_v$ from Stone&Metzger2016 (virial)')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/0_v_1.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.title('0 vs 2')
plt.plot(MLs_0[idx_0_2], MLs_2[idx_2_0], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_v$ from Stone&Metzger2016 (virial)')
plt.ylabel('M/L$_r$ from McDermid2015')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/0_v_2.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


plt.figure(dpi=500)
plt.title('0 vs 3')
plt.plot(MLs_0[idx_0_3], MLs_3[idx_3_0], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_v$ from Stone&Metzger2016 (virial)')
plt.ylabel('M/L$_i$ from Taylor2011 Relation with V-I color')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/0_v_3.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================
# plt.figure(dpi=500)
# plt.title('0 vs 4')
# plt.plot(MLs_0[idx_0_4], MLs_4[idx_4_0], linestyle='', marker='.')
# plt.xlabel('M/L$_v$ from Stone&Metzger2016 (luminosity)')
# plt.ylabel('M/L$_I$ from Hoyer+22')
# =============================================================================

plt.figure(dpi=500)
plt.title('0 vs 5')
plt.plot(MLs_0[idx_0_5], MLs_5[idx_5_0], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_v$ from Stone&Metzger2016 (virial)')
plt.ylabel('M/L$_I$ from Taylor2011 Relation with g-i color')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/0_v_5.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


plt.figure(dpi=500)
plt.title('1 vs 2')
plt.plot(MLs_1[idx_1_2], MLs_2[idx_2_1], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_v$ from Stone&Metzger2016 (luminosity)')
plt.ylabel('M/L$_r$ from McDermid2015')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/1_v_2.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.title('1 vs 3')
plt.plot(MLs_1[idx_1_3], MLs_3[idx_3_1], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_v$ from Stone&Metzger2016 (luminosity)')
plt.ylabel('M/L$_i$ from Taylor2011 Relation with V-I color')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/1_v_3.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.title('1 vs 4')
plt.plot(MLs_1[idx_1_4], MLs_4[idx_4_1], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_v$ from Stone&Metzger2016 (luminosity)')
plt.ylabel('M/L$_I$ from Hoyer+22')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/1_v_4.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.title('1 vs 5')
plt.plot(MLs_1[idx_1_5], MLs_5[idx_5_1], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_v$ from Stone&Metzger2016 (luminosity)')
plt.ylabel('M/L$_I$ from Taylor2011 Relation with g-i color')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/1_v_5.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# =============================================================================
# plt.figure(dpi=500)
# plt.title('2 vs 3')
# plt.plot(MLs_2[idx_2_3], MLs_3[idx_3_2], linestyle='', marker='.')
# plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
# plt.xlabel('M/L$_r$ from McDermid2015')
# plt.ylabel('M/L$_i$ from Taylor2011 Relation with V-I color')
# =============================================================================

# =============================================================================
# plt.figure(dpi=500)
# plt.title('2 vs 4')
# plt.plot(MLs_2[idx_2_4], MLs_4[idx_4_2], linestyle='', marker='.')
# plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
# plt.xlabel('M/L$_r$ from McDermid2015')
# plt.ylabel('M/L$_I$ from Hoyer+22')
# =============================================================================

plt.figure(dpi=500)
plt.title('2 vs 5')
plt.plot(MLs_2[idx_2_5], MLs_5[idx_5_2], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_r$ from McDermid2015')
plt.ylabel('M/L$_I$ from Taylor2011 Relation with g-i color')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/2_v_5.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================
# plt.figure(dpi=500)
# plt.title('3 vs 4')
# plt.plot(MLs_3[idx_3_4], MLs_4[idx_4_3], linestyle='', marker='.')
# plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
# plt.xlabel('M/L$_i$ from Taylor2011 Relation with V-I color')
# plt.ylabel('M/L$_I$ from Hoyer+22')
# =============================================================================

plt.figure(dpi=500)
plt.title('3 vs 5')
plt.plot(MLs_3[idx_3_5], MLs_5[idx_5_3], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_i$ from Taylor2011 Relation with V-I color')
plt.ylabel('M/L$_I$ from Taylor2011 Relation with g-i color')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/3_v_5.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.title('4 vs 5')
plt.plot(MLs_4[idx_4_5], MLs_5[idx_5_4], linestyle='', marker='.')
plt.plot([1.,6.],[1.,6.],linestyle='--', color='r')
plt.xlabel('M/L$_I$ from Hoyer+22')
plt.ylabel('M/L$_I$ from Taylor2011 Relation with g-i color')
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/ML_comparisons/4_v_5.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


