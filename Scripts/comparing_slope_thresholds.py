#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 14:56:25 2022

@author: christian
"""

import matplotlib.pyplot as plt
import numpy as np
import mge1d_util as u
from astropy.io import fits
from astropy.io import ascii
import csv
import pdb
from astropy.stats import median_absolute_deviation as mad
import re
from astropy import modeling

# Specify extentions for different slope calculation ranges 
slope_ext = '2x_pixel_scale'
#slope_ext = '1.5x_pixel_scale'
#lope_ext = 'pixel_scale'
#slope_ext = '3pts'
#slope_ext = ''

#phys_ext = '_or_5pc'
phys_ext = '_or_10pc'
#phys_ext = '_or_15pc'
#phys_ext = '_or_20pc'
#phys_ext = '_or_30pc'
#phys_ext = '_or_40pc'
#phys_ext = ''
#phys_ext = '10pc'


# use file extensions to specify the scale multiplier and physical scale
if re.split('x_',slope_ext)[0] == 'pixel_scale':
    scale_mult = 1
elif re.split('x_',slope_ext)[0] == '3pts' or re.split('x_',slope_ext)[0] == '':
    scale_mult = 0
else:
    scale_mult = float(re.split('x_',slope_ext)[0])
    
if len(re.split('_or_|pc', phys_ext)) == 1:
    phys_scale = 0
elif len(re.split('_or_|pc', phys_ext)) == 2:
    phys_scale = float(re.split('_or_|pc', phys_ext)[0])
else:
    phys_scale = float(re.split('_or_|pc', phys_ext)[1])


# Read in our data
names_c, slopes_c, cen_dens_c, vmags_c, stone_slopes_c, \
    stone_dens_at_5pc_c, NSC_comp_c, lograds_c, logdens_c, \
        dists_c = u.get_our_data(slope_ext,phys_ext,True)


# get all hst pixel scales in pc
# conversion factor between radians and arcsec
arcsec_to_radian = 4.8481e-6
all_scales = np.zeros(len(dists_c))
for i in range(len(all_scales)):
    all_scales[i] = (0.0456*arcsec_to_radian)*dists_c[i]*10**6
    
# =============================================================================
# plt.figure(dpi=500)
# plt.hist(all_scales,bins=30)
# plt.xlabel('Physical Pixel Scale 0.0456" [pc]')   
# plt.ylabel('# of Galaxies') 
# =============================================================================

# recompute the powerlaw slope of each density profile to check residuals
slopes_recomp = []
res_mads = []
residuals = np.array([])
for i in range(len(names_c)):
    s, cd, res_mad, res = u.get_density_and_slope(lograds_c[i,:], logdens_c[i,:], all_scales[i],
                                    phys_scale, scale_mult)
    if np.abs(s-slopes_c[i]) < 0.001:
        slopes_recomp.append(s)
        res_mads.append(res_mad)
        #residuals.append(res)
        residuals = np.concatenate((residuals,res))
    else:
        print('Slope Calc Failure for '+names_c[i])

slopes_recomp = np.array(slopes_recomp)
res_mads = np.array(res_mads)
#residuals = np.array(residuals)


SB_slope = u.get_density_and_slope(lograds_c[i,:], logdens_c[i,:], all_scales[i],
                                    phys_scale, scale_mult)




# plot the slopes from our galaxies vs stone profile slopes
k = np.where(stone_slopes_c != -99)
p = np.where(NSC_comp_c == 0)
w = np.intersect1d(k,p)

stone_slopes_corrected = stone_slopes_c[k]
slopes_c_corrected = slopes_c[k]
stone_slopes_no_nsc = stone_slopes_c[w]
slopes_no_nsc = slopes_c[w]

res_mads_corrected = res_mads[k]
res_mads_no_nsc = res_mads[w]

differences = stone_slopes_no_nsc - slopes_no_nsc
slope_mad = mad(np.abs(differences))
slope_med = np.median(differences)

print('STD = ',np.std(differences))
# old plot without the distribution
# =============================================================================
# plt.figure(dpi=500)
# plt.plot(stone_slopes_corrected,slopes_c_corrected,linestyle='',marker='.',
#          color='b')
# plt.plot(stone_slopes_no_nsc,slopes_no_nsc,linestyle='',marker='d',
#          color='c', label='No NSC')
# plt.plot([-2,0],[-2,0],linestyle='--',color='r')
# plt.xlabel('$\gamma$ from Stone&Metzger 2016')
# plt.ylabel('$\gamma$ from This Work')
# plt.ylim(-3.6,0.1)
# #plt.text(-0.37, np.min(slopes_c_corrected)+0.2,'MAD = {:.2f}'.format(slope_mad))
# plt.text(-0.37, -3.2,'MAD = {:.2f}'.format(slope_mad))
# plt.legend()
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Slope_Calc_Investigation/Stone_vs_Us_comparisons/slopes_'+slope_ext+phys_ext+'.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================

org = (0.9,0.3,0.)
nav = (0.,0.,0.52)

fig, ax1 = plt.subplots(dpi=500)
ax1.plot(stone_slopes_corrected,slopes_c_corrected,linestyle='',marker='o',
         color=nav, label='NSC')
ax1.plot(stone_slopes_no_nsc,slopes_no_nsc,linestyle='',marker='o',
         color=org, label='No NSC')
ax1.plot([-2,0],[-2,0],linestyle='--',color='r')
ax1.set_xlabel('$\gamma$ from Stone&Metzger 2016')
ax1.set_ylabel('$\gamma$ from This Work')
ax1.set_ylim(-3.8,0.1)
ax1.legend(loc='lower left')

if re.split('_|x', slope_ext)[0] =='pi':
    slope_part = '0.0456"'
elif re.split('_|x', slope_ext)[0] =='3pts':
    slope_part = '3 inner points'
else:
    slope_part = '{:.3f}"'.format(float(re.split('_|x', slope_ext)[0])*0.0456)
if len(phys_ext.split('_')) == 1:
    phys_part = ''
else:
    phys_part = phys_ext.split('_')[1]+' '+phys_ext.split('_')[2]
# =============================================================================
# slope_part = phys_ext
# phys_part = ''
# =============================================================================
ax1.set_title('Slope Range: '+slope_part+' '+phys_part)
#ax1.text(-0.37, -2,'MAD = {:.2f}'.format(slope_mad), size=8)
le, bo, hi, wi = .68, .19, .15, .2
ax2 = fig.add_axes([le, bo, wi, hi])
ax2.hist(differences,bins=10)
ax2.set_title('Median = {:.3f}'.format(slope_med)+
              '\n MAD = {:.3f}'.format(slope_mad), size=8)
ax2.set_xlim(-1.2,1.2)

fig.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Slope_Calc_Investigation/Stone_vs_Us_comparisons/slopes_'+slope_ext+phys_ext+'.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)



# plot the MADs of residuals from the power law fits to determine which range is best

plt.figure(dpi=500)
plt.title('Slope Range: '+slope_part+' '+phys_part)
plt.hist(res_mads,bins=30,color='0.4',alpha=0.7,label='ALL')
plt.hist(res_mads_no_nsc,bins=20,color='b',alpha=0.7,label='No NSC')
plt.xlabel('MAD of Power-law Fit Residuals')
plt.legend()
plt.xlim(0,0.25)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Slope_Calc_Investigation/Stone_vs_Us_comparisons/fit_residual_mads_'+slope_ext+phys_ext+'.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# plot all the residuals from power law fits
plt.figure(dpi=500)
plt.title('Slope Range: '+slope_part+' '+phys_part)
H = plt.hist(residuals,bins=500)
plt.xlabel('Power-law Fit Residuals')
plt.xlim(-0.5,0.5)

# add gaussian fit to the histogram
y_hist = H[0]
x_hist = np.zeros_like(y_hist)
for i in range(len(y_hist)):
    x_hist[i] = H[1][i]+np.abs(H[1][i+1]-H[1][i])/2

fitter = modeling.fitting.LevMarLSQFitter()
model = modeling.models.Gaussian1D()   # depending on the data you need to give some initial values
fit_gauss = fitter(model, x_hist, y_hist)

plt.plot(x_hist,fit_gauss(x_hist),color='r',alpha=0.7)

# =============================================================================
# plt.text(0.15, 0.8, 'Mean: {:.3f}'.format(fit_gauss.mean.value), transform=plt.gcf().transFigure, size=10)
# plt.text(0.167, 0.76, 'STD: {:.3f}'.format(fit_gauss.stddev.value), transform=plt.gcf().transFigure, size=10)
# =============================================================================
plt.text(0.15, 0.8, 'Median: {:.3f}'.format(np.median(residuals)), transform=plt.gcf().transFigure, size=10)
plt.text(0.167, 0.76, 'MAD: {:.3f}'.format(mad(residuals)), transform=plt.gcf().transFigure, size=10)

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Slope_Calc_Investigation/Stone_vs_Us_comparisons/fit_residual_distributions_'+slope_ext+phys_ext+'.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


#%%
# plot the densities from our galaxies vs stone profile densities
stone_dens_at_5pc_c_corrected = stone_dens_at_5pc_c[k]
cen_dens_c_corrected = cen_dens_c[k]
stone_dens_at_5pc_no_nsc = stone_dens_at_5pc_c[w]
cen_dens_no_nsc = cen_dens_c[w]

# =============================================================================
# plt.figure(dpi=500)
# plt.plot(np.log10(stone_dens_at_5pc_c_corrected),np.log10(cen_dens_c_corrected),
#          linestyle='',marker='.',color='b')
# plt.plot([2,6],[2,6],linestyle='--',color='r')
# plt.plot(np.log10(stone_dens_at_5pc_no_nsc),np.log10(cen_dens_no_nsc),
#          linestyle='',marker='d',color='c', label='No NSC')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from Stone&Metzger 2016')
# plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from This Work')
# plt.legend()
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Slope_Calc_Investigation/Stone_vs_Us_comparisons/densities_'+slope_ext+phys_ext+'.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================

den_differences = np.log10(stone_dens_at_5pc_no_nsc) - np.log10(cen_dens_no_nsc)
den_mad = mad(np.abs(den_differences))
den_med = np.median(den_differences)

fig, ax1 = plt.subplots(dpi=500)
ax1.plot(np.log10(stone_dens_at_5pc_c_corrected),np.log10(cen_dens_c_corrected),
         linestyle='',marker='o',color=nav, label='NSC')
ax1.plot([2,6],[2,6],linestyle='--',color='r')
ax1.plot(np.log10(stone_dens_at_5pc_no_nsc),np.log10(cen_dens_no_nsc),
         linestyle='',marker='o',color=org, label='No NSC')
ax1.set_xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from Stone&Metzger 2016')
ax1.set_ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from This Work')
slope_part = slope_ext.split('_')[0]
# =============================================================================
# if slope_ext.split('_')[0] =='pixel':
#     slope_part = '0.05"'
# elif slope_ext.split('_')[0] =='3pts':
#     slope_part = '3 inner points'
# else:
#     slope_part = slope_ext.split('_')[0]+'0.05"'
# if len(phys_ext.split('_')) == 1:
#     phys_part = ''
# else:
#     phys_part = phys_ext.split('_')[1]+' '+phys_ext.split('_')[2]
# ax1.set_title('Slope Range: '+slope_part+' '+phys_part)
# =============================================================================
ax1.set_ylim(2,6)
ax1.set_xlim(2,6)
#ax1.text(-0.37, -2,'MAD = {:.2f}'.format(slope_mad), size=8)
le, bo, hi, wi = .68, .19, .15, .2
ax2 = fig.add_axes([le, bo, wi, hi])
ax2.hist(den_differences,bins=10)
#ax2.set_title('MAD = {:.3f}'.format(den_mad), size=8)
ax2.set_title('Median = {:.3f}'.format(den_med)+
              '\n MAD = {:.3f}'.format(den_mad), size=8)
ax2.set_xlim(-1.2,1.2)

fig.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Slope_Calc_Investigation/Stone_vs_Us_comparisons/densities.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

#%%

# stand alone MAD and median plots
# order of data: mad(#pixelscale) -> [5pc,10pc,15pc,20pc,30pc,40pc]
mad1 = np.array([0.166,0.150,0.105,0.087,0.098,0.119])
mad15 = np.array([0.118,0.148,0.105,0.087,0.098,0.119])
mad2 = np.array([0.104,0.123,0.095,0.087,0.098,0.119])

med1 = np.array([-.352,-.1,-.141,-.120,-.1,-.094])
med15 = np.array([-.186,-.1,-.141,-.120,-.1,-.094])
med2 = np.array([-.166,-.1,-.141,-.120,-.1,-.094])

x = np.array([5,10,15,20,30,40])

plt.figure(dpi=500)
plt.plot(x,mad1, linestyle='', marker='d', color='b', label='0.046"')#, alpha=0.6)
plt.plot(x,mad15, linestyle='', marker='P', color='c', label='0.068"')#, alpha=0.6)
plt.plot(x,mad2, linestyle='', marker='*', color='m',label='0.091"')#, alpha=0.6)

# =============================================================================
# plt.plot([np.median(all_scales),np.median(all_scales)],[np.min(mad1),np.max(mad1)],
#          color='b',linestyle='--')
# plt.plot([np.median(1.5*all_scales),np.median(1.5*all_scales)],[np.min(mad1),np.max(mad1)],
#          color='g',linestyle='--')
# plt.plot([np.median(2*all_scales),np.median(2*all_scales)],[np.min(mad1),np.max(mad1)],
#          color='m',linestyle='--')
# =============================================================================

plt.xlabel('Max Physical Range [pc]')
plt.ylabel('MAD of Residuals')
plt.legend()

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Slope_Calc_Investigation/Stone_vs_Us_comparisons/all_mads.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.plot(x,med1, linestyle='', marker='d', color='b', label='0.046"')#, alpha=0.6)
plt.plot(x,med15, linestyle='', marker='P', color='c', label='0.068"')#, alpha=0.6)
plt.plot(x,med2, linestyle='', marker='*', color='m',label='0.091"')#, alpha=0.6)

# =============================================================================
# plt.plot([np.median(all_scales),np.median(all_scales)],[np.min(mad1),np.max(mad1)],
#          color='b',linestyle='--')
# plt.plot([np.median(1.5*all_scales),np.median(1.5*all_scales)],[np.min(mad1),np.max(mad1)],
#          color='g',linestyle='--')
# plt.plot([np.median(2*all_scales),np.median(2*all_scales)],[np.min(mad1),np.max(mad1)],
#          color='m',linestyle='--')
# =============================================================================

plt.xlabel('Max Physical Range [pc]')
plt.ylabel('Median of Residuals')
plt.legend()

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Slope_Calc_Investigation/Stone_vs_Us_comparisons/all_medians.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)









