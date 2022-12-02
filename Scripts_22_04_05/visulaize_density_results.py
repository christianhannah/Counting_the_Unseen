#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 10:20:02 2021

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

import warnings
warnings.filterwarnings("ignore")

data = ascii.read('/Users/christian/OneDrive/Desktop/TDE Code/Pechetti_Table_A3.txt', format='latex') 
data_1 = ascii.read('/Users/christian/OneDrive/Desktop/TDE Code/Pechetti_Table_1.txt', format='latex') 
data_2 = ascii.read('/Users/christian/OneDrive/Desktop/TDE Code/Pechetti_Table_A1.txt', format='latex') 

#%%


# from Pechetti+20
gal_name = ['Circinus','ESO274-1','IC5052','IC5332','NGC2784','NGC2787',
             'NGC2903','NGC3115','NGC3115B','NGC3184','NGC3274','NGC3344',
             'NGC3593','NGC4242','NGC4460','NGC4517','NGC4592','NGC4600',
             'NGC4605','NGC4941','NGC5055','NGC5068','NGC5194','NGC5195',
             'NGC5236','NGC5238','NGC5457','NGC6503','NGC7713']

# in M_sol
gal_mass = np.array([3.06e10,4.28e8,1.35e9,7.59e9,5.13e10,1.19e10,4.57e10,
                     6.76e10,8.9e8,1.95e10,3.24e8,6.31e9,6.03e9,2.0e9,2.75e9,
                     1.01e10,5.3e8,1.19e9,3.55e9,3.17e9,4.9e10,6.31e9,6.03e10,
                     1.95e10,4.47e10,1.17e8,4.47e10,4.37e9,3.16e9])




# in M_sol
# =============================================================================
# NSC_mass_1 = np.array([3.71e7,2.73e6,4.03e6,6.96e6,-1,2.71e7,7.83e7,-1,7.19e6,-1,
#                        1.34e6,1.8e7,1.58e8,7.98e5,-1,2.86e6,4.71e6,1.03e7,2.4e7,
#                 6.02e7,5.12e7,1.66e6,1.61e8,1.09e7,1.67e6,7.21e5,2.18e7,
#                 1.19e7,6.53e5])
# =============================================================================
NSC_mass = np.array(data['NSC mass'][2:])
for i in range(len(NSC_mass)):
    if NSC_mass[i] != '-':
        s = NSC_mass[i].split('$\\times$$10^')
        s[1] = s[1].split('$')[0]
        NSC_mass[i] = s[0]+'e'+s[1]
    else: 
        NSC_mass[i] = -1
NSC_mass = np.array(NSC_mass).astype(np.float)

# in pc
# =============================================================================
# NSC_reff_1 = np.array([8.,2.18,3.73,23.18,20.37,5.12,10.32,26.89,6.61,2.05,3.32,4.79,5.5,
#                      1.74,4.67,5.73,3.45,8.27,2.19,14.52,14.61,5.18,24.15,9.07,4.91,2.15,
#                      10.88,1.62,1.80])
# 
# =============================================================================
NSC_reff = np.array(data_1['R$_e$ (pc)'][1:])
NSC_r_perr = np.zeros_like(NSC_reff)
NSC_r_nerr = np.zeros_like(NSC_reff)
for i in range(len(NSC_reff)):
    s = NSC_reff[i].split('$_{')[0]
    NSC_reff[i] = s
NSC_reff = np.array(NSC_reff).astype(np.float)

# in log(M_sol/pc^3)
# =============================================================================
# log_dens_5pc_1 = np.array([2.31,2.21,2.43,2.6,-1,3.04,3.,-1,2.37,-1,1.97,2.5,-1,3.51,-1,
#                            2.78,2.42,3.95,2.55,2.11,2.77,3.34,-1,3.17,3.07,3.01,3.18,3.25,
#                            3.96])
# =============================================================================
log_dens_5pc = data['log($\\rho_{5pc}$)'][2:]
renuka_bmags = data_2['M$_B$'][1:]
renuka_ttypes = data_2['T-type'][1:]

for i in range(len(log_dens_5pc)):
    if i == 3: #get rid of IC5332 which is not in published version
        log_dens_5pc[i] = '-'
    if log_dens_5pc[i] != '-':
        log_dens_5pc[i] = float(log_dens_5pc[i])
        renuka_bmags[i] = float(renuka_bmags[i])
        renuka_ttypes[i] = int(renuka_ttypes[i])
    else:
        log_dens_5pc[i] = -1
        renuka_bmags[i] = -99
        renuka_ttypes[i] = -99
log_dens_5pc = np.array(log_dens_5pc).astype(np.float)
renuka_bmags = np.array(renuka_bmags).astype(np.float)
renuka_ttypes = np.array(renuka_ttypes).astype(int)

gamma_1 = -1*np.array([2.77,2.69,2.19,2.87,-1,1.86,2.04,-1,2.10,-1,2.84,2.99,-1,1.81,-1,
                       2.77,2.9,1.44,1.62,2.53,2.6,1.44,-1,1.56,0.93,2.48,1.79,1.59,1.89])
gamma = data['$\\gamma$'][2:]
for i in range(len(gamma)):
    if gamma[i] != '-':
        gamma[i] = float(gamma[i])
    else:
        gamma[i] = 1
gamma = np.array(gamma).astype(np.float)


morphs_r = np.zeros(len(renuka_ttypes))
for i in range(len(morphs_r)):
    if renuka_ttypes[i] <= 0:
        morphs_r[i] = 0
    else:
        morphs_r[i] = 1

# =============================================================================
# Read in our data
slope_ext = '2x_pixel_scale'
#slope_ext = '1.5x_pixel_scale'
#slope_ext = 'pixel_scale'
#slope_ext = '3pts'
phys_ext = '_or_10pc'
names_c, slopes_c, cen_dens_c, vmags_c, stone_slopes_c, \
    stone_dens_at_5pc_c, NSC_comp_c, lograds_c, logdens_c, dists_c = u.get_our_data(slope_ext,phys_ext,True)

# get the log(masses) from David for the galaxies he has in his table
logmass_c = []
for i in range(len(names_c)):
    n = names_c[i].split(' ')
    if len(n) == 1:
        logmass_c.append(u.get_galaxy_logmass(n))
    else:
        logmass_c.append(u.get_galaxy_logmass(n[0]+n[1]))
logmass_c = np.array(logmass_c)

# get index list for galaxies that had logmasses
mass_inds = []
for i in range(len(logmass_c)):
    if logmass_c[i] != 0:
        mass_inds.append(i)


# =============================================================================
# gal_file = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data_100000pc.fits'
# hdul = fits.open(gal_file)
# head = hdul[0].data
# dat = hdul[1].data
# hdul.close()
# 
# names = dat['name']
# vmags = dat['vmag']
# dists = dat['dist'] # Mpc
# MLs = dat['ml']
# ML_types = dat['ml_type']
# slopes = dat['slope']
# cen_dens = dat['dens_at_5pc'] # M_sol/pc^3
# lograds = dat['lograd'] # log(pc)
# logdens = dat['logdens'] # log(M_sol/pc^3)
# stone_slopes = dat['stone_slope']
# stone_dens_at_5pc = dat['stone_dens_at_5pc']
# NSC_comp = dat['NSC_comp']
# 
# names_c = names[np.where(lograds[:,0] <= np.log10(5))]
# slopes_c = slopes[np.where(lograds[:,0] <= np.log10(5))]
# cen_dens_c = cen_dens[np.where(lograds[:,0] <= np.log10(5))]
# vmags_c = vmags[np.where(lograds[:,0] <= np.log10(5))]
# stone_slopes_c = stone_slopes[np.where(lograds[:,0] <= np.log10(5))]
# stone_dens_at_5pc_c = stone_dens_at_5pc[np.where(lograds[:,0] <= np.log10(5))]
# NSC_comp_c = NSC_comp[np.where(lograds[:,0] <= np.log10(5))]
# 
# =============================================================================
# =============================================================================

#%%
# =============================================================================
# read in the data for quick and dirty mag conversion
file = open('B-V_vs_Hubble_type_data_RobertsHaynes94.csv', 'r')
csvreader = csv.reader(file)
type_int = []
B_V = []
for row in csvreader:
    type_int.append(int(float(row[0])))
    B_V.append(float(row[1]))

# conversion arrays between integer morphology and lauer designation
# [1,2,3,4,5,6,7,8,9,10,11,12]
# [E, S0, S0/a, Sa, Sab, Sb, Sbc, Sc, Scd, Sd, Sm, Im]

# read in fits table with general galaxy properties including R_eff from Lauer+2007
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2007_data/gal_properties_Lauer_2007.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMag','logR','r1','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data_07 = hdul[1].data 
hdul.close()

morphs_c = np.zeros(len(names_c)).astype(int)
type_ind = np.zeros(len(names_c)).astype(int)
for i in range(len(names_c)):
    gal_ind = np.where(data_07['Name'] == names_c[i])
    morph = data_07['Morph'][gal_ind]
    if morph == 'BCG' or morph == 'E' or morph == 'S0-' or morph == 'cE':
        morphs_c[i] = 0
        type_ind[i] = 0
    elif morph == 'S0':
        morphs_c[i] = 0
        type_ind[i] = 1
    elif morph == 'S0+':
        morphs_c[i] = 0
        type_ind[i] = 2
    elif morph == 'Sa':
        morphs_c[i] = 1
        type_ind[i] = 3
    elif morph == 'Sb':
        morphs_c[i] = 1
        type_ind[i] = 5
    else:
        morphs_c[i] = -1
        type_ind[i] = -1

# temporary morphology assignment
morphs_c[np.where(morphs_c == -1)] = 0

bmags_c = np.zeros_like(vmags_c)
for i in range(len(vmags_c)):
    c_fac = B_V[type_ind[i]]
    bmags_c[i] = c_fac + vmags_c[i]
#bmags_c = bmags_c[np.where(lograds[:,0] >= np.log10(5))]

# =============================================================================
#%%


# =============================================================================
# Plot slopes and densities v. bmags from renuka and me colored by early and late types

# modify arrays to only include gals with the appropriate measurments
a = np.where(log_dens_5pc != -1)
#slopes = gamma[a]
gal_mass = gal_mass[a]
#cen_dens = log_dens_5pc[a]
bmags = renuka_bmags[a]
morphs_r = morphs_r[a]
#%%
slopes = -1*np.array([.92,2.16,2.06,1.42,1.56,1.86,2.63,1.43,2.80,2.44,2.87,
                      2.09,2.75,1.8,1.87,1.68,1.23,2.48,2.76,1.82,2.75,3.22])

cen_dens = np.array([3.84,2.78,3.06,3.93,4.03,3.21,2.35,3.83,2.03,2.61,2.75,
                     3.22,3.56,3.75,3.6,2.68,3.93,3.66,2.04,3.47,3.42,1.83])
#%%
plt.figure(dpi=500)
plt.scatter(bmags,slopes,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
plt.scatter(bmags_c, slopes_c,c=morphs_c,marker='d',cmap='viridis',label='This Work')
#plt.plot(bmags,slopes,linestyle='',marker='.',color='b',label='Pechetti+20')
#plt.plot(bmags_c, slopes_c,linestyle='',marker='d',color='c',label='This Work')
plt.ylabel('$\gamma$')
plt.xlabel('M$_B$')
plt.colorbar(label='0 = Early-type, 1 = Late-type')
#plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.legend()
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_and_Our_Sample/slopes_v_bmags_colored_by_morph.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.scatter(bmags,cen_dens,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
plt.scatter(bmags_c, np.log10(cen_dens_c),c=morphs_c,marker='d',cmap='viridis',label='This Work')
#plt.plot(bmags,cen_dens,linestyle='',marker='.',color='b',label='Pechetti+20')
#plt.plot(bmags_c, np.log10(cen_dens_c),linestyle='',marker='.',color='c',label='This Work')
#plt.ylabel('$\gamma$')
plt.xlabel('M$_B$')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.colorbar(label='0 = Early-type, 1 = Late-type')
plt.legend()
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_and_Our_Sample/densities_v_bmags_colored_by_morph.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# =============================================================================

#%%
# =============================================================================
# plot the slopes vs galaxy mass

#q = np.where(NSC_comp_c == 1)
#mass_inds = np.intersect1d(mass_inds,q)

logmass_corrected = logmass_c[mass_inds]
slopes_corrected = slopes_c[mass_inds]
cen_dens_corrected = cen_dens_c[mass_inds]
morphs_corrected = morphs_c[mass_inds]


plt.figure(dpi=500)
plt.scatter(np.log10(gal_mass),slopes,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
plt.plot(np.sort(np.log10(gal_mass)),0.42*np.log10(np.sort(gal_mass)/10**9)-2.39,linestyle='--',color='r')
plt.ylabel('$\gamma$')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.colorbar(label='0 = Early-type, 1 = Late-type')
#plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.legend()
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_and_Our_Sample/slopes_v_galmass_colored_by_morph.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.scatter(np.log10(gal_mass),cen_dens,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
plt.scatter(logmass_corrected, np.log10(cen_dens_corrected),c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
plt.plot(np.sort(np.log10(gal_mass)),0.61*np.log10(np.sort(gal_mass)/10**9)+2.78,linestyle='--',color='r')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.colorbar(label='0 = Early-type, 1 = Late-type')
#plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.legend()
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_and_Our_Sample/densities_v_galmass_colored_by_morph.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)



# =============================================================================

#%%


# =============================================================================
# plot the slopes vs galaxy mass

q = np.where(NSC_comp_c == 1)
#mass_nsc_inds = np.intersect1d(mass_inds,q)

# =============================================================================
# logmass_corrected = logmass_c[mass_inds]
# slopes_corrected = slopes_c[mass_inds]
# cen_dens_corrected = cen_dens_c[mass_inds]
# morphs_corrected = morphs_c[mass_inds]
# =============================================================================

# divide samples based on type for plotting
lt_inds_r = np.where(morphs_r == 1.)
et_inds_r = np.where(morphs_r == 0.)
gal_mass_lt = gal_mass[lt_inds_r]
gal_mass_et = gal_mass[et_inds_r]
slopes_lt = slopes[lt_inds_r]
slopes_et = slopes[et_inds_r]
cen_dens_lt = cen_dens[lt_inds_r]
cen_dens_et = cen_dens[et_inds_r]

lt_inds_c = np.where(morphs_corrected == 1.)
et_inds_c = np.where(morphs_corrected == 0.)
logmass_c_lt = logmass_corrected[lt_inds_c]
logmass_c_et = logmass_corrected[et_inds_c]
slopes_c_lt = slopes_corrected[lt_inds_c]
slopes_c_et = slopes_corrected[et_inds_c]
cen_dens_c_lt = cen_dens_corrected[lt_inds_c]
cen_dens_c_et = cen_dens_corrected[et_inds_c]

logmass_c_ltnsc = logmass_corrected[np.intersect1d(lt_inds_c,q)]
logmass_c_etnsc = logmass_corrected[np.intersect1d(et_inds_c,q)]
all_ltnsc_mass = np.concatenate((np.log10(gal_mass_lt),logmass_c_ltnsc))
all_etnsc_mass = np.concatenate((np.log10(gal_mass_et),logmass_c_etnsc))

def power_law(x,a,b):
    return a*x+b
# get linear relationships for early-types and late-types 
all_masses = np.concatenate((np.log10(gal_mass),logmass_corrected))
all_lt_mass = np.concatenate((np.log10(gal_mass_lt),logmass_c_lt))
all_lt_slopes = np.concatenate((slopes_lt, slopes_c_lt))
all_dens_lt = np.concatenate((cen_dens_lt, np.log10(cen_dens_c_lt)))
all_et_mass = np.concatenate((np.log10(gal_mass_et),logmass_c_et))
all_et_slopes = np.concatenate((slopes_et, slopes_c_et))
all_dens_et = np.concatenate((cen_dens_et, np.log10(cen_dens_c_et)))


lt_sort = np.argsort(all_lt_mass)
et_sort = np.argsort(all_et_mass)
ltm = all_lt_mass[lt_sort]
etm = all_et_mass[et_sort]
lts = all_lt_slopes[lt_sort]
ets = all_et_slopes[et_sort]
ltd = all_dens_lt[lt_sort]
etd = all_dens_et[et_sort]

max_lt = u.find_nearest(ltm, 10.5) + 1
max_et = u.find_nearest(etm, 10.5) + 1

slope_pars_lt, blah = curve_fit(f=power_law, xdata=ltm[0:max_lt], ydata=lts[0:max_lt])
slope_pars_et, blah = curve_fit(f=power_law, xdata=etm[0:max_et], ydata=ets[0:max_et])
dens_pars_lt, blah = curve_fit(f=power_law, xdata=ltm[0:max_lt], ydata=ltd[0:max_lt])
dens_pars_et, blah = curve_fit(f=power_law, xdata=etm[0:max_et], ydata=etd[0:max_et])


color_lt = (0.9,0.3,0.)
color_et = (0.,0.,0.52)

plt.figure(dpi=500)
plt.plot(logmass_c_et,slopes_c_et,linestyle='',marker='d',color=color_et)#,label='This Work')
plt.plot(logmass_c_lt,slopes_c_lt,linestyle='',marker='d',color=color_lt)#,label='This Work')
plt.plot(np.log10(gal_mass_lt),slopes_lt,linestyle='',marker='o',color=color_lt)#,label='Late-TypePechetti+20')
plt.plot(np.log10(gal_mass_et),slopes_et,linestyle='',marker='o',color=color_et)#,label='Pechetti+20')

mass_array = np.arange(8,10.5,0.1)
plt.plot(mass_array, slope_pars_lt[0]*mass_array+slope_pars_lt[1],
         linestyle='--',color=color_lt)
plt.plot(mass_array, slope_pars_et[0]*mass_array+slope_pars_et[1],
         linestyle='--',color=color_et)

plt.plot(2,2,marker='o',color='k',linestyle='',label='Pechetti+20')
plt.plot(2,2,marker='d',color='k',linestyle='',label='This Work')

plt.text(0.18, 0.2, 'Early-Type',color=color_et, transform=plt.gcf().transFigure, size=10)
plt.text(0.18, 0.16, 'Late-Type',color=color_lt, transform=plt.gcf().transFigure, size=10)

#plt.plot(8.0,-4.2,linestyle='',marker="$Early$",color='m',markersize=20)
#plt.plot(8.0,-4.3,linestyle='',marker="$Late$",color='c',markersize=20)

#plt.scatter(np.log10(gal_mass),slopes,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
#plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
#plt.plot(np.sort(np.log10(gal_mass)),0.42*np.log10(np.sort(gal_mass)/10**9)-2.39,linestyle='--',color='k')
plt.ylabel('$\gamma$')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(7.9,11.9)
plt.ylim(-4.8,0.2)
#plt.colorbar(label='0 = Early-type, 1 = Late-type')
#plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.legend()
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_and_Our_Sample/slopes_v_galmass_colored_by_morph.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

plt.figure(dpi=500)
plt.plot(logmass_c_et,np.log10(cen_dens_c_et),linestyle='',marker='d',color=color_et)#,label='This Work')
plt.plot(logmass_c_lt,np.log10(cen_dens_c_lt),linestyle='',marker='d',color=color_lt)#,label='This Work')
plt.plot(np.log10(gal_mass_lt),cen_dens_lt,linestyle='',marker='o',color=color_lt)#,label='Pechetti+20')
plt.plot(np.log10(gal_mass_et),cen_dens_et,linestyle='',marker='o',color=color_et)#,label='Pechetti+20')

plt.plot(mass_array, dens_pars_lt[0]*mass_array+dens_pars_lt[1],
         linestyle='--',color=color_lt)
plt.plot(mass_array, dens_pars_et[0]*mass_array+dens_pars_et[1],
         linestyle='--',color=color_et)
# =============================================================================
# plt.plot(np.sort(all_masses), dens_pars_lt[0]*np.sort(all_masses)+dens_pars_lt[1],
#          linestyle='--',color=color_lt)
# plt.plot(np.sort(all_masses), dens_pars_et[0]*np.sort(all_masses)+dens_pars_et[1],
#          linestyle='--',color=color_et)
# =============================================================================

plt.plot(2,-2,marker='o',color='k',linestyle='',label='Pechetti+20')
plt.plot(2,-2,marker='d',color='k',linestyle='',label='This Work')

plt.text(0.18, 0.2, 'Early-Type',color=color_et, transform=plt.gcf().transFigure, size=10)
plt.text(0.18, 0.16, 'Late-Type',color=color_lt, transform=plt.gcf().transFigure, size=10)

#plt.scatter(np.log10(gal_mass),cen_dens,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
#plt.scatter(logmass_corrected, np.log10(cen_dens_corrected),c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
#plt.plot(np.sort(np.log10(gal_mass)),0.61*np.log10(np.sort(gal_mass)/10**9)+2.78,linestyle='--',color='k')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(7.9,11.9)
plt.ylim(1.6,5.0)
#plt.colorbar(label='0 = Early-type, 1 = Late-type')
#plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.legend()
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_and_Our_Sample/densities_v_galmass_colored_by_morph.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)



# =============================================================================

#%%

# =============================================================================
# let's look at the distributions of scatter around our trendlines

lts_scatter = -(slope_pars_lt[0]*ltm[0:max_lt]+slope_pars_lt[1]) + lts[0:max_lt]
ets_scatter = -(slope_pars_et[0]*etm[0:max_et]+slope_pars_et[1]) + ets[0:max_et]

ltd_scatter = -(dens_pars_lt[0]*ltm[0:max_lt]+dens_pars_lt[1]) + ltd[0:max_lt]
etd_scatter = -(dens_pars_et[0]*etm[0:max_et]+dens_pars_et[1]) + etd[0:max_et]

width = 0.1
scat_bins = np.arange(-3.4,2+2*width, width) 

plt.figure(dpi=500)
plt.hist(lts_scatter,bins=scat_bins,alpha=0.7,label='Late')
plt.hist(ets_scatter,bins=scat_bins,alpha=0.7,label='Early')
plt.legend()
plt.xlabel('$\gamma_{measured}$ - $\gamma_{predicted}$')
plt.ylabel('Counts')

width = 0.05
scat_bins = np.arange(-1.3,1+2*width, width) 

plt.figure(dpi=500)
plt.hist(ltd_scatter,bins=scat_bins,alpha=0.7,label='Late')
plt.hist(etd_scatter,bins=scat_bins,alpha=0.7,label='Early')
plt.legend()
plt.xlabel('$\\rho_{5pc, measured}$ - $\\rho_{5pc, predicted}$')
plt.ylabel('Counts')

# =============================================================================



#%%
# =============================================================================
# let's plot the distribution of our galaxies

plt.figure(dpi=500)
w = 0.2
bs = np.arange(8,11.7+2*w,w)
plt.hist(all_masses, bins=bs, color='0.2',label='All Galaxies',alpha=0.3)
plt.hist(all_et_mass, bins=bs, histtype='step',color=color_et,label='Early-Types',alpha=0.7)
plt.hist(all_lt_mass, bins=bs, histtype='step',color=color_lt,label='Late-Types',alpha=0.7)
#plt.hist(all_etnsc_mass, bins=bs, histtype='step',color='c',label='Early-Types NSC', alpha=0.6)
         #linestyle='--')
#plt.hist(all_ltnsc_mass, bins=bs, histtype='step',color='r',label='Late-Types NSC', alpha=0.6)
         #linestyle='--')
plt.legend()
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.ylabel('# of Galaxies')


plt.figure(dpi=500)
w = 0.2
bs = np.arange(8,11.7+2*w,w)
plt.hist(all_masses, bins=bs, color='0.2',label='All Galaxies',alpha=0.3)
plt.hist(all_etnsc_mass, bins=bs, histtype='step',color=color_et,label='Early-Types w/ NSC', alpha=0.7)
         #linestyle='--')
plt.hist(all_ltnsc_mass, bins=bs, histtype='step',color=color_lt,label='Late-Types w/ NSC', alpha=0.7)
         #linestyle='--')
plt.legend()
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.ylabel('# of Galaxies')



# =============================================================================

#%%

# =============================================================================
# recreate Renuka's plots

fig1, ax1 = plt.subplots(nrows=1,ncols=1,dpi=500)

#ax1.scatter(np.log10(gal_mass),cen_dens,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
ax1.plot(np.log10(gal_mass),cen_dens,linestyle='',marker='.',color='b',label='Pechetti+20')
#plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
ax1.plot(np.sort(np.log10(gal_mass)),0.61*np.log10(np.sort(gal_mass)/10**9)+2.78,linestyle='--',color='r')
ax1.plot(np.log10(gal_mass[8]), cen_dens[8],marker='d',color='m')
ax1.set_ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$]')
ax1.set_xlim(8,11)
ax1.set_ylim(1.5,4.2)
ax1.set_aspect('equal', adjustable='box')


fig2, ax2 = plt.subplots(nrows=1,ncols=1,dpi=500)

#ax2.scatter(np.log10(gal_mass),slopes,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
ax2.plot(np.log10(gal_mass),slopes,linestyle='',marker='.',color='b',label='Pechetti+20')
#plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
ax2.plot(np.sort(np.log10(gal_mass)),0.42*np.log10(np.sort(gal_mass)/10**9)-2.39,linestyle='--',color='r')
ax2.plot(np.log10(gal_mass[8]), slopes[8],marker='d',color='m')
ax2.set_ylabel('Power-law ($\gamma$)')
ax2.set_xlabel('log(M$_*$ [M$_\odot$]')
ax2.set_xlim(8,11)
ax2.set_ylim(-3.3,-0.9)
ax2.set_aspect('equal', adjustable='box')


#%%

plt.figure(dpi=500)
plt.scatter(np.log10(gal_mass),cen_dens,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
#plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
plt.plot(np.log10(gal_mass),0.592*np.log10(gal_mass/10**9)+2.819,linestyle='--',color='r')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$]')
plt.xlim(8,11)
plt.ylim(1.5,4.2)

plt.figure(dpi=500)
plt.scatter(np.log10(gal_mass),slopes,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
#plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
plt.plot(np.log10(gal_mass),0.963*np.log10(gal_mass/10**9)-2.762,linestyle='--',color='r')
plt.ylabel('$\gamma$')
plt.xlabel('log(M$_*$ [M$_\odot$]')
plt.xlim(8,11)
plt.ylim(-3.1,-0.9)

# =============================================================================

#%%

# plot the slopes from our galaxies vs stone profile slopes
k = np.where(stone_slopes_c != -99)
p = np.where(NSC_comp_c == 0)
w = np.intersect1d(k,p)

stone_slopes_corrected = stone_slopes_c[k]
slopes_c_corrected = slopes_c[k]
stone_slopes_no_nsc = stone_slopes_c[w]
slopes_no_nsc = slopes_c[w]

plt.figure(dpi=500)
plt.plot(stone_slopes_corrected,slopes_c_corrected,linestyle='',marker='.',
         color='b',label='ALL')
plt.plot([-2,0],[-2,0],linestyle='--',color='r')
plt.xlabel('$\gamma$ from Stone&Metzger 2016')
plt.ylabel('$\gamma$ from This Work')
plt.ylim(-3.6,0.1)
plt.legend()

plt.figure(dpi=500)
plt.plot(stone_slopes_no_nsc,slopes_no_nsc,linestyle='',marker='.',
         color='b', label='No NSC')
plt.plot([-2,0],[-2,0],linestyle='--',color='r')
plt.xlabel('$\gamma$ from Stone&Metzger 2016')
plt.ylabel('$\gamma$ from This Work')
plt.legend()

# plot the densities from our galaxies vs stone profile densities
stone_dens_at_5pc_c_corrected = stone_dens_at_5pc_c[k]
cen_dens_c_corrected = cen_dens_c[k]
stone_dens_at_5pc_no_nsc = stone_dens_at_5pc_c[w]
cen_dens_no_nsc = cen_dens_c[w]

plt.figure(dpi=500)
plt.plot(np.log10(stone_dens_at_5pc_c_corrected),np.log10(cen_dens_c_corrected),
         linestyle='',marker='.',color='b',label='ALL')
plt.plot([2,6],[2,6],linestyle='--',color='r')
plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from Stone&Metzger 2016')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from This Work')
plt.legend()

plt.figure(dpi=500)
plt.plot(np.log10(stone_dens_at_5pc_no_nsc),np.log10(cen_dens_no_nsc),
         linestyle='',marker='.',color='b', label='No NSC')
plt.plot([2,6],[2,6],linestyle='--',color='r')
plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from Stone&Metzger 2016')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from This Work')
plt.legend()



#%%
# =============================================================================
# plot the galaxy distribution







# =============================================================================




#%%
# =============================================================================
# Plot galaxy mass vs NSC Reff colored by density
# =============================================================================
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(log_dens_5pc != -1)
# x = NSC_reff[a]
# y = np.log10(gal_mass[a])
# z = log_dens_5pc[a]
# 
# plt.figure(dpi=500)
# plt.scatter(x,y,c=z,cmap='viridis')
# plt.xlabel('R$_{eff,NSC}$ [pc]')
# plt.ylabel('log(M$_*$ [M$_\odot$])')
# plt.colorbar(label='log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_gal_v_R_NSC_color_by_density.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# 
# # =============================================================================
# # Plot galaxy mass vs NSC Reff colored by density slope (gamma)
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(gamma != 1)
# x = NSC_reff[a]
# y = np.log10(gal_mass[a])
# z = gamma[a]
# 
# plt.figure(dpi=500)
# plt.scatter(x,y,c=z,cmap='viridis')
# plt.xlabel('R$_{eff,NSC}$ [pc]')
# plt.ylabel('log(M$_*$ [M$_\odot$])')
# plt.colorbar(label='$\gamma$')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_gal_v_R_NSC_color_by_slope.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# 
# # =============================================================================
# # Plot NSC_mass v Reff colored by density
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(log_dens_5pc != -1)
# b = np.where(NSC_mass != -1)
# k = np.intersect1d(a[0],b[0])
# 
# x = NSC_reff[k]
# y = np.log10(NSC_mass[k])
# z = log_dens_5pc[k]
# 
# plt.figure(dpi=500)
# plt.scatter(x,y,c=z,cmap='viridis')
# plt.xlabel('R$_{eff,NSC}$ [pc]')
# plt.ylabel('log(M$_{NSC}$ [M$_\odot$])')
# plt.colorbar(label='log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_v_R_NSC_color_by_density.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# 
# # =============================================================================
# # Plot NSC_mass v Reff colored by gamma
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(gamma != 1)
# b = np.where(NSC_mass != -1)
# k = np.intersect1d(a[0],b[0])
# 
# x = NSC_reff[k]
# y = np.log10(NSC_mass[k])
# z = gamma[k]
# 
# plt.figure(dpi=500)
# plt.scatter(x,y,c=z,cmap='viridis')
# plt.xlabel('R$_{eff,NSC}$ [pc]')
# plt.ylabel('log(M$_{NSC}$ [M$_\odot$])')
# plt.colorbar(label='$\gamma$')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_v_R_NSC_color_by_slope.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# 
# 
# # =============================================================================
# # Plot NSC_mass/r^3 vs density
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(log_dens_5pc != -1)
# b = np.where(NSC_mass != -1)
# k = np.intersect1d(a[0],b[0])
# 
# x = NSC_reff[k]
# y = NSC_mass[k]
# z = log_dens_5pc[k]
# 
# plt.figure(dpi=500)
# 
# plt.plot(z, np.log10(y/x**3), linestyle='', marker='o', color='b')
# #plt.scatter(x,y,c=z,cmap='viridis')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# plt.ylabel('log(M$_{NSC}$/R$_{eff,NSC}$$^3$ [M$_\odot$/pc$^3$])')
# #plt.colorbar(label='log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_over_R_NSC_cubed_v_density.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# 
# 
# # =============================================================================
# # Plot NSC_mass/r^3 vs density
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(gamma!= 1)
# b = np.where(NSC_mass != -1)
# k = np.intersect1d(a[0],b[0])
# 
# x = NSC_reff[k]
# y = NSC_mass[k]
# z = gamma[k]
# 
# plt.figure(dpi=500)
# 
# plt.plot(z, np.log10(y/x**3), linestyle='', marker='o', color='b')
# #plt.scatter(x,y,c=z,cmap='viridis')
# plt.xlabel('$\gamma$')
# plt.ylabel('log(M$_{NSC}$/R$_{eff,NSC}$$^3$ [M$_\odot$/pc$^3$])')
# #plt.colorbar(label='log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_over_R_NSC_cubed_v_slope.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# 
# 
# # =============================================================================
# # Plot NSC_mass/r^2 vs density
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(log_dens_5pc != -1)
# b = np.where(NSC_mass != -1)
# k = np.intersect1d(a[0],b[0])
# 
# x = NSC_reff[k]
# y = NSC_mass[k]
# z = log_dens_5pc[k]
# 
# plt.figure(dpi=500)
# 
# plt.plot(z, np.log10(y/x**2), linestyle='', marker='o', color='b')
# #plt.scatter(x,y,c=z,cmap='viridis')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# plt.ylabel('log(M$_{NSC}$/R$_{eff,NSC}$$^2$ [M$_\odot$/pc$^3$])')
# #plt.colorbar(label='log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_over_R_NSC_squared_v_density.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# 
# 
# # =============================================================================
# # Plot NSC_mass/r^2 vs density
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(gamma!= 1)
# b = np.where(NSC_mass != -1)
# k = np.intersect1d(a[0],b[0])
# 
# x = NSC_reff[k]
# y = NSC_mass[k]
# z = gamma[k]
# 
# plt.figure(dpi=500)
# 
# plt.plot(z, np.log10(y/x**2), linestyle='', marker='o', color='b')
# #plt.scatter(x,y,c=z,cmap='viridis')
# plt.xlabel('$\gamma$')
# plt.ylabel('log(M$_{NSC}$/R$_{eff,NSC}$$^2$ [M$_\odot$/pc$^3$])')
# #plt.colorbar(label='log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_over_R_NSC_csquared_v_slope.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# 
# #%%
# # =============================================================================
# # Plot NSC_mass/r^2 vs density
# 
# # modify arrays to only include gals with the appropriate measurments
# a = np.where(gamma != 1)
# b = np.where(NSC_mass != -1)
# #c = np.where(log_dens_5pc != -1)
# k = np.intersect1d(a[0],b[0])
# 
# 
# x = NSC_reff[k]
# y = NSC_mass[k]
# z = gamma[k]
# w = log_dens_5pc[k]
# 
# plt.figure(dpi=500)
# #plt.plot(z, np.log10(y/x**2), linestyle='', marker='o', color='b')
# plt.scatter(w,np.log10(y),c=x,cmap='viridis')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# plt.ylabel('log(M$_{NSC}$ [M$_\odot$])')
# plt.colorbar(label='R$_{eff,NSC}$ [pc]')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_v_density_color_by_reff.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# plt.figure(dpi=500)
# #plt.plot(z, np.log10(y/x**2), linestyle='', marker='o', color='b')
# plt.scatter(z,np.log10(y),c=x,cmap='viridis')
# plt.xlabel('$\gamma$')
# plt.ylabel('log(M$_{NSC}$ [M$_\odot$])')
# plt.colorbar(label='R$_{eff,NSC}$ [pc]')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_v_slope_color_by_reff.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# plt.figure(dpi=500)
# #plt.plot(z, np.log10(y/x**2), linestyle='', marker='o', color='b')
# plt.scatter(w, np.log10(y/x**2),c=x,cmap='viridis')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# plt.ylabel('log(M$_{NSC}$/R$_{eff,NSC}$$^2$ [M$_\odot$/pc$^3$])')
# plt.colorbar(label='R$_{eff,NSC}$ [pc]')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_over_R_NSC_csquared_v_density_color_by_reff.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# #plt.plot(z, np.log10(y/x**2), linestyle='', marker='o', color='b')
# plt.scatter(z, np.log10(y/x**2),c=x,cmap='viridis')
# plt.xlabel('$\gamma$')
# plt.ylabel('log(M$_{NSC}$/R$_{eff,NSC}$$^2$ [M$_\odot$/pc$^3$])')
# plt.colorbar(label='R$_{eff,NSC}$ [pc]')
# 
# plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Pechetti_cross_sections_pngs/M_NSC_over_R_NSC_csquared_v_slope_color_by_reff.png',
#             bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# # =============================================================================
# 
# =============================================================================
