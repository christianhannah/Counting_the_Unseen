#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 11:13:04 2022

@author: christian
"""

from mgefit.mge_fit_1d import mge_fit_1d
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import mge1d_util as u
import sys
import csv

import warnings
warnings.filterwarnings("ignore")

#%%

def get_density_profile(gal_name):


#%%

    # ========================= READ IN DATA ======================================

    all_names = []
    all_sersic_inds = []
    all_filts = []
    all_notes = []
    all_r_effs = []
    all_M_Ls = []
    
    table_b2_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Hoyer22_data/tab_publication.ascii.csv'
    # fields: ['galaxy', 'notes', 'filter', 'pa', 'l_pa', 'u_pa', 'ell', 'l_ell', 'u_ell', 'n', 'l_n', 'u_n', 'reff', 'u_reff', 'm', 'l_m', 'u_m', 'v', 'u_v', 'i', 'u_i', 'mlr', 'u_mlr', 'logm', 'u_logm']
    with open(table_b2_filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0][0] != '#' and row[0][0] != 'g': # marks the start of data row
                all_names.append(row[0])
                all_sersic_inds.append(row[9])
                all_filts.append(row[2])
                all_notes.append(row[1])
                all_r_effs.append(row[12])
                all_M_Ls.append(row[21])
    
    all_names = np.array(all_names)
    all_sersic_inds = np.array(all_sersic_inds).astype(float)
    all_filts = np.array(all_filts)
    all_notes = np.array(all_notes)
    all_r_effs = np.array(all_r_effs).astype(float)
    all_M_Ls = np.array(all_M_Ls).astype(float)
            

    # store the data for specific galaxy
    gal_inds = np.where(all_names == gal_name)[0]
    sersic_inds = all_sersic_inds[gal_inds]
    r_effs = all_r_effs[gal_inds]
    M_L = all_M_Ls[gal_inds[0]]
    notes = all_notes[gal_inds]
    filts = all_filts[gal_inds]
    
    
    # ensure the galaxy name returned data
    if len(gal_inds) == 0:
        print('*** ERROR: No data found for '+gal_name+'. ***')
        return

# =============================================================================


###############################################################################


# ======================== REPARE DATA FOR FIT ================================

    # store distance to convert sigma from arcsec to pc
    distance = gal_data['Dist'][0] * 10**6 # in pc
    # conversion factor between radians and arcsec
    arcsec_to_radian = 4.8481e-6

    # convert radii from arcsec to pc
    init_rad = (SB_data['MajAxis']*arcsec_to_radian)*distance # in pc    
    # find index corresponding to maximum radius desired for fit
    max_rad_ind = u.find_nearest(init_rad, max_rad)

    # ensure the radii are logarithmically spaced 
    num_rad = 100
    #rad = np.geomspace(np.min(SB_data['MajAxis']), np.max(SB_data['MajAxis']), num_rad)
    rad = np.geomspace(np.min(SB_data['MajAxis']), SB_data['MajAxis'][max_rad_ind], num_rad)
    SBs_i = 10**np.interp(np.log10(rad), np.log10(SB_data['MajAxis']), np.log10(SB_data['SuBr']))
    
    # can switch interpolaton to cubic spline
    #f = interpolate.CubicSpline(np.log10(SB_data['MajAxis']), np.log10(SB_data['SuBr']))
    #SBs_i = 10**f(np.log10(rad))

# =============================================================================
#     # plot interpolated data with the original
#     plt.figure(dpi=500)
#     #plt.plot(rad,SBs, marker='.')
#     plt.plot(rad,SBs_i, marker='.', label='Interpolated')
#     plt.plot(SB_data['MajAxis'], SB_data['SuBr'], marker='+', alpha=0.7, label='Original')
#     plt.ylim(np.max(SB_data['SuBr'])+0.2, np.min(SB_data['SuBr'])-0.2)
#     plt.ylabel('$\mu$ [mag/arcsec$^2$]')
#     plt.xlabel('Radius [arcsec]')
#     plt.text(np.max(rad)-7.5, np.min(SBs_i)+0.5, gal_name+' ('+str(SB_data['Filt'][0])+')')
#     plt.legend()
#     #plt.show()
#     plt.clf()
# =============================================================================
    
    
    # convert SBs from mag/arcsec^2 to L_sol/pc^2
    if SB_data['Filt'][0] == 'F555W':
        M_sol = 4.84 # Vega system
    elif SB_data['Filt'][0] == 'F814W':
        M_sol = 4.12 # Vega system
    elif SB_data['Filt'][0] == 'F606W':
        M_sol = 4.62 # Vega system
    SBs = 10**((SBs_i - M_sol - 21.572)/(-2.5))

    SBs_true = 10**((SB_data['SuBr'] - M_sol - 21.572)/(-2.5))

# =============================================================================
#     # plot converted SB profile
#     plt.figure(dpi=500)
#     plt.plot(rad, SBs)
#     plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
#     plt.xlabel('Radius [arcsec]')
#     plt.text(np.max(rad)-5.5, np.max(SBs)-0.5, gal_name+' ('+str(SB_data['Filt'][0])+')')
#     #plt.show()
#     plt.clf()
# =============================================================================

# =============================================================================


###############################################################################


# ======================== DEPROJECT SB PROFILE================================

    mge = mge_fit_1d(rad, SBs, ngauss=20, plot=True)
    plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_mge_fit.pdf',
                bbox_inches='tight', pad_inches=0.1, dpi=500)
    plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_mge_fit.png',
                bbox_inches='tight', pad_inches=0.1, dpi=500)
    plt.clf()

# =============================================================================
#     # store distance to convert sigma from arcsec to pc
#     distance = gal_data['Dist'][0] * 10**6 # in pc
#     # conversion factor between radians and arcsec
#     arcsec_to_radian = 4.8481e-6
# =============================================================================

    # store RA and DEC for cross matching with Sloan atlas
    RA = gal_data_w_reff['_RA'][0]
    DEC = gal_data_w_reff['_DE'][0]
    
    # compute M/L ratio for the galaxy
    #M_L = 11 # temporary
    vmag = gal_data['VMAG'][0]
    L_v = 10**(0.4*(-vmag + 4.83)) # L_sol_v
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
    M_L_1 = 4.9*(L_v/10**10)**0.18
    print('M/L = ', M_L)
    
    # convert radii from arcsec to pc
    radii = (rad*arcsec_to_radian)*distance # in pc
    
    # compute the total luminosity to compare with total luminosity computed with mge gaussians
    L_tot_SB = u.integrate_SB(radii, SBs)
    
    # define array to store radial densities
    densities = np.zeros_like(rad)
    
    # store the total luminosity from mge gaussians 
    L_tot_mge = 0

    heights = []
    sigs = []
    tots = []
    L_tots = []

    # compute 3-d density assuming spherical symmetry (so, as f(r))
    for j in range(len(rad)):
        for i in range(len(mge.sol[0])):
            # store total luminosity and width of each gaussian component
            L_tot_i = mge.sol[0][i] # in solar luminosities
            height = L_tot_i/(np.sqrt(2*np.pi)*mge.sol[1][i]) #height of 1-d Gaussian
            sig_i = mge.sol[1][i] # in arcsec
            # convert sigma from arcsec to pc using distance
            sig = (sig_i*arcsec_to_radian)*distance # in pc
            L_tot = (height)*2*np.pi*(sig)**2 # total area under 2-d gaussian
            #pdb.set_trace()
            if j == 0:
                tots.append(L_tot_i)
                heights.append(height)
                sigs.append(sig_i)
                L_tots.append(L_tot)
                L_tot_mge += L_tot
            # convert total luminosity to total mass using M/L ratio
            M_tot = L_tot*M_L
            # compute desity contribution from one gaussian component
            dens_comp = (M_tot/((2*np.pi)**(3/2)*sig**3))*np.exp(-(radii[j]**2/(2*sig**2))) # M_sol/pc^3
            # add density from each gaussian to total density
            #print(dens_comp, radii[j])
            densities[j] += dens_comp
        #pdb.set_trace()

    print('L_tot from SB data =  ', L_tot_SB, ' L_sol')
    print('L_tot from MGE data = ', L_tot_mge, ' L_sol')
    print('Ratio (SB/MGE) = ',L_tot_SB/L_tot_mge)

#%%  
    # read data from Stone & Metzger (2016) who deprojected densities using nuker law 
    # SB profiles (essentially double power-laws) and assuming spherical symmetry
    stone_data =  u.get_stone_data(gal_name_nospace)   
    if stone_data != 0:
        stone_rad = stone_data[1]
        stone_dens = stone_data[2]
    else:
        stone_rad = 0
        stone_dens = 0

#%%
# =============================================================================
#     # plot density profile
#     plt.figure(dpi=500)
#     plt.plot(np.log10(radii), np.log10(densities), label='This Work')
#     if stone_data != 0:
#         plt.plot(np.log10(stone_rad), np.log10(stone_dens), alpha=0.7, label='Stone&Metzger2016')
#     plt.xlabel('log(Radius [pc])')
#     plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
#     if stone_data != 0:
#         plt.legend()
#     plt.xlim(np.min(np.log10(radii))-0.2,np.max(np.log10(radii))+0.2)
#     if stone_data != 0:
#         plt.ylim(min([np.min(np.log10(densities)),np.min(np.log10(stone_dens))])-0.2,
#                  max([np.max(np.log10(densities)),np.max(np.log10(stone_dens))])-1)
#     else:
#         plt.ylim(np.min(np.log10(densities))-0.2,
#                  np.max(np.log10(densities))+0.2)
#     x0,xmax = plt.xlim()
#     y0,ymax = plt.ylim()
#     width = np.abs(xmax-x0)
#     height = np.abs(ymax-y0)
#     plt.text(x0+width*.05, y0+height*0.05, gal_name+', M$_v$:'+str(gal_data['VMAG'][0])+
#              ', '+SB_data['Filt'][0]) 
# 
#     plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_density_profile.pdf',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#     plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_density_profile.png',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================

# =============================================================================
#%%

    def gauss(height, sig, r):
        return height*np.exp((-0.5)*r**2/sig**2)

    gausses = [] 
    for i in range(len(mge.sol[0])):
        gausses.append(gauss(heights[i],sigs[i],rad))

    summed_gauss = np.zeros_like(rad)
    for i in range(len(rad)):
        for j in range(len(mge.sol[0])):
            summed_gauss[i] += gausses[j][i]

# =============================================================================
#     plt.figure(dpi=500)
#     plt.plot(rad, SBs, color='b', linestyle='', marker='o')
#     for i in range(len(mge.sol[0])):
#         plt.plot(rad,gausses[i])
# 
#     plt.plot(rad, summed_gauss, color='orange',)#, alpha=0.7)    
# 
#     plt.yscale('log')
#     plt.xscale('log')
#     plt.ylim(min(SBs)-20,max(SBs)+1000)
#     plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
#     plt.xlabel('Radius ["]')
# 
# 
#     plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/reconstructed_MGEs/'+gal_name_nospace+'_reconstructed_MGE_plot.pdf',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#     plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/reconstructed_MGEs/'+gal_name_nospace+'_reconstructed_MGE_plot.png',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================

    L_tot_gauss = u.integrate_SB(radii, summed_gauss)

    print()
    print('L_tot from SB data =     ', L_tot_SB, ' L_sol')
    print('L_tot from summed MGE  = ', L_tot_gauss, ' L_sol')
    
    return densities, radii, vmag, M_L, ML_type, distance, RA, DEC, SB_data['SuBr'],\
        SB_data['MajAxis'], SBs_i, rad, heights, sigs, SBs, SBs_true




get_density_profile('blah')





