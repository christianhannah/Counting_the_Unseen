#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 16:17:03 2022

@author: christian
"""

#from mgefit.mge_fit_1d import mge_fit_1d
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import mge1d_util as u
import sys

import astropy.units as uni
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery


import warnings
warnings.filterwarnings("ignore")

#%%

def get_density_profile(gal_name):

    # check to see if a galaxy is in the Stone & Metzger (2016) paper sample
# =============================================================================
#     if u.check_for_galaxy_in_stone(desc+number) == False:
#         print('Galaxy not in Stone&Metzger(2016).')
# =============================================================================
    all_MLs = []
    all_ML_types = []

#%%

    # specify name of galaxy for fit
    #gal_name = desc+' '+number
    #gal_name_nospace = u.format_gal_name(gal_name)
    gal_name_nospace = gal_name

    # ========================= READ IN DATA ======================================

    # read in the MGE parameters from Pechetti et al. (2017) 
    mge_filename = '../Data_Sets/Pechetti17_data/pechetti_2017_data.fit'
    hdul = fits.open(mge_filename)
    data1 = hdul[1].data
    
    
    # store the data for specific galaxy
    gal_inds = np.where(data1['Name']==gal_name_nospace)[0]
    L_tots = 10**data1['log_I_'][gal_inds]
    sigmas = 10**data1['log_Sig_'][gal_inds]
    ax_ratios = data1['q'][gal_inds]
    filt = np.unique(data1['Filt'][gal_inds])[0]
    
    # ensure the galaxy name returned data
    if len(gal_inds) == 0:
        print('*** ERROR: No data found for '+gal_name+'. ***')
        return

# =============================================================================


###############################################################################


# ===================== COMPUTE 3D DENSITY PROFILE ============================

    # store distance to convert sigma from arcsec to pc
    distance = u.get_galaxy_distance(gal_name_nospace)*10**6
    #pdb.set_trace()
    # conversion factor between radians and arcsec
    arcsec_to_radian = 4.8481e-6

    # store RA and DEC for cross matching with Sloan atlas
    RA, DEC = u.get_galaxy_RA_and_Dec(gal_name_nospace)
    #pdb.set_trace()
    
    # correct the peak surface brightnesses from each gaussian for extinction
    coords = SkyCoord(ra=RA*uni.degree, dec=DEC*uni.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    
    wfc3_names = np.array(['NGC3522','PGC028887','PGC050395'])
    
    # correct the A_v value for the appropriate bandpass 
    # (values from Table 6 in Schlafly & Finkbeiner (2011))
    if filt == 'F475W': #ACS
        A_v = E_BV*3.268
    elif filt == 'F606W': #WFPC2
        A_v = E_BV*2.415
    elif filt == 'F702W': #WFPC2
        A_v = E_BV*1.948
    elif filt == 'F814W':
        # decide between WFPC2 and WFC3 correction
        if gal_name in wfc3_names: #WFC3
            A_v = E_BV*1.536
        else: #WFPC2
            A_v = E_BV*1.549
        
    L_tots = L_tots - A_v
    
    # store v-band magnitude for galaxy
    vmag = u.get_galaxy_vmag(gal_name_nospace)
    
    # get M/L ratio for the galaxy
    M_L = u.get_renuka_ML(gal_name_nospace)
    print('M/L = ', M_L)
    ML_type = 2
    
    # store all M/L's possible for comparison
    all_MLs.append(M_L)
    all_ML_types.append(ML_type)
    
    L_v = 10**(0.4*(-vmag + 4.83)) # L_sol_v
    
    # get M/L_i from the color-M/L relation of Taylor+2011
    #ML_type = 5
    M_L_temp, ML_type_temp, all_NSA_MLs, all_NSA_ML_types, gicolor = u.get_nsa_M_L(gal_name, RA, DEC, distance)

    # store all possible M/L values for comparison
    if M_L != -999:
        for i in range(len(all_NSA_MLs)):
            all_ML_types.append(all_NSA_ML_types[i])
            all_MLs.append(all_NSA_MLs[i])
    M_L_temp = 4.9*(L_v/10**10)**0.18  
    ML_type_temp = 1
    all_ML_types.append(ML_type_temp)
    all_MLs.append(M_L_temp) 
    
    
    # construct a radial sampling array for the density profiles to be evaluated
    num_rad = 100
    rad = np.geomspace(0.0456/2., 10, num_rad)
    radii = (rad*arcsec_to_radian)*distance # in pc  
    
    # define array to store radial densities
    densities = np.zeros_like(radii)

    heights = []
    sigs = []

    # compute 3-d density assuming spherical symmetry (so, as f(r))
    for j in range(len(radii)):
        for i in range(len(L_tots)):
            
            # store total luminosity and width of each gaussian component
            L_tot_i = L_tots[i] # in solar luminosities
            #height = L_tot_i/(np.sqrt(2*np.pi)*sigmas[i]) # height of 1-d Gaussian
            height = L_tots[i]
            sig_i = sigmas[i] # in arcsec
            
            # convert sigma from arcsec to pc using distance
            sig = (sig_i*arcsec_to_radian)*distance # in pc
            #L_tot = (height)*2*np.pi*(sig)**2 # total area under 2-d gaussian
            L_tot = L_tots[i]*2*np.pi*(sig)**2#*ax_ratios[i] # total area under 2-d gaussian
            if j == 0:
                heights.append(height)
                sigs.append(sig_i)
            
            # convert total luminosity to total mass using M/L ratio
            M_tot = L_tot*M_L
            
            # compute desity contribution from one gaussian component
            dens_comp = (M_tot/((2*np.pi)**(3/2)*sig**3))*np.exp(-(radii[j]**2/(2*sig**2))) # M_sol/pc^3
            
            # add density from each gaussian to total density
            densities[j] += dens_comp

    #pdb.set_trace()
#%%
# =============================================================================
#     # plot density profile
#     plt.figure(dpi=500)
#     plt.plot(np.log10(radii), np.log10(densities), label='This Work')
#     plt.xlabel('log(Radius [pc])')
#     plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
#     plt.xlim(np.min(np.log10(radii))-0.2,np.max(np.log10(radii))+0.2)
#     x0,xmax = plt.xlim()
#     y0,ymax = plt.ylim()
#     width = np.abs(xmax-x0)
#     height = np.abs(ymax-y0)
#     plt.text(x0+width*.05, y0+height*0.05, gal_name+', M$_v$:{:.2f}'.format(vmag)+
#              ', '+filt) 
# =============================================================================

    #plt.savefig('../Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_density_profile.pdf',
    #            bbox_inches='tight', pad_inches=0.1, dpi=500)
    #plt.savefig('../Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_density_profile.png',
    #            bbox_inches='tight', pad_inches=0.1, dpi=500)

# =============================================================================
#%%

    def gauss(height, sig, r):
        return height*np.exp((-0.5)*r**2/sig**2)

    gausses = [] 
    for i in range(len(sigs)):
        gausses.append(gauss(heights[i],sigs[i],rad))

    summed_gauss = np.zeros_like(rad)
    for i in range(len(rad)):
        for j in range(len(sigs)):
            summed_gauss[i] += gausses[j][i]

# =============================================================================
#     plt.figure(dpi=500)
#     for i in range(len(sigs)):
#         plt.plot(rad,gausses[i],linestyle='--')
# 
#     plt.plot(rad, summed_gauss, color='orange',linewidth=2.0)#, alpha=0.7)    
#     plt.ylim(np.min(summed_gauss)-0.1*np.min(summed_gauss),
#              np.max(summed_gauss)+0.1*np.max(summed_gauss))
#     plt.yscale('log')
#     plt.xscale('log')
#     plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
#     plt.xlabel('Radius ["]')
# =============================================================================


    #pdb.set_trace()
    return densities, radii, vmag, M_L, ML_type, distance, RA, DEC, filt, heights, \
        sigs, rad, summed_gauss, np.array(all_MLs), np.array(all_ML_types), gicolor, summed_gauss


#get_density_profile('PGC028887')


