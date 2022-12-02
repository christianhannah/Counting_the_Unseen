#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program holds multiple functions for use in deprojecting density profiles 
from 1-D SB profiles of galaxies using mge_fit_1d.py

Created on Wed Jul 28 14:14:57 2021

@author: christian
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
import pdb
from astropy.stats import median_absolute_deviation as mad
from astroquery.ned import Ned

# =============================================================================
# function to get the log stellar mass for a galaxy from david's table
def get_galaxy_logmass(gal_name): # use no space gal name 
    # read in fits table with galaxy masses from David
    filename = '/Users/christian/OneDrive/Desktop/TDE Code/Galaxy Masses/sample_3.fits'
    # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
    hdul = fits.open(filename)  # open a FITS file
    data = hdul[1].data 
    hdul.close()
    
    # edit some names to match the ID in David's table
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1545':
        gal_name = 'IC3509'
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1627':
        gal_name = 'PGC041180'
    if gal_name == 'VCC1199':
        gal_name = 'PGC041959'
    
        
    
    ind = np.where(data['objname'] == gal_name)
    if ind[0].size == 0:
        print(gal_name+" not found in David's table")
        return 0
    else:
        return data['logmass'][ind][0]#, data['logmass_src'][ind]
# =============================================================================

# =============================================================================
# function to get the distance for a galaxy from david's table
def get_galaxy_distance(gal_name): # use no space gal name 
    # read in fits table with galaxy masses from David
    filename = '/Users/christian/OneDrive/Desktop/TDE Code/Galaxy Masses/sample_3.fits'
    # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
    hdul = fits.open(filename)  # open a FITS file
    data = hdul[1].data 
    hdul.close()
    
    # edit some names to match the ID in David's table
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1545':
        gal_name = 'IC3509'
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1627':
        gal_name = 'PGC041180'
    if gal_name == 'VCC1199':
        gal_name = 'PGC041959'
    if gal_name == 'BTS76':
        gal_name = 'PGC2832100'
    if gal_name == 'DDO084':
        gal_name = 'UGC05829'
    if gal_name == 'DDO133':
        gal_name = 'UGC07698'
    if gal_name == 'LEG09':
        gal_name = 'PGC083321'
    
    
    
    
    ind = np.where(data['objname'] == gal_name)
    if ind[0].size == 0:
        print(gal_name+" not found in David's table")
        return 0
    else:
        return data['bestdist'][ind][0]
# =============================================================================

# =============================================================================
# function to get the v-band mag for a galaxy from David's table
def get_galaxy_vmag(gal_name): # use no space gal name 
    # read in fits table with galaxy masses from David
    filename = '/Users/christian/OneDrive/Desktop/TDE Code/Galaxy Masses/sample_3.fits'
    # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
    hdul = fits.open(filename)  # open a FITS file
    data = hdul[1].data 
    hdul.close()
    
    # edit some names to match the ID in David's table
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1545':
        gal_name = 'IC3509'
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1627':
        gal_name = 'PGC041180'
    if gal_name == 'VCC1199':
        gal_name = 'PGC041959'
    
    ind = np.where(data['objname'] == gal_name)
    if ind[0].size == 0:
        print(gal_name+" not found in David's table")
        return 0
    else:
        return data['VMag'][ind][0]
# =============================================================================

# =============================================================================
# function to get the v-band mag for a galaxy from David's table
def get_galaxy_RA_and_Dec(gal_name): # use no space gal name 
    # read in fits table with galaxy masses from David
    filename = '/Users/christian/OneDrive/Desktop/TDE Code/Galaxy Masses/sample_3.fits'
    # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
    hdul = fits.open(filename)  # open a FITS file
    data = hdul[1].data 
    hdul.close()
    
    # edit some names to match the ID in David's table
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1545':
        gal_name = 'IC3509'
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1440':
        gal_name = 'IC0798'
    if gal_name == 'VCC1627':
        gal_name = 'PGC041180'
    if gal_name == 'VCC1199':
        gal_name = 'PGC041959'
    
    ind = np.where(data['objname'] == gal_name)
    if ind[0].size == 0:
        print(gal_name+" not found in David's table")
        return 0
    else:
        return data['RA'][ind][0], data['Dec'][ind][0]
# =============================================================================


# =============================================================================
# function to get the M/L for a galaxy from Renuka_ML.dat file
def get_renuka_ML(gal_name): # use no space gal name 
    
    file_path = '/Users/christian/OneDrive/Desktop/TDE Code/Pechetti17_data/Renuka_ML.dat'
    f = open(file_path, 'r')
    all_lines = f.readlines()
    f.close()
    
    names = []
    MLs = []
    for i in range(len(all_lines)-1):
        split = all_lines[i+1].split()
        names.append(split[0])
        MLs.append(float(split[7]))
        
    names = np.array(names)
    MLs = np.array(MLs)
    gal_ind = np.where(names == gal_name)[0]
    
    if len(gal_ind) == 0:
        print("*** ERROR: "+gal_name+" not found in Renuka_ML.dat ***")
        return 0
    else:
        return MLs[gal_ind]

# =============================================================================


# =============================================================================
# function to check that a galaxy exists in the deprojected Stone & Metzger 
# (2016) sample
def check_for_galaxy_in_stone(gal_name): # use no space gal name 

    file_path = '/Users/christian/OneDrive/Desktop/TDE Code/Stone_Metzger_data/RenukaAnil.dat'
    f = open(file_path, 'r')
    all_lines = f.readlines()
    f.close()
    
    names = []
    for i in range(len(all_lines)):
        split = all_lines[i].split()
        names.append(split[0])    
    names = np.array(names)
    #pdb.set_trace()
    names_f = names[np.where(names == gal_name)]  
    
    if len(names_f) != 0:
        return True
    else:
        return False
# =============================================================================

# =============================================================================
# function to format the galaxy name similar to stone (no spaces, etc.)
def format_gal_name(gal_name): 
    
    final_name = ''
    parts = gal_name.split()
    if len(parts) > 2:
        print('Update the name format function for stone data.')
    if parts[0] == 'NGC':
        if len(parts[1]) == 3:
            final_name = parts[0]+'0'+parts[1]
        else:
            final_name = parts[0]+parts[1]
    elif len(parts) == 1:
        final_name = parts[0]
    else:
        final_name = parts[0]+parts[1]
    
    return final_name
# =============================================================================

# =============================================================================
# function to read in 3-D radius and 3-D densities from Stone & Metzger (2016)
def get_stone_data(gal_name): # use no space gal name 
    
    file_path = '/Users/christian/OneDrive/Desktop/TDE Code/Stone_Metzger_data/RenukaAnil.dat'
    f = open(file_path, 'r')
    all_lines = f.readlines()
    f.close()
    
    names = []
    radii = np.zeros(len(all_lines))
    densities = np.zeros(len(all_lines))

    for i in range(len(all_lines)):
        split = all_lines[i].split()
        names.append(split[0])
        radii[i] = float(split[1])
        densities[i] = float(split[2])
        
    names = np.array(names)
    names_f = names[np.where(names == gal_name)]
    radii_f = radii[np.where(names == gal_name)]
    densities_f = densities[np.where(names == gal_name)]
    
    if len(names_f) == 0:
        return 0
    else:
        return [names_f, radii_f, densities_f]
# =============================================================================

# =============================================================================
def integrate_SB(r,y):
    rsum = 0
    rad = 0
    dr_pre = 0
    for i in range(len(r)-1):
        width = np.abs(r[i]-r[i+1])
        rad += dr_pre+width/2
        #print(rad)
        dr_pre = width/2
        height = (min(y[i],y[i+1])+(np.abs(y[i]-y[i+1])/2))
        rsum += height*rad*width
    return rsum*2*np.pi

# =============================================================================

# =============================================================================
# define function to find value in an array nearest to a supplied value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# =============================================================================

# =============================================================================
# function to compute and return the power-law slope of the density profile at 
# the center 
def power_law(x,a,b):
    return a*x+b
def get_density_and_slope_simple(x,y):
    max_ind = 2
    x1 = x[0:max_ind]
    y1 = y[0:max_ind]
    pars, cov = curve_fit(f=power_law, xdata=x1, ydata=y1)
    par = np.polyfit(x1,y1,deg=1)
    #pdb.set_trace()
    #plt.figure(dpi=500)
    #plt.plot(x1,y1)
    #plt.plot(x1, power_law(x1,pars[0],pars[1]))
    ind_5pc = find_nearest(x, np.log10(5))
    inner_dens = y[ind_5pc]
    #pdb.set_trace()
    #print(par[0],pars[0])
    #plt.figure(dpi=500)
    #plt.plot(x1,y1,linestyle='',marker='.')
    #plt.plot(x1,pars[0]*x1+pars[1])
    
    return pars[0], 10**inner_dens

def get_density_and_slope(x,y,scale,phys_scale,scale_mult):
    
    # add a bit of logic to skip galaxies with no coverage inside the phys_scale
    if x[0] >= np.log10(phys_scale):
        return 0,0,0,0
    
    # decide between using the HST scale (fixed angle) or a fixed physical scale
    if scale_mult*scale < phys_scale: 
        radius = phys_scale
    else:
        radius = scale_mult*scale
    max_ind = find_nearest(x, np.log10(radius))  
    x1 = x[0:max_ind+1]
    y1 = y[0:max_ind+1]
    pars, cov = curve_fit(f=power_law, xdata=x1, ydata=y1)
    #par = np.polyfit(x1,y1,deg=1)
    
    # compute the MAD of residuals from the fit
    y_fit = pars[0]*x1+pars[1]
    resids = y1-y_fit
    
    #plt.figure(dpi=500)
    #plt.plot(x1,y1,label='Original')
    #plt.plot(x1,y_fit, label='Fit')
    #plt.legend()
    
    res_mad = mad(resids)
    
    #pdb.set_trace()
    #plt.figure(dpi=500)
    #plt.plot(x1,y1)
    #plt.plot(x1, power_law(x1,pars[0],pars[1]))
    ind_5pc = find_nearest(x, np.log10(5))
    inner_dens = y[ind_5pc]
    #pdb.set_trace()
    #print(par[0],pars[0])
    #plt.figure(dpi=500)
    #plt.plot(x1,y1,linestyle='',marker='.')
    #plt.plot(x1,pars[0]*x1+pars[1])
    
    return pars[0], 10**inner_dens, res_mad, resids


# =============================================================================

# =============================================================================

def get_our_data(slope_ext,phys_ext,cut):
    gal_file = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data_'+slope_ext+phys_ext+'.fits'
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
    
    names_c = names[np.where(lograds[:,0] <= np.log10(5))]
    slopes_c = slopes[np.where(lograds[:,0] <= np.log10(5))]
    cen_dens_c = cen_dens[np.where(lograds[:,0] <= np.log10(5))]
    vmags_c = vmags[np.where(lograds[:,0] <= np.log10(5))]
    stone_slopes_c = stone_slopes[np.where(lograds[:,0] <= np.log10(5))]
    stone_dens_at_5pc_c = stone_dens_at_5pc[np.where(lograds[:,0] <= np.log10(5))]
    NSC_comp_c = NSC_comp[np.where(lograds[:,0] <= np.log10(5))]
    lograds_c = lograds[np.where(lograds[:,0] <= np.log10(5))[0],:]
    logdens_c = logdens[np.where(lograds[:,0] <= np.log10(5))[0],:]
    dists_c = dists[np.where(lograds[:,0] <= np.log10(5))[0]]
    #pdb.set_trace()
    
    if cut:
        return names_c, slopes_c, cen_dens_c, vmags_c, stone_slopes_c, \
            stone_dens_at_5pc_c, NSC_comp_c, lograds_c, logdens_c, dists_c
    else:
        return names, slopes, cen_dens, vmags, stone_slopes, \
            stone_dens_at_5pc, NSC_comp, lograds, logdens, dists
            
# =============================================================================



# =============================================================================
# a little function for finding outliers
def get_n_max(arr, n):
    temp_arr = np.copy(arr)
    maxes = np.zeros((n))
    maxes_idx = np.zeros((n)).astype(int)
    #pdb.set_trace()
    for i in range(n):
        maxes[i] = np.max(temp_arr)
        maxes_idx[i] = np.argmax(temp_arr)
        #pdb.set_trace()
        temp_arr[maxes_idx[i]] = -999999
    return maxes, maxes_idx

# =============================================================================

# =============================================================================
# a little function for finding outliers
def get_n_min(arr, n):
    temp_arr = np.copy(arr)
    mins = np.zeros((n))
    mins_idx = np.zeros((n)).astype(int)
    #pdb.set_trace()
    for i in range(n):
        mins[i] = np.max(temp_arr)
        mins_idx[i] = np.argmax(temp_arr)
        #pdb.set_trace()
        temp_arr[mins_idx[i]] = -999999
    return mins, mins_idx

# =============================================================================



dist = get_galaxy_distance('DDO084')

