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
import astropy.units as uni
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import csv
from   numpy.linalg import eig, inv

#%%
# =============================================================================
# function to return M/L ratio using NSA colors with Taylor+2011 color-M/L relation
def get_nsa_M_L(gal_name, ra1, dec1, dist):

    all_MLs = []
    all_ML_types = []    

    # Read in table for V-I colors for galaxies from Lauer+05
    names = []
    v_i = []
    lauer_color_filename = '../Data_Sets/Lauer2005_data/color_data_Lauer_05.txt'
    # fields: name	instrument	channel	filter	pa	l_pa	u_pa	ell	l_ell	u_ell	n	l_n	u_n	reff	l_reff	
    #         u_reff	ieff	l_ieff	u_ieff	m	l_m	u_m	v	u_v	i	u_i	logm	u_logm
    with open(lauer_color_filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            #pdb.set_trace()
            if row[0] != 'Name': # marks the start of data row
                names.append(row[0])
                v_i.append(float(row[1]))
       
    names = np.array(names)
    v_i = np.array(v_i)
    ind_v_i = np.where(names == gal_name)[0]
        
    # read in data from David's Table
    filename = '../Data_Sets/Nasa_Sloan_Atlas/nsa_v1_0_1.fits'
    hdul = fits.open(filename)  # open a FITS file
    data = hdul[1].data
    hdul.close()
    
    ra2 = np.array(data['RA'])
    dec2 = np.array(data['Dec'])

    c = SkyCoord(ra=ra1*uni.degree, dec=dec1*uni.degree)
    catalog = SkyCoord(ra=ra2*uni.degree, dec=dec2*uni.degree)
    idx, d2d, d3d = match_coordinates_sky(c,catalog)

    max_sep = (1./60.) * uni.degree
    #pdb.set_trace()
    if d2d > max_sep:
        print('*** ERROR ***: Max Separation Reached.')
        print(gal_name+" not found in David's table.")
        return -999, -999, np.array([-999]), np.array([-999]), -999

    if data['DFLAGS'][idx,3] > 32:
        print('** ERROR **: "bad" flag for g-band')
        return -999, -999, np.array([-999]), np.array([-999]), -999
    if data['DFLAGS'][idx,5] > 32:
        print('** ERROR **: "bad" flag for i-band')
        return -999, -999, np.array([-999]), np.array([-999]), -999

    name = data['IAUNAME'][idx]
    petros = data['FIBER_FLUX'][idx]

    if petros.shape == (1,7):
        petros = petros[0]

    # store petrosian fluxes from SDSS in nanomaggies
    pflux_g = petros[3]
    pflux_i = petros[5]
    #use formula from SDSS site to calculate pogson magnitudes from petrosian flux
    pmag_g = 22.5 - (2.5*np.log10(pflux_g))
    pmag_i = 22.5 - (2.5*np.log10(pflux_i))
    
    # calculate array of distance moduli from my calculated best distance
    DM = 5*np.log10(dist) - 5 # dist in pc
    
    # store given extinction values
    # multiplying extinction by 0.86 to convert from Schlegel, Finkbeiner & 
    #  Davis extincion to Schlafly and Finkbeiner
    ext_g = data['EXTINCTION'][idx, 3] * 0.86
    ext_i = data['EXTINCTION'][idx, 5] * 0.86
    
    # calculate absolute magnitudes for each band
    gMag = pmag_g - DM - ext_g
    iMag = pmag_i - DM - ext_i

    # convert colors from vega system to AB (m_AB - m_vega = 0.34) (not necessary)
    #gMag_AB = 0.34 + gMag
    #iMag_AB = 0.34 + iMag

    # true g-i color
    gicolor = gMag - iMag
    ML_type = 5
    
    # use relation from Taylor+2011 for g-i to M/L_i
    logM_L = -0.68 + 0.70*gicolor
    M_L = 10**logM_L
    
    # store both versions of M/L for comparison
    all_MLs.append(M_L)
    all_ML_types.append(ML_type)
    
    if len(ind_v_i) != 0:
        ML_type = 3
        if v_i[ind_v_i] <= 1.8:
            gicolor = 1.481*v_i[ind_v_i] - 0.536
        else:
            gicolor = 0.83*v_i[ind_v_i] + 0.600
        # use relation from Taylor+2011 for g-i to M/L_i
        logM_L = -0.68 + 0.70*gicolor
        M_L = 10**logM_L
        # store both versions of M/L for comparison
        all_MLs.append(M_L)
        all_ML_types.append(ML_type)

    # convert from AB sytem to solar units (incorrect)
    #M_L = 10**logM_L * 4.58
    #pdb.set_trace()
    return M_L, ML_type, np.array(all_MLs), np.array(all_ML_types), gicolor 
 
# =============================================================================


# =============================================================================
# Function to parse David's Table for galaxy data based on coordinates
def get_David_gal_data(gal_name, ra1, dec1):
    
    # read in data from David's Table
    filename = '../Data_Sets/David_Tables/catalog.fits'
    # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
    hdul = fits.open(filename)  # open a FITS file
    data = hdul[1].data
    hdul.close()
    
    ra2 = np.array(data['RA'])
    dec2 = np.array(data['Dec'])
    
    c = SkyCoord(ra=ra1*uni.degree, dec=dec1*uni.degree)
    catalog = SkyCoord(ra=ra2*uni.degree, dec=dec2*uni.degree)
    idx, d2d, d3d = match_coordinates_sky(c,catalog)

    max_sep = (1./60.) * uni.degree
    if d2d > max_sep:
        print('*** ERROR ***: Max Separation Reached.')
        print(gal_name+" not found in David's table.")
        return '', -999,-999,-999, -999,-999
    
    #print('Match Success: Separation = {:.9f}'.format(d2d[0]))
    #pdb.set_trace()
    name = data['objname'][idx]
    dist = data['bestdist'][idx]
    logmass = data['logmass'][idx]
    vmag = -999#data['i_lum_nsa'][idx]
    gi_color = data['gi_color'][idx]
    if data['best_type'][idx] == 'early':
        galtype = 0
    elif data['best_type'][idx] == 'late':
        galtype = 1
    
    
    # read in data from David's Table with vmags
    filename1 = '../Data_Sets/David_Tables/sample_4.fits'
    # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
    hdul = fits.open(filename1)  # open a FITS file
    data1 = hdul[1].data
    hdul.close()
    ra3 = np.array(data1['RA'])
    dec3 = np.array(data1['Dec'])
    catalog1 = SkyCoord(ra=ra3*uni.degree, dec=dec3*uni.degree)
    idx, d2d, d3d = match_coordinates_sky(c,catalog1)
    vmag = data1['VMag'][idx]
    
    
    return name, dist, logmass, vmag, gi_color, galtype

# =============================================================================


# =============================================================================
# function to get the log stellar mass for a galaxy from david's table
def get_galaxy_logmass(gal_name): # use no space gal name 
    # read in fits table with galaxy masses from David
    filename = '../Data_Sets/David_Tables/catalog.fits'
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
    if gal_name == 'bts76':
        gal_name = 'SDSSJ115844.10+273506.0'
    if gal_name == 'ddo084':
        gal_name = 'UGC05829'
    if gal_name == 'kk2000-03':
        gal_name = 'PGC009140'
    if gal_name == 'kk2000-53':
        gal_name = '[KK2000] 53'
    if gal_name == 'kk96':
        gal_name = 'KK96'
    if gal_name == 'leg09':
        gal_name = 'LeG09'
    if gal_name == 'lvj1217+47':
        gal_name = 'LV J1217+4703'
    if gal_name == 'ngc5011c':
        gal_name = 'ESO269-068'
    if gal_name == 'pgc4310323':
        gal_name = 'SDSSJ120531.04+310434.1'
    if gal_name == 'ugc07242':
        gal_name = 'UGC07242'
    
        
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
    filename = '../Data_Sets/David_Tables/catalog.fits'
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
    if gal_name == 'bts76':
        gal_name = 'SDSSJ115844.10+273506.0'
    if gal_name == 'ddo084':
        gal_name = 'UGC05829'
    if gal_name == 'kk2000-03':
        gal_name = 'PGC009140'
    if gal_name == 'kk2000-53':
        gal_name = '[KK2000] 53'
    if gal_name == 'kk96':
        gal_name = 'KK96'
    if gal_name == 'leg09':
        gal_name = 'LeG09'
    if gal_name == 'lvj1217+47':
        gal_name = 'LV J1217+4703'
    if gal_name == 'ngc5011c':
        gal_name = 'ESO269-068'
    if gal_name == 'pgc4310323':
        gal_name = 'SDSSJ120531.04+310434.1'
    if gal_name == 'ugc07242':
        gal_name = 'UGC07242'
    
    
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
    filename = '../Data_Sets/David_Tables/sample_4.fits'
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
    if gal_name == 'bts76':
        gal_name = 'SDSSJ115844.10+273506.0'
    if gal_name == 'ddo084':
        gal_name = 'UGC05829'
    if gal_name == 'kk2000-03':
        gal_name = 'PGC009140'
    if gal_name == 'kk2000-53':
        gal_name = '[KK2000] 53'
    if gal_name == 'kk96':
        gal_name = 'KK96'
    if gal_name == 'leg09':
        gal_name = 'LeG09'
    if gal_name == 'lvj1217+47':
        gal_name = 'LV J1217+4703'
    if gal_name == 'ngc5011c':
        gal_name = 'ESO269-068'
    if gal_name == 'pgc4310323':
        gal_name = 'SDSSJ120531.04+310434.1'
    if gal_name == 'ugc07242':
        gal_name = 'UGC07242'
    
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
    filename = '../Data_Sets/David_Tables/catalog.fits'
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
        gal_name = 'PGC041959'
    if gal_name == 'VCC1199':
        gal_name = 'PGC041180'
    if gal_name == 'bts76':
        gal_name = 'SDSSJ115844.10+273506.0'
    if gal_name == 'ddo084':
        gal_name = 'UGC05829'
    if gal_name == 'kk2000-03':
        gal_name = 'PGC009140'
    if gal_name == 'kk2000-53':
        gal_name = '[KK2000] 53'
    if gal_name == 'kk96':
        gal_name = 'KK96'
    if gal_name == 'leg09':
        gal_name = 'LeG09'
    if gal_name == 'lvj1217+47':
        gal_name = 'LV J1217+4703'
    if gal_name == 'ngc5011c':
        gal_name = 'ESO269-068'
    if gal_name == 'pgc4310323':
        gal_name = 'SDSSJ120531.04+310434.1'
    if gal_name == 'ugc07242':
        gal_name = 'UGC07242'
    
    
    ind = np.where(data['objname'] == gal_name)
    if ind[0].size == 0:
        print(gal_name+" not found in David's table")
        return 0, 0
    else:
        return data['RA'][ind][0], data['Dec'][ind][0]
# =============================================================================


# =============================================================================
# function to get the M/L for a galaxy from Renuka_ML.dat file
def get_renuka_ML(gal_name): # use no space gal name 
    
    file_path = '../Data_Sets/Pechetti17_data/Renuka_ML.dat'
    f = open(file_path, 'r')
    all_lines = f.readlines()
    f.close()
    
    names = []
    MLs = []
    for i in range(len(all_lines)-1):
        split = all_lines[i+1].split()
        names.append(split[0])
        MLs.append(float(split[3]))
        
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

    file_path = '../Data_Sets/Stone_Metzger_data/RenukaAnil.dat'
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
    
    file_path = '../Data_Sets/Stone_Metzger_data/RenukaAnil.dat'
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
    gal_file = '../Result_Tables/all_gal_data_'+slope_ext+phys_ext+'.fits'
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
    all_MLs = dat['all_mls']
    all_ML_types = dat['all_ml_types']
    dist_flags = dat['dist_flag']
    
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
    all_MLs_c = all_MLs[np.where(lograds[:,0] <= np.log10(5))[0],:]
    all_ML_types_c = all_ML_types[np.where(lograds[:,0] <= np.log10(5))[0],:]
    MLs_c = MLs[np.where(lograds[:,0] <= np.log10(5))]
    ML_types_c = ML_types[np.where(lograds[:,0] <= np.log10(5))]
    dist_flags_c = dist_flags[np.where(lograds[:,0] <= np.log10(5))]
    #pdb.set_trace()
    
    if cut:
        return names_c, slopes_c, cen_dens_c, vmags_c, stone_slopes_c, \
            stone_dens_at_5pc_c, NSC_comp_c, lograds_c, logdens_c, dists_c, \
            all_MLs_c, all_ML_types_c, MLs_c, ML_types_c, dist_flags_c
    else:
        return names, slopes, cen_dens, vmags, stone_slopes, \
            stone_dens_at_5pc, NSC_comp, lograds, logdens, dists, \
            all_MLs, all_ML_types, MLs, ML_types, dist_flags
            
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

# =============================================================================
# function to retrieve the F814W zero-point magnitude and exposure time 
# for a galaxy from Hoyer+22
def get_zp_and_texp(name):
    mydict = {"bts76": (25.509,1030),
              "ddo084": (25.507,760),
              "eso553-046": (25.512,900),
              "kk2000-03": (25.510,1200),
              "kk2000-53": (25.510,1000),
              "kk96": (25.509,1096),
              "leg09": (25.509,1096),
              "lvj1217+4703": (25.509,1164),
              "ngc5011c": (25.512,900),
              "pgc4310323": (25.507,760),
              "ugc07242": (25.525,900)}
    
    return mydict[name]
 
# =============================================================================


# =============================================================================
# functions for fitting ellipsoid from the web
def ls_ellipsoid(xx,yy,zz):                                  
    #finds best fit ellipsoid. Found at http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
    #least squares fit to a 3D-ellipsoid
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz  = 1
    #
    # Note that sometimes it is expressed as a solution to
    #  Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz  = 1
    # where the last six terms have a factor of 2 in them
    # This is in anticipation of forming a matrix with the polynomial coefficients.
    # Those terms with factors of 2 are all off diagonal elements.  These contribute
    # two terms when multiplied out (symmetric) so would need to be divided by two
    
    # change xx from vector of length N to Nx1 matrix so we can use hstack
    x = xx[:,np.newaxis]
    y = yy[:,np.newaxis]
    z = zz[:,np.newaxis]
    
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
    J = np.hstack((x*x,y*y,z*z,x*y,x*z,y*z, x, y, z))
    K = np.ones_like(x) #column of ones
    
    #np.hstack performs a loop over all samples and creates
    #a row in J for each x,y,z sample:
    # J[ix,0] = x[ix]*x[ix]
    # J[ix,1] = y[ix]*y[ix]
    # etc.
    
    JT=J.transpose()
    JTJ = np.dot(JT,J)
    InvJTJ=np.linalg.inv(JTJ);
    ABC= np.dot(InvJTJ, np.dot(JT,K))

    # Rearrange, move the 1 to the other side
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
    #    or
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
    #  where J = -1
    eansa=np.append(ABC,-1)

    return (eansa)

def polyToParams3D(vec,printMe):                             
    #gets 3D parameters of an ellipsoid. Found at http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
    # convert the polynomial form of the 3D-ellipsoid to parameters
    # center, axes, and transformation matrix
    # vec is the vector whose elements are the polynomial
    # coefficients A..J
    # returns (center, axes, rotation matrix)
    
    #Algebraic form: X.T * Amat * X --> polynomial form
    
    if printMe: print('\npolynomial\n',vec)
    
    Amat=np.array(
    [
    [ vec[0],     vec[3]/2.0, vec[4]/2.0, vec[6]/2.0 ],
    [ vec[3]/2.0, vec[1],     vec[5]/2.0, vec[7]/2.0 ],
    [ vec[4]/2.0, vec[5]/2.0, vec[2],     vec[8]/2.0 ],
    [ vec[6]/2.0, vec[7]/2.0, vec[8]/2.0, vec[9]     ]
    ])
    
    if printMe: print('\nAlgebraic form of polynomial\n',Amat)
    
    #See B.Bartoni, Preprint SMU-HEP-10-14 Multi-dimensional Ellipsoidal Fitting
    # equation 20 for the following method for finding the center
    A3=Amat[0:3,0:3]
    A3inv=inv(A3)
    ofs=vec[6:9]/2.0
    center=-np.dot(A3inv,ofs)
    if printMe: print('\nCenter at:',center)
    
    # Center the ellipsoid at the origin
    Tofs=np.eye(4)
    Tofs[3,0:3]=center
    R = np.dot(Tofs,np.dot(Amat,Tofs.T))
    if printMe: print('\nAlgebraic form translated to center\n',R,'\n')
    
    R3=R[0:3,0:3]
    R3test=R3/R3[0,0]
    # print('normed \n',R3test)
    s1=-R[3, 3]
    R3S=R3/s1
    (el,ec)=eig(R3S)
    
    recip=1.0/np.abs(el)
    axes=np.sqrt(recip)
    if printMe: print('\nAxes are\n',axes  ,'\n')
    
    inve=inv(ec) #inverse is actually the transpose here
    if printMe: print('\nRotation matrix\n',inve)
    return (center,axes,inve)
# =============================================================================



# ['bts76', 'ddo084', 'kk2000-03', 'kk2000-53', 'kk96', 'leg09',
# 'lvj1217+4703', 'ngc5011c', 'pgc4310323', 'ugc07242']
# =============================================================================
# a = get_David_gal_data('bts76', 179.68375, 27.585)
# b = get_David_gal_data('ddo084', 160.67458, 34.44889)
# c = get_David_gal_data('kk2000-03', 36.17792, -73.51278)
# d = get_David_gal_data('kk2000-53', 197.80917, -38.90611)
# e = get_David_gal_data('kk96', 162.61292, 12.36083)
# f = get_David_gal_data('leg09', 160.64417, 12.15056)
# g = get_David_gal_data('lvj1217+4703', 184.29208, 47.06361)
# h = get_David_gal_data('ngc5011c', 198.29958, -43.26556)
# i = get_David_gal_data('pgc4310323', 181.37917, 31.07611)
# j = get_David_gal_data('ugc07242', 183.53083, 66.09222)
# k = get_David_gal_data('eso553-046', 81.77375, -20.67806)
# =============================================================================

#a = get_galaxy_distance('A1020-M1')
#b = get_galaxy_RA_and_Dec('NGC0821')
#a = get_nsa_M_L('NGC 821', 32.08767, 10.99475, 22400000.0)
#a = get_nsa_M_L('NGC 821', 32.08803, 10.9949848, 22400000.0)
#a = get_David_gal_data('kk2000-03', 358.829292, 5.915833)
#b = get_galaxy_RA_and_Dec('NGC1700')



