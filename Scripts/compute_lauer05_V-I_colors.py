#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 23:59:20 2022

Code to use the SB profiles from Lauer+05 in F555W and F814W to compute the 
V-I colors of these galaxies within a radius of 0.5"

@author: christian
"""

import numpy as np 
from astropy.io import fits
import mge1d_util as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as uni
from dustmaps.sfd import SFDQuery
from dustmaps.config import config
config['data_dir'] = '../Data_Sets/Dustmaps/'

# conversion factor between radians and arcsec
arcsec_to_radian = 4.8481e-6


# define the names of the galaxies for which colors must be computed
names = np.array(['NGC 1023','NGC 1427','NGC 1439','NGC 2434','NGC 3115',
                  'NGC 3384','NGC 3945','NGC 4026','NGC 4278','NGC 4365',
                  'NGC 4552','NGC 7213','IC 1459','NGC 584','NGC 821',
                  'NGC 2300','NGC 3377','NGC 3379','NGC 3607','NGC 3608',
                  'NGC 4291','NGC 4382','NGC 4472','NGC 4473','NGC 4486B',
                  'NGC 4494','NGC 4589','NGC 4621','NGC 4649','NGC 5576',
                  'NGC 5813','NGC 5982','NGC 7727'])


# get the F555W SB data
hdul = fits.open('../Data_Sets/Lauer2005_data/SBs_Lauer_2005.fit')
data = hdul[1].data
hdul.close()

# read in table with RA and DEC information to get David's distances
hdul = fits.open('../Data_Sets/Lauer2005_data/gal_properties_Lauer_2005.fit')
data1 = hdul[1].data 
hdul.close()


# define the number of radii desired in interpolated profiles
num_rad = 100

# store and plot V-band SB data
SB_V = np.zeros((len(names),num_rad))
rad_V = np.zeros((len(names),num_rad))
distances = np.zeros(len(names))
plt.figure(dpi=500)
for i in range(len(names)):
    coord_ind = np.where(data1['name'] == names[i])[0]
    inds = np.where(data['name'] == names[i])[0]
    filt = np.unique(data['filt'][inds])
    # let's get the distance from David's table to use instead
    distances[i] = u.get_David_gal_data(names[i], data1['_RA'][coord_ind], data1['_DE'][coord_ind])[1]*10**6 # in pc
    SB_init = data['SuBr'][inds]
    rad_init = (data['MajAxis'][inds]*arcsec_to_radian)*distances[i] # convert radii from arcsec to pc
    # interpolate the SB data
    rad_V[i,:] = np.geomspace(rad_init[0], rad_init[-1], num_rad)
    SB_V[i,:] = 10**np.interp(np.log10(rad_V[i,:]), np.log10(rad_init), np.log10(SB_init))
    
    # correct the surface brightness profiles for extinction
    coords = SkyCoord(ra=data1['_RA'][coord_ind]*uni.degree, dec=data1['_DE'][coord_ind]*uni.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass 
    # (values from Table 6 in Schlafly & Finkbeiner (2011))
    # all WFPC2
    if filt == 'F555W':
        A_v = E_BV*2.755
    elif filt == 'F606W':
        A_v = E_BV* 2.415
    elif filt == 'F814W':
        A_v = E_BV*1.549
    SB_V[i,:] = SB_V[i,:] - A_v
    
    plt.plot(np.log10(rad_V[i,:]),SB_V[i,:])
plt.ylim(25,9)
plt.plot([np.log10(5),np.log10(5)],[25,9],'--k', label='r = 5pc')
plt.xlabel('log(Radius [pc])')   
plt.ylabel('Surface Brightness [mag/arcsec$^2$]')
plt.legend()

    
    
    
    
    
    
    
    
    
    
    