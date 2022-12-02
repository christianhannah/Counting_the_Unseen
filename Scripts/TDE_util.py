#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 15:05:13 2022

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
from scipy import integrate

# =============================================================================
# function to retreive the data for the model galaxy sample
def get_model_data(filename):
    hdul = fits.open(filename)
    #head = hdul[0].data
    dat = hdul[1].data
    hdul.close()

    gal_numbers = dat['no.']
    gal_types = dat['type']
    gal_dists = dat['dist']
    gal_logm = dat['galmass']
    gal_logm_BH = dat['bhmass']
    gal_gamma = dat['gamma']
    gal_rho = dat['rho5pc']
    
    
    return gal_numbers, gal_types, gal_dists, gal_logm, gal_logm_BH, \
        gal_gamma, gal_rho

# =============================================================================



# =============================================================================
# function to return power-law inner density profile
    def rho_r(r,slope,dens5pc):
        return dens5pc*(r/5.0)**(slope)
    
# =============================================================================

# =============================================================================
# function to return power-law inner density profile
    def get_log_rho(r,slope,dens5pc,delta):
        return np.log10(dens5pc*(r/5.0)**(slope)+delta)
    
# =============================================================================


# =============================================================================
# derive the total gravitational potential (psi(r)) as a function of r
    def psi_r(r,M_BH,slope,dens5pc):
        G = 4.301e-3 # pc (km/s)^2 M_sol^-1
        psi_1 = G*M_BH/r
        
        M_enc = np.zeros_like(r)
        for i in range(len(M_enc)):
            M_enc[i] = get_enc_mass(r,slope,dens5pc,i)
        psi_2 = G*M_enc/r
        
        psi_3 = np.zeros_like(r)
        for i in range(len(psi_3)):
            psi_3[i] = get_ext_potential(r,slope,dens5pc,i)
            
        return psi_1+psi_2+psi_3

# =============================================================================



# =============================================================================
# function to compute the mass enclosed from density profile 
def get_enc_mass(r,slope,dens5pc,max_ind):
    y = dens5pc*(r[0:max_ind+1]/5)**(slope+2)
    return 4*np.pi*integrate.trapezoid(y, r[0:max_ind+1])
    
# =============================================================================
    


# =============================================================================
# function to compute the contribution to the potential of the galaxy at 
# larger radii
def get_ext_potential(r,slope,dens5pc,min_ind):
    G = 4.301e-3 # pc (km/s)^2 M_sol^-1
    y = dens5pc*(r[min_ind:]/5)**(slope+1)  
    return 4*np.pi*G*integrate.trapezoid(y,r[min_ind:])

# =============================================================================



# =============================================================================
# define function to find value in an array nearest to a supplied value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    