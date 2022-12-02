#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 17:54:52 2022

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
from tqdm import tqdm
import scipy
import warnings
import TDE_util as tu
warnings.filterwarnings("ignore")



# read in the model galaxy sample builder fits file
gsmf_ext = 'driver22'
BH_mass_ext = '_reines15'
vol_ext = '_R_30.00Mpc'
filename = '/Users/christian/OneDrive/Desktop/TDE Code/model_galaxy_sample_'+\
            gsmf_ext+BH_mass_ext+vol_ext+'.fits'
gal_numbers, gal_types, gal_dists, gal_logm, \
    gal_logm_BH, gal_gamma, gal_logrho = tu.get_model_data(filename)

# for a single run
gal_numbers = gal_numbers[0]


#let's define our radial sampling
min_rad = 0
max_rad = 10**3
rad_step = 5
r = np.arange(min_rad,max_rad,rad_step) # in pc
# logarithmically spaced grid as well
num_rad = 10**3
r_ls = np.geomspace(min_rad,max_rad,num_rad)


# CONSTANTS
G = 4.301e-3 # pc (km/s)^2 M_sol^-1

#%%

for i in tqdm(range(len(gal_numbers)), position=0, leave=True):
    
    # specify density power-law slope and density at 5pc to construct rho(r)
    slope = gal_gamma[i]
    dens5pc = 10**gal_logrho[i] # M_sol/pc^3
    M_BH = 10**gal_logm_BH[i] # M_sol
    M_gal = 10**gal_logm # M_sol
    

    # get the total gravitational potential (psi) as a function of r
    psi_r = tu.psi_r(r,M_BH,slope,dens5pc)


    # =============================================================================
    # Compute the stellar distribution function (DF) as a function of specific 
    # orbital energy (f(epsilon))
    # =============================================================================
    # =============================================================================
    # =============================================================================
    # STEP 1: Define radial grid that is logarithmically spaced and compute
    # psi(r) vs. log(rho(r)+delta) where delta is some small density offset 
    # to avoid taking the logarithm of zero. 
    
    # interpolate psi_r onto logarithmically spaced radii
    psi_r_ls = np.interp(r_ls, r, psi_r)
    
    # compute log(rho(r)+delta) for logarithmically spaced radii
    delta = 10**(-6)
    log_rho = tu.get_log_rho(r_ls,slope,dens5pc,delta)

    # use standard finite-difference derivatives to compute 
    # d(log(rho(r)+delta))/d(psi(r))
    dlrho_dpsi = np.zeros(len(log_rho)-1)
    r_der = np.zeros_like(dlrho_dpsi)
    rho_p_delta = np.zeros_like(dlrho_dpsi)
    for i in range(len(dlrho_dpsi)):
        diff_1 = log_rho[i+1] - log_rho[i]
        diff_2 = psi_r_ls[i+1] - psi_r_ls[i]
        r_der = r_ls[i+1]-r_ls[i]
        rho_p_delta = log_rho[i] + diff_1/2
        dlrho_dpsi[i] = diff_1/diff_2
    
    # get d(rho(r))/d(psi(r)) via 
    # d(rho(r))/d(psi(r)) = [rho(r)+delta]*d(log(rho(r)+delta))/d(psi(r))
    drho_dpsi_r = rho_p_delta*dlrho_dpsi
    # =============================================================================
    # =============================================================================
    # STEP 2: Change variables in DF equation from sqrt(psi-epsilon) to Q to 
    # remove square-root divergence
    
    
    # =============================================================================


    ### TO-DO ###
    # use the DF to calculate the orbit-averaged angular momentum diffusion 
    # coefficient for highly eccentric orbits
    # Note: I believe that here we want the local diffusion coefficient expressed
    # in terms of the DF moments



    ### TO-DO ###
    # Compute the flux of stars that scatter into the loss cone per unit time 
    # and energy




    ### TO-DO ###
    # Compute the TDE rate by integrating the total flux into the loss cone for 
    # for stars of a given mass, and then by integrating over the stellar mass function.








