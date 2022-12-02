#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 12:53:28 2022

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
from scipy import integrate
warnings.filterwarnings("ignore")


# ============================ FUNCTIONS ======================================

# =============================================================================
# function to compute the mass enclosed from density profile 
def get_enc_mass(r,slope,rho_infl,r_infl,max_ind):
    y = rho_infl*(r[0:max_ind+1]/r_infl)**(slope+2)
    return 4*np.pi*integrate.trapezoid(y, r[0:max_ind+1])
# =============================================================================

# =============================================================================
# function to return power-law inner density profile
def get_rho_r(r,slope,rho_infl,r_infl):
    return rho_infl*(r/r_infl)**(-slope)
# =============================================================================

# =============================================================================
# function to return power-law inner density profile
def get_log_rho(r,slope,rho_infl,r_infl,delta):
    return np.log10(rho_infl*(r/r_infl)**(slope)+delta) 
# =============================================================================

# =============================================================================
# function to compute the contribution to the potential of the galaxy at 
# larger radii
def get_ext_potential(r,slope,rho_infl,r_infl,min_ind):
    G = 4.301e-3 # pc (km/s)^2 M_sol^-1
    y = rho_infl*(r[min_ind:]/r_infl)**(slope+1)  
    return 4*np.pi*G*integrate.trapezoid(y,r[min_ind:])
# =============================================================================

# =============================================================================
# derive the total gravitational potential (psi(r)) as a function of r
def get_psi_r(r,M_BH,slope,rho_infl,r_infl):
    G = 4.301e-3 # pc (km/s)^2 M_sol^-1
    psi_1 = G*M_BH/r
     
    M_enc = np.zeros_like(r)
    for i in range(len(M_enc)):
        M_enc[i] = get_enc_mass(r,slope,rho_infl,r_infl,i)
    psi_2 = G*M_enc/r
    
    psi_3 = np.zeros_like(r)
    for i in range(len(psi_3)):
        psi_3[i] = get_ext_potential(r,slope,rho_infl,r_infl,i)
            
    return psi_1+psi_2+psi_3
# =============================================================================

# =============================================================================
# define function to find value in an array nearest to a supplied value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# =============================================================================   

# =========================== END OF FUNCTIONS ================================


#let's define our radial sampling
min_rad = 1e-6
max_rad = 10**6
n_rad = 10**6
r = np.linspace(min_rad,max_rad,n_rad) # in pc
# logarithmically spaced grid as well
num_rad = 10**6
r_ls = np.geomspace(min_rad,max_rad,num_rad)


# CONSTANTS
#G = 4.301e-3 # pc (km/s)^2 M_sol^-1
G = 4.301e-3 * 10**6 # pc (m/s)^2 M_sol^-1
#%%

    
# specify density power-law slope, influence radius, and density at the 
# influence radius
slope = 1.6
rho_infl = 3e3 # M_sol/pc^3
r_infl = 15 # pc
r_infl_ind = u.find_nearest(r,r_infl)

# use the enclosed mass of the density profile out to r_infl to define the BH mass
M_BH = get_enc_mass(r,slope,rho_infl,r_infl,r_infl_ind) # M_sol
    
#%%

# get the total gravitational potential (psi) as a function of r
# In this toy model, psi(r) is just GM_BH/r
psi_r_init = G*M_BH/r # m^2/s^2
rho_r_init = get_rho_r(r,slope,rho_infl,r_infl)



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
psi_r_ls = np.interp(r_ls, r, psi_r_init)
    
# compute log(rho(r)+delta) for logarithmically spaced radii
delta = 10**(-6)
log_rho = get_log_rho(r_ls,slope,rho_infl,r_infl,delta)

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
# STEP 1: tabulate rho as a function of psi
#rho_r = rho_infl*((G*M_BH)/(psi_r*r_infl))**(-slope)
 
# function to return d^2rho/dpsi^2
def get_dsqrho_dpsisqr(psi,slope,rho_infl,r_infl,M_BH,G):
    return ((slope-1)*slope*rho_infl*((G*M_BH)/(r_infl*psi))**(-slope))/psi**2
  
# =============================================================================
# plt.figure(dpi=500)
# plt.plot(psi_r,rho_r)
# plt.xlabel('$\psi$',fontsize=16)
# plt.ylabel('$\\rho$',fontsize=16)
# =============================================================================
     
e = np.linspace(100,np.max(psi_r_init)*10,10**6)
#e = np.linspace(0.01,10**20,1000)

f = np.zeros(len(e))

for i in range(len(e)):
    
    # first define new array of psi values from 0 to epsilon for integration in DF
    psi_r = np.linspace(0,e[i],10**6)

    # interpolate to get rho values for the new psi grid
    rho_r = np.interp(psi_r,psi_r_init,rho_r_init)
    
    
    psi_r_ls = np.interp(r_ls, r, psi_r_init)
    
    # compute log(rho(r)+delta) for logarithmically spaced radii
    delta = 10**(-6)
    log_rho = get_log_rho(r_ls,slope,rho_infl,r_infl,delta)

    # use standard finite-difference derivatives to compute 
    # d(log(rho(r)+delta))/d(psi(r))
    dlrho_dpsi = np.zeros(len(log_rho)-1)
    #r_der = np.zeros_like(dlrho_dpsi)
    psi_der = np.zeros_like(dlrho_dpsi)
    rho_p_delta = np.zeros_like(dlrho_dpsi)
    for i in range(len(dlrho_dpsi)):
        diff_1 = log_rho[i+1] - log_rho[i]
        diff_2 = psi_r_ls[i+1] - psi_r_ls[i]
        #r_der = r_ls[i+1]-r_ls[i]
        psi_der[i] = psi_r_ls[i]+ (psi_r_ls[i+1] - psi_r_ls[i])/2
        rho_p_delta = log_rho[i] + diff_1/2
        dlrho_dpsi[i] = diff_1/diff_2
    
    # get d(rho(r))/d(psi(r)) via 
    # d(rho(r))/d(psi(r)) = [rho(r)+delta]*d(log(rho(r)+delta))/d(psi(r))
    drho_dpsi_r = rho_p_delta*dlrho_dpsi
    
    plt.figure(dpi=500)
    plt.plot(psi_der,drho_dpsi_r)
    
    pdb.set_trace()
    # perform change of variables so our integral goes from -1 to 1 instead of 
    # 0 to epsilon so I can use the double exponential transfom to remove the 
    # boundary singularity
    psi_prime_r = 2*(psi_r/e[i])-1
    Q = 1/(e[i]-psi_r)
    
    
    plt.figure(dpi=500)
    plt.plot(psi_prime_r,rho_r,linestyle='',marker='.',color='k')
    plt.xlabel('$\psi$'+"'")
    plt.ylabel('$\\rho$')
    
    
    # change of variable for the double exponential transform
    t_top = 2
    t_bot = -t_top
    t = np.linspace(t_bot,t_top,10**6)
    rho_t = rho_infl*(-(2*G*M_BH)/(r_infl*e[i]*(-np.tanh(np.pi/2*np.sinh(t)))-1))**(-slope)
    
    plt.figure(dpi=500)
    plt.plot(t,rho_t)
    plt.xlabel('t')
    plt.ylabel('$\\rho$')
    
    
    #rho_new = np.interp()
    
    #t = np.arcsinh((2*np.arcsinh(psi_prime_r))/np.pi)
    
    
    pdb.set_trace()
    second_deriv = np.zeros(len(psi_prime_r)-2)    
    psi_p_sec_der = np.zeros_like(second_deriv)
    for k in range(len(second_deriv)):
        k+=1
        second_deriv[k] = (rho_r[k+1]-2*rho_r[k]+rho_r[k-1])/(psi_prime_r[k+1]-
                                                            psi_prime_r[k])**2
        psi_p_sec_der[k] = psi_prime_r[k]
                           
        
    
    
    
    part_1 = (4/e[i]**2)
    
    
    part_1 = get_dsqrho_dpsisqr(psi_r,slope,rho_infl,r_infl,M_BH,G)
    
    
    
    
    part_2 = np.sqrt(e[i]-psi_r)
    y = part_1/part_2
    f[i] = (1/(np.sqrt(8)*np.pi**2))*integrate.trapezoid(y,psi_r)
    pdb.set_trace()

plt.figure(dpi=500)
plt.plot(e,f)
plt.xlabel('$\epsilon$')
plt.ylabel('f($\epsilon$)')

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








