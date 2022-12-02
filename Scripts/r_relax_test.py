#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:46:10 2022

@author: christian
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import pdb
from scipy import integrate
import mge1d_util as u

# =============================================================================
# function to compute the mass enclosed from density profile 
# =============================================================================
# def get_enc_mass(r,slope,dens5pc,max_ind):
#     y = dens5pc*(r[0:max_ind+1]/5)**(slope+2)
#     return 4*np.pi*integrate.trapezoid(y, r[0:max_ind+1])
# =============================================================================

def get_enc_mass(r,y,max_ind):
    #y = dens5pc*(r[0:max_ind+1]/5)**(slope+2)
    return 4*np.pi*integrate.trapezoid(r[0:max_ind+1]**2*y[0:max_ind+1], r[0:max_ind+1])
    
# =============================================================================



# CONSTANTS
pc_to_m = 3.08567758128e16 
M_sol_to_kg = 1.989e30
years_to_sec = 3.154e7
M_sol =  1.989e30 # kg
#M_sol = 1
R_sol = 696.34e6 #m
#R_sol = 696.34e6/pc_to_m

G = 6.67e-11 
#G = 4.301e-3 * 10**6 # pc (m/s)^2 M_sol^-1
M_star = 1*M_sol

# =============================================================================
# # read in fits table with model galaxy parameters
# model_gal_filename = '../Result_Tables/model_galaxy_sample_driver22_reines15_R_30.00Mpc.fits'
# hdul = fits.open(model_gal_filename)
# model_data = hdul[1].data 
# hdul.close()
# 
# types = model_data['type']
# gammas = model_data['gamma']
# rhos = 10**model_data['rho5pc'] * M_sol_to_kg/pc_to_m**3
# log_rhos = model_data['rho5pc']
# bh_masses = 10**model_data['bhmass'] * M_sol_to_kg
# log_bh_masses = model_data['bhmass']
# gal_masses = model_data['galmass']
# 
# et_age = 10
# lt_age = 5
# t_ages = np.zeros(len(types))
# for i in range(len(t_ages)):
#     if types[i] == 0:
#         t_ages[i] = et_age*10**9*years_to_sec
#     else:
#         t_ages[i] = lt_age*10**9*years_to_sec
# 
# kappa = 1
# r5pc = 5*pc_to_m
# r_relaxes = np.zeros(len(types))
# for i in range(len(r_relaxes)):
#     r_relaxes[i] = ((G**2*rhos[i]*r5pc**(-gammas[i])*M_star*np.log(0.4*bh_masses[i]/M_star)*t_ages[i])/(kappa*(G*bh_masses[i])**(3/2)))**(1/(-gammas[i]-3/2))
#     r_relaxes[i] /= pc_to_m
# 
# removal_inds_i = np.concatenate((np.where(np.isnan(r_relaxes) == True)[0],np.where(np.isinf(r_relaxes) == True)[0]))
# removal_inds = np.concatenate((removal_inds_i, np.where(log_bh_masses<6)[0]))
# removal_inds = np.unique(removal_inds)
# 
# r_relaxes_nonans = np.delete(r_relaxes, removal_inds)
# log_bh_masses_nonans = np.delete(log_bh_masses, removal_inds)
# log_rhos_nonans = np.delete(log_rhos, removal_inds)
# gammas_nonans = np.delete(gammas, removal_inds)
# #%%
# 
# # calculate the percentage of r_relax above 5pc
# perc_above = len(np.where(r_relaxes_nonans>5)[0])/len(r_relaxes_nonans)
# 
# plt.figure(dpi=500)
# plt.title('Early = {} Gyr, Late = {} Gyr'.format(et_age,lt_age))
# a = plt.hist(r_relaxes_nonans, bins=50, range=(0,10))
# plt.xlabel('r$_{relax}$ [pc]')
# #plt.xlim(0,10)
# plt.plot([5,5],[0,np.max(a[0])],linestyle='--',color='r')
# plt.text(0.65, 0.75, '{:.2f}% > 5pc'.format(perc_above*100), transform=plt.gcf().transFigure, size=10)
# plt.savefig('../Plots/r_relax_investigation/r_relax_distribution_ET_{}_LT_{}.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# #%%
# n_outliers = 30
# 
# plt.figure(dpi=500)
# plt.plot(log_bh_masses_nonans,r_relaxes_nonans,linestyle='',marker='.')
# outlier_inds = u.get_n_max(r_relaxes_nonans,n_outliers)[1]
# plt.plot(log_bh_masses_nonans[outlier_inds],r_relaxes_nonans[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(log_bh_masses_nonans),np.max(log_bh_masses_nonans)],[5,5],linestyle='--',color='r')
# plt.xlabel('log(M$_{BH}$ [M$_\odot$])')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_BH_mass_ET_{}_LT_{}.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# plt.plot(log_rhos_nonans,r_relaxes_nonans,linestyle='',marker='.')
# plt.plot(log_rhos_nonans[outlier_inds],r_relaxes_nonans[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(log_rhos_nonans),np.max(log_rhos_nonans)],[5,5],linestyle='--',color='r')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/$pc^3$])')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_rho_ET_{}_LT_{}.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# plt.plot(gammas_nonans,r_relaxes_nonans,linestyle='',marker='.')
# plt.plot(gammas_nonans[outlier_inds],r_relaxes_nonans[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(gammas_nonans),np.max(gammas_nonans)],[5,5],linestyle='--',color='r')
# plt.xlabel('$\gamma$')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_gamma_ET_{}_LT_{}.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# #%%
# 
# r_relax_sl = np.delete(r_relaxes_nonans, np.where(gammas_nonans > -2.25))
# bh_masses_sl = np.delete(log_bh_masses_nonans, np.where(gammas_nonans > -2.25))
# gammas_sl = np.delete(gammas_nonans, np.where(gammas_nonans > -2.25))
# rhos_sl = np.delete(log_rhos_nonans, np.where(gammas_nonans > -2.25))
# 
# 
# # calculate the percentage of r_relax above 5pc
# perc_above = len(np.where(r_relax_sl>5)[0])/len(r_relax_sl)
# 
# plt.figure(dpi=500)
# plt.title('Early = {} Gyr, Late = {} Gyr, $\gamma$ < -2.25'.format(et_age,lt_age))
# a = plt.hist(r_relax_sl, bins=50, range=(0,10))
# plt.xlabel('r$_{relax}$ [pc]')
# plt.xlim(0,10)
# plt.plot([5,5],[0,np.max(a[0])],linestyle='--',color='r')
# plt.text(0.65, 0.75, '{:.2f}% > 5pc'.format(perc_above*100), transform=plt.gcf().transFigure, size=10)
# plt.savefig('../Plots/r_relax_investigation/r_relax_distribution_ET_{}_LT_{}_slopes_above_225.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# n_outliers = 30
# 
# plt.figure(dpi=500)
# plt.plot(bh_masses_sl,r_relax_sl,linestyle='',marker='.')
# outlier_inds = u.get_n_max(r_relax_sl,n_outliers)[1]
# plt.plot(bh_masses_sl[outlier_inds],r_relax_sl[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(bh_masses_sl),np.max(bh_masses_sl)],[5,5],linestyle='--',color='r')
# plt.xlabel('log(M$_{BH}$ [M$_\odot$])')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_BH_mass_ET_{}_LT_{}_slopes_above_225.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# plt.plot(rhos_sl,r_relax_sl,linestyle='',marker='.')
# plt.plot(rhos_sl[outlier_inds],r_relax_sl[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(rhos_sl),np.max(rhos_sl)],[5,5],linestyle='--',color='r')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/$pc^3$])')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_rho_ET_{}_LT_{}_slopes_above_225.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# plt.plot(gammas_sl,r_relax_sl,linestyle='',marker='.')
# plt.plot(gammas_sl[outlier_inds],r_relax_sl[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(gammas_sl),np.max(gammas_sl)],[5,5],linestyle='--',color='r')
# plt.xlabel('$\gamma$')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_gamma_ET_{}_LT_{}_slopes_above_225.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# 
# #%% 
# 
# # let's add in the code needed to perform the same investigation but for our 
# # real galaxies
# 
# # read in fits table with model galaxy parameters
# model_gal_filename = '../Result_Tables/all_gal_data_2x_pixel_scale_or_10pc_extinction_corr_nsa_ml_w_vi_color.fits'
# hdul = fits.open(model_gal_filename)
# gal_data = hdul[1].data 
# hdul.close()
# #%%
# types = gal_data['type']
# gammas = gal_data['gamma']
# rhos = 10**gal_data['rho5pc'] * M_sol_to_kg/pc_to_m**3
# log_rhos = gal_data['rho5pc']
# bh_masses = 10**gal_data['bhmass'] * M_sol_to_kg
# log_bh_masses = gal_data['bhmass']
# gal_masses = gal_data['galmass']
# 
# et_age = 10
# lt_age = 5
# t_ages = np.zeros(len(types))
# for i in range(len(t_ages)):
#     if types[i] == 0:
#         t_ages[i] = et_age*10**9*years_to_sec
#     else:
#         t_ages[i] = lt_age*10**9*years_to_sec
# 
# kappa = 1
# r5pc = 5*pc_to_m
# r_relaxes = np.zeros(len(types))
# for i in range(len(r_relaxes)):
#     r_relaxes[i] = ((G**2*rhos[i]*r5pc**(-gammas[i])*M_star*np.log(0.4*bh_masses[i]/M_star)*t_ages[i])/(kappa*(G*bh_masses[i])**(3/2)))**(1/(-gammas[i]-3/2))
#     r_relaxes[i] /= pc_to_m
# 
# 
# 
# 
# 
# =============================================================================


#%%

# =============================================================================
# OUR REAL GALAXIES
# =============================================================================


# read in fits table with model galaxy parameters
model_gal_filename = '../Result_Tables/final_gal_data.fits'
hdul = fits.open(model_gal_filename)
gal_data = hdul[1].data 
hdul.close()

names = gal_data['name']
RAs = gal_data['RA']
types = gal_data['type']
gammas = gal_data['slope']
rhos = 10**gal_data['cen_dens'] * M_sol_to_kg/pc_to_m**3
rhos_1 = 10**gal_data['cen_dens']
log_rhos = gal_data['cen_dens']
gal_masses = 10**gal_data['logmass']
dists = gal_data['dist']

#%%
##### BH Mass Assignment #####
BH_mass_ext = '_reines15'

# functions using relations from Reines et al. 2015
def get_BH_mass_lt(m_gal):
    mean = 7.45 + 1.05*np.log10(m_gal/10**11) # gives log(M_BH)
    scat = 0.55 #dex
    #return np.random.normal(loc=mean,scale=scat)
    return mean
def get_BH_mass_et(m_gal):
    mean = 8.95 + 1.40*np.log10(m_gal/10**11) # gives log(M_BH)
    scat = 0.47 #dex
    #return np.random.normal(loc=mean,scale=scat)
    return mean

print('Assigning BH masses...')

log_bh_masses = np.zeros_like(gammas)
for i in range(len(log_bh_masses)):
    if types[i] == 0:
        log_bh_masses[i] = get_BH_mass_et(gal_masses[i])
    else:
        log_bh_masses[i] = get_BH_mass_lt(gal_masses[i])
    #log_bh_masses[i] = 6
bh_masses = 10**(log_bh_masses) * M_sol_to_kg
#log_bh_masses = np.log10(bh_masses)

et_age = 1
lt_age = 1
t_ages = np.zeros(len(types))
for i in range(len(t_ages)):
    if types[i] == 0:
        t_ages[i] = et_age*10**9*years_to_sec
    else:
        t_ages[i] = lt_age*10**9*years_to_sec

kappa = 0.34 # from binney & tremaine
r5pc = 5*pc_to_m
r_relaxes = np.zeros(len(types))
r_relaxes_106 = np.zeros(len(types))
for i in range(len(r_relaxes)):
    r_relaxes[i] = ((G**2*rhos[i]*r5pc**(-gammas[i])*M_star*np.log(0.4*bh_masses[i]/M_star)*t_ages[i]*(1-gammas[i])**(3/2))/(kappa*(G*bh_masses[i])**(3/2)))**(1/(-gammas[i]-3/2))
    r_relaxes[i] /= pc_to_m
    r_relaxes_106[i] = ((G**2*rhos[i]*r5pc**(-gammas[i])*M_star*np.log(0.4*(10**6*M_sol_to_kg)/M_star)*t_ages[i]*(1-gammas[i])**(3/2))/(kappa*(G*(10**6*M_sol_to_kg))**(3/2)))**(1/(-gammas[i]-3/2))
    r_relaxes_106[i] /= pc_to_m

removal_inds_i = np.concatenate((np.where(np.isnan(r_relaxes) == True)[0],np.where(np.isinf(r_relaxes) == True)[0]))
removal_inds = np.concatenate((removal_inds_i, np.where(log_bh_masses<1)[0]))
removal_inds = np.unique(removal_inds)

r_relaxes_nonans = np.log10(np.delete(r_relaxes, removal_inds))
r_relaxes_106_nonans = np.log10(np.delete(r_relaxes_106, removal_inds))
log_bh_masses_nonans = np.delete(log_bh_masses, removal_inds)
rhos_nonans = np.delete(rhos, removal_inds)
rhos_1_nonans = np.delete(rhos_1, removal_inds)
log_rhos_nonans = np.delete(log_rhos, removal_inds)
gammas_nonans = np.delete(gammas, removal_inds)
names_nonans = np.delete(names, removal_inds)
gal_masses_nonans = np.delete(gal_masses, removal_inds)
RAs_nonans = np.delete(RAs, removal_inds)
dists_nonans = np.delete(dists, removal_inds)
gal_masses_nonans = np.delete(gal_masses, removal_inds)
#%%

high_inds = np.where(r_relaxes>5)[0]

# calculate the percentage of r_relax above 5pc
perc_above = len(r_relaxes[high_inds])/len(r_relaxes)

ymax = 200
ymin = -5

r_relaxes_temp = np.zeros_like(r_relaxes_nonans)
for i in range(len(r_relaxes_temp)):
    if r_relaxes_nonans[i] < -5:
        r_relaxes_temp[i] = -4.99
    else:
        r_relaxes_temp[i] = r_relaxes_nonans[i]


plt.figure(dpi=500)
plt.title('Early = {} Gyr, Late = {} Gyr'.format(et_age,lt_age))
a = plt.hist(r_relaxes_temp, bins=100, range=(ymin,ymax))
plt.xlabel('log(r$_{relax}$ [pc])')
#plt.xlim(0,10)
plt.plot([np.log10(5),np.log10(5)],[0,np.max(a[0])],linestyle='--',color='r')
plt.text(0.65, 0.75, '{:.2f}% > 5pc'.format(perc_above*100), transform=plt.gcf().transFigure, size=10)
plt.savefig('../Plots/r_relax_investigation/r_relax_distribution_ET_{}_LT_{}.png'.format(et_age,lt_age),
                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#%%
# =============================================================================
# n_outliers = 1
# 
# plt.figure(dpi=500)
# plt.plot(log_bh_masses_nonans,r_relaxes_nonans,linestyle='',marker='.')
# outlier_inds = high_inds
# #= u.get_n_max(r_relaxes_nonans,n_outliers)[1]
# plt.plot(log_bh_masses_nonans[outlier_inds],r_relaxes_nonans[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(log_bh_masses_nonans),np.max(log_bh_masses_nonans)],[5,5],linestyle='--',color='r')
# plt.xlabel('log(M$_{BH}$ [M$_\odot$])')
# #plt.ylabel('r$_{relax}$ [pc]')
# plt.ylabel('log(r$_{relax}$ [pc])')
# plt.ylim(0,ymax)
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_BH_mass_ET_{}_LT_{}.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# plt.plot(log_rhos_nonans,r_relaxes_nonans,linestyle='',marker='.')
# plt.plot(log_rhos_nonans[outlier_inds],r_relaxes_nonans[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(log_rhos_nonans),np.max(log_rhos_nonans)],[5,5],linestyle='--',color='r')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/$pc^3$])')
# #plt.ylabel('r$_{relax}$ [pc]')
# plt.ylabel('log(r$_{relax}$ [pc])')
# plt.ylim(0,ymax)
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_rho_ET_{}_LT_{}.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# plt.plot(gammas_nonans,r_relaxes_nonans,linestyle='',marker='.')
# plt.plot(gammas_nonans[outlier_inds],r_relaxes_nonans[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(gammas_nonans),np.max(gammas_nonans)],[5,5],linestyle='--',color='r')
# plt.xlabel('$\gamma$')
# #plt.ylabel('r$_{relax}$ [pc]')
# plt.ylabel('log(r$_{relax}$ [pc])')
# plt.ylim(0,ymax)
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_gamma_ET_{}_LT_{}.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# =============================================================================

#%%

fig, ax = plt.subplots(dpi=500)
ax.plot(log_bh_masses_nonans,r_relaxes_nonans,linestyle='',marker='.')
outlier_inds = high_inds
steep_outlier_inds = np.intersect1d(outlier_inds, np.where(gammas<=-1.75))

gamma_inds = np.where(gammas<=-1.75)[0]
gamma_inds = gamma_inds[16:] # remove pechetti gals for now
#gamma_inds = gamma_inds[-10:] # just the Hoyer+22 sample


#pdb.set_trace()
ax.plot(log_bh_masses_nonans[steep_outlier_inds],r_relaxes_nonans[steep_outlier_inds],linestyle='',marker='.')
ax.plot([np.min(log_bh_masses_nonans),np.max(log_bh_masses_nonans)],[np.log10(5),np.log10(5)],linestyle='--',color='r')
ax.set_xlabel('log(M$_{BH}$ [M$_\odot$])')
plt.ylabel('log(r$_{relax}$ [pc])')
#ax.set_ylim(ymin,ymax)
ax.set_ylim(-3,10)
ax2 = ax.twinx()
ax2.hist(log_bh_masses_nonans, bins=100, color='0.6', alpha=0.2)
ax2.set_ylabel('# of Galaxies')
plt.savefig('../Plots/r_relax_investigation/r_relax_v_BH_mass_ET_{}_LT_{}_w_hist.png'.format(et_age,lt_age),
                 bbox_inches='tight', pad_inches=0.1, dpi=500)

fig, ax = plt.subplots(dpi=500)
ax.plot(log_rhos,r_relaxes_nonans,linestyle='',marker='.')
ax.plot(log_rhos[steep_outlier_inds],r_relaxes_nonans[steep_outlier_inds],linestyle='',marker='.')
ax.plot([np.min(log_rhos),np.max(log_rhos)],[np.log10(5),np.log10(5)],linestyle='--',color='r')
ax.set_xlabel('log($\\rho_{5pc}$ [M$_\odot$/$pc^3$])')
plt.ylabel('log(r$_{relax}$ [pc])')
#ax.set_ylim(ymin,ymax)
ax.set_ylim(-3,10)
ax2 = ax.twinx()
ax2.hist(log_rhos, bins=100, color='0.6', alpha=0.2)
ax2.set_ylabel('# of Galaxies')
plt.savefig('../Plots/r_relax_investigation/r_relax_v_rho_ET_{}_LT_{}_w_hist.png'.format(et_age,lt_age),
                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#%%

fig, ax = plt.subplots(dpi=500)
ax.plot(gammas,r_relaxes_nonans,linestyle='',marker='.')
ax.plot(gammas[steep_outlier_inds],r_relaxes_nonans[steep_outlier_inds],linestyle='',marker='.')
ax.plot([np.min(gammas),np.max(gammas)],[np.log10(5),np.log10(5)],linestyle='--',color='r')
ax.set_xlabel('$\gamma$')
plt.ylabel('log(r$_{relax}$ [pc])')
ax.set_ylim(ymin,ymax)
#ax.set_ylim(-3,10)
ax.set_xlim(-4.1,-1.75)
ax2 = ax.twinx()
ax2.hist(gammas, bins=100, color='0.6', alpha=0.2)
ax2.set_ylabel('# of Galaxies')
plt.savefig('../Plots/r_relax_investigation/r_relax_v_gamma_ET_{}_LT_{}_w_hist.png'.format(et_age,lt_age),
                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#%%
fig, ax = plt.subplots(dpi=500)
ax.plot(np.log10(gal_masses_nonans),r_relaxes_nonans,linestyle='',marker='.')
ax.plot(np.log10(gal_masses_nonans[steep_outlier_inds]),r_relaxes_nonans[steep_outlier_inds],linestyle='',marker='.')
ax.plot([np.min(np.log10(gal_masses_nonans)),np.max(np.log10(gal_masses_nonans))],[np.log10(5),np.log10(5)],linestyle='--',color='r')
ax.set_xlabel('log(M$_{gal} [M_\odot]$)')
plt.ylabel('log(r$_{relax}$ [pc])')
#ax.set_ylim(ymin,ymax)
ax.set_ylim(-3,10)
#ax.set_xlim(-4.1,-1.75)
ax2 = ax.twinx()
ax2.hist(np.log10(gal_masses_nonans), bins=100, color='0.6', alpha=0.2)
ax2.set_ylabel('# of Galaxies')
plt.savefig('../Plots/r_relax_investigation/r_relax_v_galmass_ET_{}_LT_{}_w_hist.png'.format(et_age,lt_age),
                 bbox_inches='tight', pad_inches=0.1, dpi=500)

#%%

# =============================================================================
# ol_MGEs = MGEs_nonans[:,steep_outlier_inds]
# ol_lograd = lograd_nonans[:,steep_outlier_inds]
# ol_dists = dists_nonans[:,steep_outlier_inds]
# 
# res_limit = lograd_nonans[0,:]
# =============================================================================
slope_ext = '2x_pixel_scale'
phys_ext = '_or_10pc_extinction_corr_nsa_ml_w_vi_color'
gal_file = '../Result_Tables/all_gal_data_'+slope_ext+phys_ext+'.fits'
hdul = fits.open(gal_file)
head = hdul[0].data
dat = hdul[1].data
hdul.close()

names_1 = dat['name']
RAs_1 = dat['RA']
MGEs = dat['MGEs']
lograd = dat['lograd']
logdens = dat['logdens']

r_relaxes_modified = r_relaxes_nonans[22:]
r_relaxes_106_modified = r_relaxes_106_nonans[22:]
RA_nonans_modified = RAs_nonans[22:]
log_bh_masses_modified = log_bh_masses_nonans[22:]
rhos_modified = rhos_1_nonans[22:]
rhos_SI = rhos_nonans[22:]
gammas_modified = gammas_nonans[22:]
gal_masses_modified = np.log10(gal_masses_nonans[22:])
sort_inds = np.zeros_like(RA_nonans_modified).astype(int)
for i in range(len(RA_nonans_modified)):
    sort_inds[i] = u.find_nearest(RAs_1, RA_nonans_modified[i])

    
MGEs_sort = MGEs#[sort_inds,:]
lograd_sort = lograd#[sort_inds,:]
RAs_1_sort = RAs_1#[sort_inds]
logdens_sort = logdens#[sort_inds,:]

lograd_interp = np.log10(np.geomspace(0.1, 10**6, 2000))
logdens_interp = np.zeros((len(lograd_sort[:,0]),2000))
#fivepc_ind = u.find_nearest(lograd_interp, np.log10())
for i in range(len(lograd_sort[:,0])):
    logdens_interp[i,:] = np.interp(lograd_interp, lograd_sort[i,:], logdens_sort[i,:])
    # adjust the radii inward of 5pc to be the inner power law fit
    inner_ind = u.find_nearest(lograd_interp, np.log10(5))
    logdens_interp[i,0:inner_ind] = np.log10(rhos_modified[i]*((10**lograd_interp[0:inner_ind])/5)**(gammas_modified[i]))
# =============================================================================
#     plt.figure(dpi=500)
#     plt.plot(lograd_sort[i,:],logdens_sort[i,:])
#     plt.plot(lograd_interp, logdens_interp[i,:])
#     #pdb.set_trace()
# =============================================================================

# galaxy order: ['bts76', 'ddo084', 'kk2000-03', 'kk2000-53', 'kk96', 'leg09',
#                'lvj1217+47', 'ngc5011c', 'pgc4310323', 'ugc07242']
R_e_hoyer = np.array([33.16, 3.81, 3.75, 6.58, 6.05, 7.03, 4.74, 15.00, 4.44, 2.772])
for i in range(len(R_e_hoyer)):
    logRe = np.log10(R_e_hoyer[i])
    ind = u.find_nearest(lograd_interp, logRe)
    logdens_interp[-len(R_e_hoyer)+i, ind:] = np.ones_like(logdens_interp[-len(R_e_hoyer)+i, ind:])*logdens_interp[-len(R_e_hoyer)+i, ind]
# =============================================================================
# lograd_interp_m5 = np.log10(np.geomspace(1e-5, 10**3, 1000))
# logdens_interp_m5 = np.zeros((len(lograd_sort[:,0]),1000))
# #fivepc_ind = u.find_nearest(lograd_interp, np.log10())
# for i in range(len(lograd_sort[:,0])):
#     logdens_interp_m5[i,:] = np.interp(lograd_interp, lograd_sort[i,:], logdens_sort[i,:])
#     # adjust the radii inward of 5pc to be the inner power law fit
#     inner_ind = u.find_nearest(lograd_interp, np.log10(10))
#     logdens_interp_m5[i,0:inner_ind+1] = np.log10(rhos_modified[i]*(10**lograd_interp[0:inner_ind+1]/5)**(gammas_modified[i]))
# =============================================================================


# =============================================================================
# # use the interpolated density profiles to include M_enc in the r_relax computation
# r_relaxes = np.zeros(len(types))
# for i in range(len(r_relaxes)):
#     r_relaxes[i] = ((G**2*rhos[i]*r5pc**(-gammas[i])*M_star*np.log(0.4*bh_masses[i]/M_star)*t_ages[i])/(kappa*(G*bh_masses[i])**(3/2)))**(1/(-gammas[i]-3/2))
#     r_relaxes[i] /= pc_to_m
# =============================================================================




res_limit = lograd_sort[:,0]

r_relaxes_steep_outlier = r_relaxes_modified[steep_outlier_inds[5:]-22]
res_limit_steep_outlier = res_limit[steep_outlier_inds[5:]-22]
log_bh_masses_steep_outlier = log_bh_masses_modified[steep_outlier_inds[5:]-22]
lograd_steep_outlier = lograd_sort[steep_outlier_inds[5:]-22,:]
logdens_steep_outlier_interp = logdens_interp[steep_outlier_inds[5:]-22,:]
logdens_steep_outlier = logdens_sort[steep_outlier_inds[5:]-22,:]

r_relaxes_gamma = r_relaxes_modified[gamma_inds-22]
r_relaxes_106_gamma = r_relaxes_106_modified[gamma_inds-22]
res_limit_gamma = res_limit[gamma_inds-22]
log_bh_masses_gamma = log_bh_masses_modified[gamma_inds-22]
lograd_gamma = lograd_sort[gamma_inds-22,:]
logdens_gamma_interp = logdens_interp[gamma_inds-22,:]
logdens_gamma = logdens_sort[gamma_inds-22,:]
gammas_gamma = gammas_modified[gamma_inds-22]
rhos_gamma = rhos_modified[gamma_inds-22]
gal_masses_gamma = gal_masses_modified[gamma_inds-22]
# =============================================================================
# soi_steep_outlier = np.zeros(len(lograd_steep_outlier[:,0]))
# for j in range(len(lograd_steep_outlier[:,0])):
#     for i in range(len(lograd_steep_outlier[0,:])-1):
#         M_enc = get_enc_mass(10**lograd_steep_outlier[j,:], 10**logdens_steep_outlier[j,:], i+1)
#         if M_enc >= 2*10**log_bh_masses_steep_outlier[j]:
#             soi_steep_outlier[j] = 10**lograd_steep_outlier[j,i+1]
#             break
#         
# soi_steep_outlier_106 = np.zeros(len(lograd_steep_outlier[:,0]))
# for j in range(len(lograd_steep_outlier[:,0])):
#     for i in range(len(lograd_steep_outlier[0,:])-1):
#         M_enc = get_enc_mass(10**lograd_steep_outlier[j,:], 10**logdens_steep_outlier[j,:], i+1)
#         if M_enc >= 2*10**6:
#             soi_steep_outlier_106[j] = 10**lograd_steep_outlier[j,i+1]
#             break
#         else:
#             soi_steep_outlier_106[j] = 10**lograd_steep_outlier[j,-1]
# =============================================================================
        

# =============================================================================
# soi = np.zeros(len(lograd_sort[:,0]))
# for j in range(len(lograd_sort[:,0])):
#     for i in range(len(lograd_sort[0,:])-1):
#         M_enc = get_enc_mass(10**lograd_sort[j,:], 10**logdens_sort[j,:], i+1)
#         if M_enc >= 2*10**log_bh_masses_modified[j]:
#             soi[j] = 10**lograd_sort[j,i+1]
#             break
# =============================================================================

soi = np.zeros(len(lograd_sort[:,0]))
for j in range(len(lograd_sort[:,0])):
    for i in range(len(lograd_interp)-1):
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_interp[j,:], i+1)
        if M_enc >= 2*10**log_bh_masses_modified[j]:
            soi[j] = 10**lograd_interp[i+1]
            break

soi_steep_outlier = np.zeros(len(lograd_steep_outlier[:,0]))
for j in range(len(lograd_steep_outlier[:,0])):
    for i in range(len(lograd_interp)-1):
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_steep_outlier_interp[j,:], i+1)
        if M_enc >= 2*10**log_bh_masses_steep_outlier[j]:
            soi_steep_outlier[j] = 10**lograd_interp[i+1]
            break

soi_gamma = np.zeros(len(lograd_gamma[:,0]))
for j in range(len(lograd_gamma[:,0])):
    for i in range(len(lograd_interp)-1):
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_gamma_interp[j,:], i+1)
        if M_enc >= 2*10**log_bh_masses_gamma[j]:
            soi_gamma[j] = 10**lograd_interp[i+1]
            break


soi_106 = np.zeros(len(lograd_sort[:,0]))
for j in range(len(lograd_sort[:,0])):
    for i in range(len(lograd_interp)-1):
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_interp[j,:], i+1)
        if M_enc >= 2*10**6:
            soi_106[j] = 10**lograd_interp[i+1]
            break
        else:
            soi_106[j] = 10**lograd_interp[-1]
    
soi_steep_outlier_106 = np.zeros(len(lograd_steep_outlier[:,0]))
for j in range(len(lograd_steep_outlier[:,0])):
    for i in range(len(lograd_interp)-1):
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_steep_outlier_interp[j,:], i+1)
        if M_enc >= 2*10**6:
            soi_steep_outlier_106[j] = 10**lograd_interp[i+1]
            break
        else:
            soi_steep_outlier_106[j] = 10**lograd_interp[-1] 

soi_gamma_106 = np.zeros(len(lograd_gamma[:,0]))
for j in range(len(lograd_gamma[:,0])):
    for i in range(len(lograd_interp)-1):
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_gamma_interp[j,:], i+1)
        if M_enc >= 2*10**6:
            soi_gamma_106[j] = 10**lograd_interp[i+1]
            break
        else:
            soi_gamma_106[j] = 10**lograd_interp[-1] 
      
soi_resids_gam = 10**res_limit_gamma - soi_gamma         
resolved_gam_inds = np.where(soi_resids_gam < 0)[0]      

# =============================================================================
# soi_106 = np.zeros(len(lograd_sort[:,0]))
# for j in range(len(lograd_sort[:,0])):
#     for i in range(len(lograd_sort[0,:])-1):
#         M_enc = get_enc_mass(10**lograd_sort[j,:], 10**logdens_sort[j,:], i+1)
#         if M_enc >= 2*10**6:
#             soi_106[j] = 10**lograd_sort[j,i+1]
#             break
#         else:
#             soi_106[j] = 10**lograd_sort[j,-1]
# 
# =============================================================================

# =============================================================================
# soi_1 = np.zeros_like(log_bh_masses_modified)
# soi_1 = (((3-gammas_modified)*10**(log_bh_masses_modified))/(2*np.pi*rhos_modified*5**(gammas_modified)))**(1/(3-gammas_modified))
# 
# =============================================================================
# =============================================================================
# soi_106 = np.zeros_like(log_bh_masses_modified)
# soi_106 = (((3-gammas_modified)*10**(6))/(2*np.pi*rhos_modified*5**(gammas_modified)))**(1/(3-gammas_modified))
# =============================================================================
# =============================================================================
# for i in range(len(soi)):
#     for j in range(len(lograd_nonans[:,i])-1):
#         M_enc = get_enc_mass(lograd_nonans[:,i],gammas_nonans[i],rhos_nonans[i],j+1)
#         if M_enc
# 
# =============================================================================
#%%
# let's plot the sphere of influence vs the resolution limit of our galaxies
plt.figure(dpi=500)
plt.plot(np.log10(soi), res_limit, '.b')
#plt.plot(np.log10(soi_steep_outlier), res_limit_steep_outlier, '.r')
plt.plot(np.log10(soi_gamma), res_limit_gamma, '.r')
#plt.plot([0.2,8],[0.2,8],'--k')
plt.plot([-0.7,2.6],[-0.7,2.6],'--k')
plt.xlabel('log(Sphere of Influence Radius [pc])')
plt.ylabel('log(Resolution Limit Radius [pc])')
#plt.xlim(-0.3,1.5)

#%%

# let's plot the sphere of influence vs the resolution limit of our galaxies
plt.figure(dpi=500)
plt.scatter(np.log10(soi), res_limit, s=15, color='0.9')
#plt.plot(np.log10(soi_steep_outlier), res_limit_steep_outlier, '.r')
#plt.plot(np.log10(soi_gamma), res_limit_gamma, '.r')
plt.scatter(np.log10(soi_gamma), res_limit_gamma, s=15, c=log_bh_masses_gamma, cmap='viridis')
#plt.scatter(np.log10(soi_gamma[resolved_gam_inds]), res_limit_gamma[resolved_gam_inds], s=15, c=log_bh_masses_modified[resolved_gam_inds], cmap='viridis')

#plt.plot(np.log10(soi_gamma[resolved_gam_inds]), res_limit_gamma[resolved_gam_inds], '.b')
plt.plot([-0.7,2.6],[-0.7,2.6],'--k')
#plt.plot([np.min([np.min(np.log10(soi_106)),np.min(res_limit)]),np.min([np.min(np.log10(soi_106)),np.min(res_limit)])], \
#         [np.max([np.max(np.log10(soi_106)),np.max(res_limit)]),np.max([np.max(np.log10(soi_106)),np.max(res_limit)])],'--k')
plt.xlabel('log(Sphere of Influence Radius [pc])')
plt.ylabel('log(Resolution Limit Radius [pc])')
plt.colorbar(label='log(M$_{BH}$ [M$_\odot$])')
#plt.xlim(-0.3,1.5)



# let's plot the sphere of influence vs the resolution limit of our galaxies
plt.figure(dpi=500)
plt.scatter(np.log10(soi), res_limit, s=15, color='0.9')
#plt.plot(np.log10(soi_steep_outlier), res_limit_steep_outlier, '.r')
#plt.plot(np.log10(soi_gamma), res_limit_gamma, '.r')
plt.scatter(np.log10(soi_gamma), res_limit_gamma, s=15, c=np.log10(rhos_gamma), cmap='viridis')
#plt.scatter(np.log10(soi_gamma[resolved_gam_inds]), res_limit_gamma[resolved_gam_inds], s=15, c=log_bh_masses_modified[resolved_gam_inds], cmap='viridis')

#plt.plot(np.log10(soi_gamma[resolved_gam_inds]), res_limit_gamma[resolved_gam_inds], '.b')
plt.plot([-0.7,2.6],[-0.7,2.6],'--k')
#plt.plot([np.min([np.min(np.log10(soi_106)),np.min(res_limit)]),np.min([np.min(np.log10(soi_106)),np.min(res_limit)])], \
#         [np.max([np.max(np.log10(soi_106)),np.max(res_limit)]),np.max([np.max(np.log10(soi_106)),np.max(res_limit)])],'--k')
plt.xlabel('log(Sphere of Influence Radius [pc])')
plt.ylabel('log(Resolution Limit Radius [pc])')
plt.colorbar(label='log($\\rho$ [M$_\odot$]/pc$^3$])')
#plt.xlim(-0.3,1.5)



# let's plot the sphere of influence vs the resolution limit of our galaxies
plt.figure(dpi=500)
plt.scatter(np.log10(soi), res_limit, s=15, color='0.9')
#plt.plot(np.log10(soi_steep_outlier), res_limit_steep_outlier, '.r')
#plt.plot(np.log10(soi_gamma), res_limit_gamma, '.r')
plt.scatter(np.log10(soi_gamma), res_limit_gamma, s=15, c=gammas_gamma, cmap='viridis')
#plt.scatter(np.log10(soi_gamma[resolved_gam_inds]), res_limit_gamma[resolved_gam_inds], s=15, c=log_bh_masses_modified[resolved_gam_inds], cmap='viridis')

#plt.plot(np.log10(soi_gamma[resolved_gam_inds]), res_limit_gamma[resolved_gam_inds], '.b')
plt.plot([-0.7,2.6],[-0.7,2.6],'--k')
#plt.plot([np.min([np.min(np.log10(soi_106)),np.min(res_limit)]),np.min([np.min(np.log10(soi_106)),np.min(res_limit)])], \
#         [np.max([np.max(np.log10(soi_106)),np.max(res_limit)]),np.max([np.max(np.log10(soi_106)),np.max(res_limit)])],'--k')
plt.xlabel('log(Sphere of Influence Radius [pc])')
plt.ylabel('log(Resolution Limit Radius [pc])')
plt.colorbar(label='$\gamma$')
#plt.xlim(-0.3,1.5)


#%%
plt.figure(dpi=500)
plt.scatter(gammas_gamma, np.log10(rhos_gamma), s=15, color='0.9')
plt.scatter(gammas_gamma[resolved_gam_inds], np.log10(rhos_gamma[resolved_gam_inds]), c=log_bh_masses_gamma[resolved_gam_inds], s=15, cmap='viridis')
plt.xlabel('$\gamma$')
plt.ylabel('log($\\rho$ [M$_\odot$]/pc$^3$])')
plt.colorbar(label='log(M$_{BH}$ [M$_\odot$])')


#%%


plt.figure(dpi=500)
#plt.scatter(res_limit,r_relaxes_modified)
plt.scatter(res_limit_gamma,r_relaxes_gamma,s=15,c=log_bh_masses_gamma)
#plt.scatter(r_relaxes_steep_outlier,res_limit_steep_outlier,s=15)
plt.xlabel('log(Resolution Limit [pc])')
plt.ylabel('log(r$_{relax}$ [pc])')
plt.plot([-0.4,6],[-0.4,6],'--k')
plt.ylim(-8,7)
plt.xlim(-0.45,0.65)
plt.colorbar(label='log(M$_{BH}$ [M$_\odot$])')

plt.figure(dpi=500)
#plt.scatter(res_limit,r_relaxes_modified)
plt.scatter(res_limit_gamma,r_relaxes_106_gamma,s=15)
#plt.scatter(r_relaxes_steep_outlier,res_limit_steep_outlier,s=15)
plt.xlabel('log(Resolution Limit [pc])')
plt.ylabel('log(r$_{relax}$ [pc])')
plt.plot([-0.4,6],[-0.4,6],'--k')
plt.ylim(-8,7)
plt.xlim(-0.45,0.65)
#plt.colorbar(label='log(M$_{BH}$ [M$_\odot$])')


#%%
plt.figure(dpi=500)
for i in range(len(logdens_sort[:,0])):
    plt.plot(lograd_interp, logdens_interp[i,:], '0.7')

plt.plot(lograd_interp, logdens_steep_outlier_interp[0,:], 'b', label='r_relax outliers')
for i in range(len(logdens_steep_outlier[:,0])):
    plt.plot(lograd_interp, logdens_steep_outlier_interp[i,:], 'b')

# =============================================================================
# plt.plot(lograd_interp, logdens_gamma_interp[0,:], 'b', label='$\gamma$ < -1.75')
# for i in range(len(logdens_gamma[:,0])):
#     plt.plot(lograd_interp, logdens_gamma_interp[i,:], 'b')
# =============================================================================

#plt.plot([0.2,8],[0.2,8],'--k')
#plt.plot([-0.7,2.6],[-0.7,2.6],'--k')
#plt.xlabel('log(Sphere of Influence Radius [pc])')
#plt.ylabel('log(Resolution Limit Radius [pc])')
#plt.xlim(-1,2)
#plt.ylim(0,5)
#plt.ylim(-7,5.5)
#plt.plot(lograd_interp, logdens_interp[-8,:], 'r')
#plt.plot(lograd_interp, logdens_interp[-7,:], 'r')
plt.ylabel('log($\\rho$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(Radius [pc])')
plt.legend()


# =============================================================================
# # let's plot the sphere of influence vs the resolution limit of our galaxies
# plt.figure(dpi=500)
# plt.plot(np.log10(soi_1), res_limit, '.b')
# #plt.plot([0.2,8],[0.2,8],'--k')
# plt.xlabel('log(Sphere of Influence Radius [pc])')
# plt.ylabel('log(Resolution Limit Radius [pc])')
# #plt.xlim(3,9)
# =============================================================================

#%%
# let's plot the sphere of influence vs the resolution limit of our galaxies
plt.figure(dpi=500)
plt.plot(np.log10(soi_106), res_limit, '.k')
#plt.plot(np.log10(soi_steep_outlier_106), res_limit_steep_outlier, '.r')
plt.plot(np.log10(soi_gamma_106), res_limit_gamma, '.r')
plt.plot(np.log10(soi_gamma_106[resolved_gam_inds]), res_limit_gamma[resolved_gam_inds], '.b')
plt.plot([-0.7,2.6],[-0.7,2.6],'--k')
#plt.plot([np.min([np.min(np.log10(soi_106)),np.min(res_limit)]),np.min([np.min(np.log10(soi_106)),np.min(res_limit)])], \
#         [np.max([np.max(np.log10(soi_106)),np.max(res_limit)]),np.max([np.max(np.log10(soi_106)),np.max(res_limit)])],'--k')
plt.xlabel('log(Sphere of Influence Radius [pc])')
plt.ylabel('log(Resolution Limit Radius [pc])')
plt.xlim(-1.3,6.3)

#%%

plt.figure(dpi=500)
a = plt.hist(log_bh_masses_modified, bins=30, color='0.8')
plt.hist(log_bh_masses_gamma, bins = a[1])
plt.hist(log_bh_masses_gamma[resolved_gam_inds], bins = a[1], alpha=0.8)
plt.xlabel('log(M$_{BH}$ [M$_\odot$])')

plt.figure(dpi=500)
a = plt.hist(gammas_modified, bins=30, color='0.8')
plt.hist(gammas_gamma, bins = a[1])
plt.hist(gammas_gamma[resolved_gam_inds], bins = a[1], alpha=0.8)
plt.xlabel('$\gamma$')

plt.figure(dpi=500)
a = plt.hist(np.log10(rhos_modified), bins=30, color='0.8')
plt.hist(np.log10(rhos_gamma), bins = a[1])
plt.hist(np.log10(rhos_gamma[resolved_gam_inds]), bins = a[1], alpha=0.8)
plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')

#%%

plt.figure(dpi=500)
plt.scatter(log_bh_masses_modified, gammas_modified, color='0.8', s=5)
plt.scatter(log_bh_masses_gamma, gammas_gamma, color='r', s=5)
plt.scatter(log_bh_masses_gamma[resolved_gam_inds], gammas_gamma[resolved_gam_inds], color='b', s=5)
plt.xlabel('log(M$_{BH}$ [M$_\odot$])')
plt.ylabel('$\gamma$')

plt.figure(dpi=500)
plt.scatter(log_bh_masses_modified, np.log10(rhos_modified), color='0.8', s=5)
plt.scatter(log_bh_masses_gamma, np.log10(rhos_gamma), color='r', s=5)
plt.scatter(log_bh_masses_gamma[resolved_gam_inds], np.log10(rhos_gamma[resolved_gam_inds]), color='b', s=5)
plt.xlabel('log(M$_{BH}$ [M$_\odot$])')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')

plt.figure(dpi=500)
plt.scatter(gammas_modified, np.log10(rhos_modified), color='0.8', s=5)
plt.scatter(gammas_gamma, np.log10(rhos_gamma), color='r', s=5)
plt.scatter(gammas_gamma[resolved_gam_inds], np.log10(rhos_gamma[resolved_gam_inds]), color='b', s=5)
plt.xlabel('$\gamma$')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')

#%%
plt.figure(dpi=500)
plt.scatter(gal_masses_modified, gammas_modified, color='0.8', s=5)
plt.scatter(gal_masses_gamma, gammas_gamma, color='r', s=5)
plt.scatter(gal_masses_gamma[resolved_gam_inds], gammas_gamma[resolved_gam_inds], color='b', s=5)
plt.xlabel('log(M$_{gal}$ [M$_\odot$])')
plt.ylabel('$\gamma$')

plt.figure(dpi=500)
plt.scatter(gal_masses_modified, np.log10(rhos_modified), color='0.8', s=5)
plt.scatter(gal_masses_gamma, np.log10(rhos_gamma), color='r', s=5)
plt.scatter(gal_masses_gamma[resolved_gam_inds], np.log10(rhos_gamma[resolved_gam_inds]), color='b', s=5)
plt.xlabel('log(M$_{gal}$ [M$_\odot$])')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')


#%%
# let's plot the relaxation time as f(r) for all galaxies.
relax_time_r_all = np.zeros_like(lograd_sort)
for j in range(len(relax_time_r_all[:,0])):
    for i in range(len(relax_time_r_all[0,:])):
        M_BH = 10**log_bh_masses_modified[j] * M_sol_to_kg
        M_enc = get_enc_mass(10**lograd_sort[j,:], 10**logdens_sort[j,:], i) * M_sol_to_kg
        #M_enc = 0
        relax_time_r_all[j,i] = kappa*((np.sqrt(G*(M_BH+M_enc)/(10**lograd_sort[j,i]*pc_to_m)))**3)/(G**2*10**logdens_sort[j,i]*(M_sol_to_kg/pc_to_m**3)*M_star*np.log(0.4*M_BH/M_star))
        relax_time_r_all[j,i] /= years_to_sec


# let's plot the relaxation time as f(r) for the outlier galaxies.
relax_time_r = np.zeros_like(lograd_steep_outlier)
for j in range(len(relax_time_r[:,0])):
    for i in range(len(relax_time_r[0,:])):
        M_BH = 10**log_bh_masses_steep_outlier[j] * M_sol_to_kg
        M_enc = get_enc_mass(10**lograd_steep_outlier[j,:], 10**logdens_steep_outlier[j,:], i) * M_sol_to_kg
        #M_enc = 0
        relax_time_r[j,i] = kappa*((np.sqrt(G*(M_BH+M_enc)/(10**lograd_steep_outlier[j,i]*pc_to_m)))**3)/(G**2*10**logdens_steep_outlier[j,i]*(M_sol_to_kg/pc_to_m**3)*M_star*np.log(0.4*M_BH/M_star))
        relax_time_r[j,i] /= years_to_sec

plt.figure(dpi=500)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(t$_{relax}$ [yr])')
for i in range(len(relax_time_r_all[:,0])):
    #plt.figure(dpi=500)
    plt.plot(lograd_sort[i,:], np.log10(relax_time_r_all[i,:]), '0.5')
plt.plot(lograd_steep_outlier[0,:], np.log10(relax_time_r[0,:]), 'b', label='$\gamma$ < -1.75')
for i in range(len(relax_time_r[:,0])):
    #plt.figure(dpi=500)
    plt.plot(lograd_steep_outlier[i,:], np.log10(relax_time_r[i,:]), 'b')
plt.legend()
#plt.ylim(6,10.13)
plt.ylim(6,20)

#%%


# let's plot the relaxation time as f(r) for all galaxies.
relax_time_r_all_106 = np.zeros_like(lograd_sort)
for j in range(len(relax_time_r_all_106[:,0])):
    for i in range(len(relax_time_r_all_106[0,:])):
        M_BH = 10**6 * M_sol_to_kg
        #M_enc = get_enc_mass(10**lograd_sort[j,:], 10**logdens_sort[j,:], i) * M_sol_to_kg
        M_enc = 0
        relax_time_r_all_106[j,i] = kappa*((np.sqrt(G*(M_BH+M_enc)/(10**lograd_sort[j,i]*pc_to_m)))**3)/(G**2*10**logdens_sort[j,i]*(M_sol_to_kg/pc_to_m**3)*M_star*np.log(0.4*M_BH/M_star))
        relax_time_r_all_106[j,i] /= years_to_sec


# let's plot the relaxation time as f(r) for the outlier galaxies.
relax_time_r_106 = np.zeros_like(lograd_steep_outlier)
for j in range(len(relax_time_r_106[:,0])):
    for i in range(len(relax_time_r_106[0,:])):
        M_BH = 10**6 * M_sol_to_kg
        #M_enc = get_enc_mass(10**lograd_steep_outlier[j,:], 10**logdens_steep_outlier[j,:], i) * M_sol_to_kg
        M_enc = 0
        relax_time_r_106[j,i] = kappa*((np.sqrt(G*(M_BH+M_enc)/(10**lograd_steep_outlier[j,i]*pc_to_m)))**3)/(G**2*10**logdens_steep_outlier[j,i]*(M_sol_to_kg/pc_to_m**3)*M_star*np.log(0.4*M_BH/M_star))
        relax_time_r_106[j,i] /= years_to_sec

plt.figure(dpi=500)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(t$_{relax}$ [yr])')
for i in range(len(relax_time_r_all_106[:,0])):
    #plt.figure(dpi=500)
    plt.plot(lograd_sort[i,:], np.log10(relax_time_r_all_106[i,:]), '0.5')
plt.plot(lograd_steep_outlier[0,:], np.log10(relax_time_r_106[0,:]), 'b', label='$\gamma$ < -1.75')
for i in range(len(relax_time_r_106[:,0])):
    #plt.figure(dpi=500)
    plt.plot(lograd_steep_outlier[i,:], np.log10(relax_time_r_106[i,:]), 'b')
plt.legend()

#plt.ylim(6,10.13)
plt.ylim(6,20)

#%%

# DO THE SAME AS ABOVE BUT WITH INTERPOLATED DATA
# let's plot the relaxation time as f(r) for all galaxies.
relax_time_r_all = np.zeros((len(lograd_sort[:,0]),len(lograd_interp)))
for j in range(len(relax_time_r_all[:,0])):
    for i in range(len(relax_time_r_all[0,:])):
        M_BH = 10**log_bh_masses_modified[j] * M_sol_to_kg
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_interp[j,:], i) * M_sol_to_kg
        relax_time_r_all[j,i] = kappa*((np.sqrt(G*(M_BH+M_enc)/(10**lograd_interp[i]*pc_to_m)))**3)/(G**2*10**logdens_interp[j,i]*(M_sol_to_kg/pc_to_m**3)*M_star*np.log(0.4*M_BH/M_star))
        relax_time_r_all[j,i] /= years_to_sec


# let's plot the relaxation time as f(r) for the outlier galaxies.
relax_time_r = np.zeros((len(lograd_steep_outlier[:,0]),len(lograd_interp)))
for j in range(len(relax_time_r[:,0])):
    for i in range(len(relax_time_r[0,:])):
        M_BH = 10**log_bh_masses_steep_outlier[j] * M_sol_to_kg
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_steep_outlier_interp[j,:], i) * M_sol_to_kg
        relax_time_r[j,i] = kappa*((np.sqrt(G*(M_BH+M_enc)/(10**lograd_interp[i]*pc_to_m)))**3)/(G**2*10**logdens_steep_outlier_interp[j,i]*(M_sol_to_kg/pc_to_m**3)*M_star*np.log(0.4*M_BH/M_star))
        relax_time_r[j,i] /= years_to_sec

plt.figure(dpi=500)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(t$_{relax}$ [yr])')
for i in range(len(relax_time_r_all[:,0])):
    #plt.figure(dpi=500)
    plt.plot(lograd_interp, np.log10(relax_time_r_all[i,:]), '0.5')
plt.plot(lograd_interp, np.log10(relax_time_r[0,:]), 'b', label='$\gamma$ < -1.75')
for i in range(len(relax_time_r[:,0])):
    #plt.figure(dpi=500)
    plt.plot(lograd_interp, np.log10(relax_time_r[i,:]), 'b')
plt.legend()
plt.xlim(-1,2)


#%%


# let's plot the relaxation time as f(r) for all galaxies.
# let's plot the relaxation time as f(r) for all galaxies.
relax_time_r_all_106 = np.zeros((len(lograd_sort[:,0]),len(lograd_interp)))
for j in range(len(relax_time_r_all_106[:,0])):
    for i in range(len(relax_time_r_all_106[0,:])):
        M_BH = 10**6 * M_sol_to_kg
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_interp[j,:], i) * M_sol_to_kg
        relax_time_r_all_106[j,i] = kappa*((np.sqrt(G*(M_BH+M_enc)/(10**lograd_interp[i]*pc_to_m)))**3)/(G**2*10**logdens_interp[j,i]*(M_sol_to_kg/pc_to_m**3)*M_star*np.log(0.4*M_BH/M_star))
        relax_time_r_all_106[j,i] /= years_to_sec


# let's plot the relaxation time as f(r) for the outlier galaxies.
relax_time_r_106 = np.zeros((len(lograd_steep_outlier[:,0]),len(lograd_interp)))
for j in range(len(relax_time_r_106[:,0])):
    for i in range(len(relax_time_r_106[0,:])):
        M_BH = 10**6 * M_sol_to_kg
        M_enc = get_enc_mass(10**lograd_interp, 10**logdens_steep_outlier_interp[j,:], i) * M_sol_to_kg
        relax_time_r_106[j,i] = kappa*((np.sqrt(G*(M_BH+M_enc)/(10**lograd_interp[i]*pc_to_m)))**3)/(G**2*10**logdens_steep_outlier_interp[j,i]*(M_sol_to_kg/pc_to_m**3)*M_star*np.log(0.4*M_BH/M_star))
        relax_time_r_106[j,i] /= years_to_sec
        
plt.figure(dpi=500)
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(t$_{relax}$ [yr])')
for i in range(len(relax_time_r_all_106[:,0])):
    #plt.figure(dpi=500)
    plt.plot(lograd_interp, np.log10(relax_time_r_all_106[i,:]), '0.5')
plt.plot(lograd_interp, np.log10(relax_time_r_106[0,:]), 'b', label='$\gamma$ < -1.75')
for i in range(len(relax_time_r_106[:,0])):
    #plt.figure(dpi=500)
    plt.plot(lograd_interp, np.log10(relax_time_r_106[i,:]), 'b')
plt.legend()
plt.xlim(-1,2)
#plt.ylim(6,10.13)




#%%

# plot bh mass vs rho with outliers highlighted
plt.figure(dpi=500)
plt.plot(log_bh_masses_nonans,log_rhos,linestyle='',marker='.')
outlier_inds = high_inds
#= u.get_n_max(r_relaxes_nonans,n_outliers)[1]
plt.plot(log_bh_masses_nonans[outlier_inds],log_rhos[outlier_inds],linestyle='',marker='.')
#plt.plot([np.min(log_bh_masses_nonans),np.max(log_bh_masses_nonans)],[5,5],linestyle='--',color='r')
plt.xlabel('log(M$_{BH}$ [M$_\odot$])')
#plt.ylabel('r$_{relax}$ [pc]')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/$pc^3$])')
#plt.ylim(0,150)
plt.savefig('../Plots/r_relax_investigation/log_rho_v_BH_mass_ET_{}_LT_{}.png'.format(et_age,lt_age),
                 bbox_inches='tight', pad_inches=0.1, dpi=500)

# plot bh mass vs gamma with outliers highlighted
plt.figure(dpi=500)
plt.plot(log_bh_masses_nonans,gammas,linestyle='',marker='.')
outlier_inds = high_inds
#= u.get_n_max(r_relaxes_nonans,n_outliers)[1]
plt.plot(log_bh_masses_nonans[outlier_inds],gammas[outlier_inds],linestyle='',marker='.')
plt.plot([np.min(log_bh_masses),np.max(log_bh_masses)],[-2.25,-2.25],linestyle='--',color='r')
#plt.plot([np.min(log_bh_masses_nonans),np.max(log_bh_masses_nonans)],[5,5],linestyle='--',color='r')
plt.xlabel('log(M$_{BH}$ [M$_\odot$])')
#plt.ylabel('r$_{relax}$ [pc]')
plt.ylabel('$\gamma$')
#plt.ylim(0,150)
plt.savefig('../Plots/r_relax_investigation/gammas_v_BH_mass_ET_{}_LT_{}.png'.format(et_age,lt_age),
                 bbox_inches='tight', pad_inches=0.1, dpi=500)

# plot gamma vs rho with outliers highlighted
plt.figure(dpi=500)
plt.plot(gammas,log_rhos,linestyle='',marker='.')
outlier_inds = high_inds
#= u.get_n_max(r_relaxes_nonans,n_outliers)[1]
plt.plot(gammas[outlier_inds],log_rhos[outlier_inds],linestyle='',marker='.')
plt.plot([-2.25,-2.25], [np.min(log_rhos),np.max(log_rhos)],linestyle='--',color='r')
#plt.plot([np.min(log_bh_masses_nonans),np.max(log_bh_masses_nonans)],[5,5],linestyle='--',color='r')
plt.xlabel('$\gamma$')
#plt.ylabel('r$_{relax}$ [pc]')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/$pc^3$])')
#plt.ylim(0,150)
plt.savefig('../Plots/r_relax_investigation/log_rho_v_gammas_ET_{}_LT_{}.png'.format(et_age,lt_age),
                 bbox_inches='tight', pad_inches=0.1, dpi=500)

# =============================================================================
# r_relax_sl = np.delete(r_relaxes_nonans, np.where(gammas_nonans > -2.25))
# bh_masses_sl = np.delete(log_bh_masses_nonans, np.where(gammas_nonans > -2.25))
# gammas_sl = np.delete(gammas_nonans, np.where(gammas_nonans > -2.25))
# rhos_sl = np.delete(log_rhos_nonans, np.where(gammas_nonans > -2.25))
# 
# 
# # calculate the percentage of r_relax above 5pc
# perc_above = len(np.where(r_relax_sl>5)[0])/len(r_relax_sl)
# 
# plt.figure(dpi=500)
# plt.title('Early = {} Gyr, Late = {} Gyr, $\gamma$ < -2.25'.format(et_age,lt_age))
# a = plt.hist(r_relax_sl, bins=50, range=(0,10))
# plt.xlabel('r$_{relax}$ [pc]')
# plt.xlim(0,10)
# plt.plot([5,5],[0,np.max(a[0])],linestyle='--',color='r')
# plt.text(0.65, 0.75, '{:.2f}% > 5pc'.format(perc_above*100), transform=plt.gcf().transFigure, size=10)
# plt.savefig('../Plots/r_relax_investigation/r_relax_distribution_ET_{}_LT_{}_slopes_above_225.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# n_outliers = 30
# 
# plt.figure(dpi=500)
# plt.plot(bh_masses_sl,r_relax_sl,linestyle='',marker='.')
# outlier_inds = u.get_n_max(r_relax_sl,n_outliers)[1]
# plt.plot(bh_masses_sl[outlier_inds],r_relax_sl[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(bh_masses_sl),np.max(bh_masses_sl)],[5,5],linestyle='--',color='r')
# plt.xlabel('log(M$_{BH}$ [M$_\odot$])')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_BH_mass_ET_{}_LT_{}_slopes_above_225.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# plt.plot(rhos_sl,r_relax_sl,linestyle='',marker='.')
# plt.plot(rhos_sl[outlier_inds],r_relax_sl[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(rhos_sl),np.max(rhos_sl)],[5,5],linestyle='--',color='r')
# plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/$pc^3$])')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_rho_ET_{}_LT_{}_slopes_above_225.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# plt.figure(dpi=500)
# plt.plot(gammas_sl,r_relax_sl,linestyle='',marker='.')
# plt.plot(gammas_sl[outlier_inds],r_relax_sl[outlier_inds],linestyle='',marker='.')
# plt.plot([np.min(gammas_sl),np.max(gammas_sl)],[5,5],linestyle='--',color='r')
# plt.xlabel('$\gamma$')
# plt.ylabel('r$_{relax}$ [pc]')
# plt.savefig('../Plots/r_relax_investigation/r_relax_v_gamma_ET_{}_LT_{}_slopes_above_225.png'.format(et_age,lt_age),
#                  bbox_inches='tight', pad_inches=0.1, dpi=500)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #%%
# 
# 
# 
# # little code for Ben
# mass_bounds = [10,11]
# num = 0
# MW_mass = 0
# other_mass = 0
# for i in range(len(gal_masses)):
#     if gal_masses[i] >= mass_bounds[0] and gal_masses[i] <= mass_bounds[1]:
#         MW_mass += 10**gal_masses[i]
#     else:
#         other_mass += 10**gal_masses[i]
#         
# 
# =============================================================================


