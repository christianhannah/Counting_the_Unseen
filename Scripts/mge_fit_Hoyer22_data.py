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
from astropy.modeling.functional_models import Sersic1D

import astropy.units as uni
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery

import warnings
warnings.filterwarnings("ignore")

#%%

def get_density_profile(gal_name):

    all_ML_types = []
    all_MLs = []
#%%

    # ========================= READ IN DATA ======================================
    
    # Read in new Hoyer+22 table with I_eff data
    all_names_n = []
    all_filts_n = []
    all_i_eff = []
    
    new_table_filename = '../Data_Sets/Hoyer22_data/nsc_parameters.ascii.csv'
    # fields: name	instrument	channel	filter	pa	l_pa	u_pa	ell	l_ell	u_ell	n	l_n	u_n	reff	l_reff	
    #         u_reff	ieff	l_ieff	u_ieff	m	l_m	u_m	v	u_v	i	u_i	logm	u_logm
    with open(new_table_filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] != 'name': # marks the start of data row
                all_names_n.append(row[0])
                all_filts_n.append(row[3])
                all_i_eff.append(row[16])
    
    # convert lists to numpy arrays
    all_names_n = np.array(all_names_n)
    all_filts_n = np.array(all_filts_n)
    all_i_eff = np.array(all_i_eff).astype(float)
    
    #------------------------------------------------------------------------#
    
    
    # Read in Hoyer+22 table b1 which has the galaxy coordinates, pixel scales, 
    # and masses
    all_names_b1 = []
    all_filts_b1 = []
    all_ras = []
    all_decs = []
    all_logmasses_b1 = []
    all_pixel_scales = []
    
    table_b1_filename = '../Data_Sets/Hoyer22_data/table_b1.ascii.csv'
    # fields: ['Name', 'RA', 'DE', 'dm', 'e_dm', 'logm', 'e_logm', 'Instrument', 'Channel', 'Filter', 'PropID', 'Vegamag', 'texp', 'ps']
    with open(table_b1_filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] != 'Name': # marks the start of data row
                all_names_b1.append(row[0])
                all_ras.append(row[1])
                all_decs.append(row[2])
                all_logmasses_b1.append(row[5])
                all_pixel_scales.append(row[13])
                all_filts_b1.append(row[9])
    
    # convert lists to numpy arrays
    all_names_b1 = np.array(all_names_b1)
    all_filts_b1 = np.array(all_filts_b1)
    all_ras = np.array(all_ras).astype(float)
    all_decs = np.array(all_decs).astype(float)
    all_logmasses_b1 = np.array(all_logmasses_b1).astype(float) 
    all_pixel_scales = np.array(all_pixel_scales).astype(float)

    #------------------------------------------------------------------------#

    all_names = []
    all_sersic_inds = []
    all_filts = []
    all_notes = []
    all_r_effs = []
    all_M_Ls = []
    all_vmags = []
    
    # Read in Hoyer+22 table b2 data for the Sersic parameters and M/L's
    table_b2_filename = '../Data_Sets/Hoyer22_data/tab_publication.ascii.csv'
    # fields: ['galaxy', 'notes', 'filter', 'pa', 'l_pa', 'u_pa', 'ell', 'l_ell', 'u_ell', 'n', 'l_n', 'u_n', 'reff', 'u_reff', 'm', 'l_m', 'u_m', 'v', 'u_v', 'i', 'u_i', 'mlr', 'u_mlr', 'logm', 'u_logm']
    with open(table_b2_filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0][0] != '#' and row[0][0] != 'g': # marks the start of data row
                all_names.append(row[0])
                if row[9] != '':
                    all_sersic_inds.append(row[9])
                else:
                    all_sersic_inds.append('-999')
                all_filts.append(row[2])
                all_notes.append(row[1])
                if row[12] != '':
                    all_r_effs.append(row[12])
                else:
                    all_r_effs.append('-999')
                if row[21] != '':
                    all_M_Ls.append(row[21])
                else:
                    all_M_Ls.append('-999')
                if row[17] != '':
                    all_vmags.append(row[17])
                else:
                    all_vmags.append('-999')
    
    
    # convert lists to numpy arrays
    all_names = np.array(all_names)
    all_sersic_inds = np.array(all_sersic_inds).astype(float)
    all_filts = np.array(all_filts)
    all_notes = np.array(all_notes)
    all_r_effs = np.array(all_r_effs).astype(float)
    all_M_Ls = np.array(all_M_Ls).astype(float)
    all_vmags = np.array(all_vmags).astype(float)
            

    for i in range(len(all_names)):
        print(all_names[i], all_notes[i], all_filts[i], all_sersic_inds[i])

    # filter the data for a specific galaxy
    # new table data
    gal_inds_n = np.where(all_names_n == gal_name)[0]
    filts_n = all_filts_n[gal_inds_n]
    i_effs = all_i_eff[gal_inds_n]
    # Table b1 data
    gal_inds_b1 = np.where(all_names_b1 == gal_name)[0]
    RA = np.unique(all_ras[gal_inds_b1])
    DEC = np.unique(all_decs[gal_inds_b1])
    filts_b1 = all_filts_b1[gal_inds_b1]
    logmass = np.unique(all_logmasses_b1[gal_inds_b1])
    pixel_scales = all_pixel_scales[gal_inds_b1]  
    vmags = all_vmags[gal_inds_b1]
    # Table b2 data
    gal_inds = np.where(all_names == gal_name)[0]
    sersic_inds = all_sersic_inds[gal_inds]
    r_effs = all_r_effs[gal_inds]
    M_L = all_M_Ls[gal_inds[0]]
    ML_type = 4
    notes = all_notes[gal_inds]
    filts = all_filts[gal_inds]
    
    
    # ensure the galaxy name returned data
    if len(gal_inds_b1) == 0:
        print('*** ERROR: No RA/Dec found for '+gal_name+'. ***')
        return
    if len(gal_inds) == 0:
        print('*** ERROR: No data found for '+gal_name+'. ***')
        return


    # isolate data for the F814W filter
    ind_814_n = np.where(filts_n == 'f814w')[0]
    i_eff = i_effs[ind_814_n]
    
    ind_814_b1 = np.where(filts_b1 == 'f814w')[0]
    pixel_scale = pixel_scales[ind_814_b1]
    vmag = vmags[ind_814_b1]
    
    ind_814 = np.where(filts == 'f814w')[0]
    sersic_ind = sersic_inds[ind_814]
    r_eff = r_effs[ind_814]
    
    # get F814W zero-point magnitude and exposure time for the galaxy 
    m_zp, t_exp = u.get_zp_and_texp(gal_name)
    
    # convert i_eff from counts/pix^2 to mag/arcsec^2
    i_eff_cpa = i_eff/(pixel_scale**2) # counts/arcsec^2
    m_inst = -2.5*np.log10(i_eff_cpa/t_exp)
    m_cal = m_inst + m_zp
    I_eff_mpa = m_cal # mag/arcsec^2
    
    # correct the surface brightness profiles for extinction
    coords = SkyCoord(ra=RA*uni.degree, dec=DEC*uni.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass F814W
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*1.526 # ACS WFC
    I_eff_mpa_corrected = I_eff_mpa - A_v
    
    
    # convert SBs from mag/arcsec^2 to L_sol/pc^2
    M_sol = 4.12 # Vega system for F814W
    I_eff = 10**((I_eff_mpa_corrected - M_sol - 21.572)/(-2.5))
 
# =============================================================================


###############################################################################


# ======================== REPARE DATA FOR FIT ================================
    #pdb.set_trace()
    # get the distance and other relevant data from David's table
    david_objname, dist, logmass_d, vmag, gal_gi_color, galtype = u.get_David_gal_data(gal_name, RA, DEC)

    # convert sdistance from Mpc to pc
    distance = dist * 10**6 # in pc
    # conversion factor between radians and arcsec
    arcsec_to_radian = 4.8481e-6

    #pdb.set_trace()
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
    rad = np.geomspace(pixel_scale[0]/2., 10, num_rad)
    radii = (rad*arcsec_to_radian)*distance # in pc  
    
    S1 = Sersic1D(amplitude=I_eff, r_eff=r_eff, n=sersic_ind)
    SBs = S1(radii).flatten()
    
# =============================================================================
#     # convert SBs from mag/arcsec^2 to L_sol/pc^2
#     M_sol = 4.12 # Vega system for F814W
#     SBs = 10**((SBs_i - M_sol - 21.572)/(-2.5))
# =============================================================================


# =============================================================================
#     # plot converted SB profile
#     plt.figure(dpi=500)
#     plt.plot(radii, SBs)
#     plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
#     plt.xlabel('Radius [pc]')
#     plt.text(np.max(rad)-5.5, np.max(SBs)-0.5, gal_name+' (F814W)')
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.show()
#     #plt.clf()
#     #pdb.set_trace()
# =============================================================================
# =============================================================================


###############################################################################


# ======================== DEPROJECT SB PROFILE================================

    # perform MGE fit
    mge = mge_fit_1d(rad, SBs, ngauss=20, inner_slope=4, outer_slope=1, plot=False)
    
    # define array to store radial densities
    densities = np.zeros_like(rad)
    
    heights = []
    sigs = []

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
                heights.append(height)
                sigs.append(sig_i)

            # convert total luminosity to total mass using M/L ratio
            M_tot = L_tot*M_L
            # compute desity contribution from one gaussian component
            dens_comp = (M_tot/((2*np.pi)**(3/2)*sig**3))*np.exp(-(radii[j]**2/(2*sig**2))) # M_sol/pc^3
            # add density from each gaussian to total density
            #print(dens_comp, radii[j])
            densities[j] += dens_comp
        #pdb.set_trace()

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
#     plt.savefig('../Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_density_profile.pdf',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#     plt.savefig('../Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_density_profile.png',
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
#     plt.savefig('../Plots/MGE_FITS_and_DENSITIES_pdfs/reconstructed_MGEs/'+gal_name_nospace+'_reconstructed_MGE_plot.pdf',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#     plt.savefig('../Plots/MGE_FITS_and_DENSITIES_pngs/reconstructed_MGEs/'+gal_name_nospace+'_reconstructed_MGE_plot.png',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================

    #pdb.set_trace()  
    return densities, radii, vmag, M_L, ML_type, distance, RA, DEC, 'F814W', pixel_scale, \
        heights, sigs, SBs, rad, np.array(all_MLs), np.array(all_ML_types), gicolor, summed_gauss
    

# =============================================================================
# # bit of code to quickly compare galaxy masses from Nils and David
# names =  ['bts76', 'ddo084', 'kk2000-03', 'kk2000-53', 'kk96', 'leg09',
#            'lvj1217+4703', 'ngc5011c', 'pgc4310323', 'ugc07242', 'eso553-046']
# 
# a = get_density_profile(names[1])
# 
# #%%
# logm_d = np.array([8.52871022,7.0731452,6.99314961,7.51864636,6.6119821,7.62864264])
# logm_n = np.array([7.65,6.68,7.08,7.24,6.55,7.68])
# 
# plt.figure(dpi=500)
# plt.plot(logm_d, logm_n, linestyle='', marker='.')
# plt.plot([6.5,8.7],[6.5,8.7],linestyle='--',color='r')
# plt.ylabel('log(M$_*$ [M$_\odot$]) from Hoyer+22')
# plt.xlabel('log(M$_*$ [M$_\odot$]) from David')
# =============================================================================


