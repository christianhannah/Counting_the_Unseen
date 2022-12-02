#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 09:44:38 2022

@author: christian
"""

from mgefit.mge_fit_1d import mge_fit_1d
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import mge1d_util as u
#%%

number = '4278'

# check to see if a galaxy is in the Stone & Metzger (2016) paper sample
if u.check_for_galaxy_in_stone('NGC'+number) == False:
    print('Galaxy not in Stone&Metzger(2016).')


#%%

# specify name of galaxy for fit
gal_name = 'NGC '+number
gal_name_nospace = u.format_gal_name(gal_name)
if len(number) == 3:
    gal_name_mod = 'NGC 0'+number
else:
    gal_name_mod = gal_name
    
    
# ========================= READ IN DATA ======================================

# read in fits table with SB data
SB_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/SBs_Lauer_2005.fit'
# fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
hdul = fits.open(SB_table_filename)  # open a FITS file
data1 = hdul[1].data 
hdul.close()
#pdb.set_trace()
# read in fits table with general galaxy properties
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/gal_properties_Lauer_2005.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMAG','sigma','N','Filt','PS96','Prof','Simbad','NED','_RA','_DE','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data2 = hdul[1].data 
hdul.close()

# read in fits table with general galaxy properties including R_eff from Lauer+2007
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2007_data/gal_properties_Lauer_2007.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMag','logR','r1','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data3 = hdul[1].data 
hdul.close()


# store the data for specific galaxy
SB_data = data1[np.where(data1['name']==gal_name)]
gal_data = data2[np.where(data2['name']==gal_name)]
gal_data_w_reff = data3[np.where(data3['name']==gal_name_mod)]

# ensure the galaxy name returned data
if len(SB_data) == 0:
    print('*** ERROR: No SB data found for '+gal_name+'. ***')
if len(gal_data) == 0:
    print('*** ERROR: No galaxy properties found for '+gal_name+' in the '+
          'Lauer2005 table. ***')
if len(gal_data_w_reff) == 0:
    print('*** ERROR: No galaxy properties found for '+gal_name+' in the '+
          'Lauer2007 table. ***')
    
# =============================================================================


###############################################################################


# ======================== REPARE DATA FOR FIT ================================

# ensure the radii are logarithmically spaced 
num_rad = 100
rad = np.geomspace(np.min(SB_data['MajAxis']), np.max(SB_data['MajAxis']), num_rad)
SBs_i = 10**np.interp(np.log10(rad), np.log10(SB_data['MajAxis']), np.log10(SB_data['SuBr']))

#%%
# =============================================================================
# # plot interpolated data with the original
# plt.figure(dpi=500)
# #plt.plot(rad,SBs, marker='.')
# plt.plot(rad,SBs_i, label='Interpolated in log')
# plt.plot(rad,SBs_1, label='Interpolated')
# plt.plot(SB_data['MajAxis'], SB_data['SuBr'], marker='+', linestyle='',alpha=0.7, label='Original')
# #plt.ylim(np.max(SB_data['SuBr'])+0.2, np.min(SB_data['SuBr'])-0.2)
# plt.ylabel('$\mu$ [mag/arcsec$^2$]')
# plt.xlabel('Radius [arcsec]')
# #plt.text(np.max(rad)-7.5, np.min(SBs_i)+0.5, gal_name+' ('+str(SB_data['Filt'][0])+')')
# plt.legend()
# 
# plt.xlim(0.02,0.06)
# plt.ylim(12.7,11.7)
# 
# #pdb.set_trace()
# #plt.show()
# plt.clf()
# =============================================================================
#%%

# convert SBs from mag/arcsec^2 to L_sol/pc^2
if SB_data['Filt'][0] == 'F555W':
    M_sol = 4.84 # Vega system
elif SB_data['Filt'][0] == 'F814W':
    M_sol = 4.12 # Vega system
elif SB_data['Filt'][0] == 'F606W':
    M_sol = 4.62 # Vega system
SBs = 10**((SBs_i - M_sol - 21.572)/(-2.5))

# =============================================================================


###############################################################################


# ======================== DEPROJECT SB PROFILE================================

mge = mge_fit_1d(rad, SBs, ngauss=20, plot=True, inner_slope=2, outer_slope=3)

#plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_mge_fit.png',
#            bbox_inches='tight', pad_inches=0.1, dpi=500)


# store distance to convert sigma from arcsec to pc
distance = gal_data['Dist'][0] * 10**6 # in pc
# conversion factor between radians and arcsec
arcsec_to_radian = 4.8481e-6


# compute M/L ratio for the galaxy
L_v = 10**(0.4*(-gal_data['VMAG'][0] + 4.83)) # L_sol_v
if len(gal_data_w_reff) > 0:
    logRe = gal_data_w_reff['logR'][0] # log(pc)
else:
    logRe = 0.0
disp = gal_data['sigma'][0] # km/s
G = 4.301e-3 # pc (km/s)^2 M_sol^-1
if logRe != 0.0:
    Re = 10**logRe
    M_L = (2*disp**2*Re)/(G*L_v)
else:
    M_L = 4.9*(L_v/10**10)**0.18
M_L_1 = 4.9*(L_v/10**10)**0.18

#pdb.set_trace()
print('M/L = ', M_L)

# convert radii from arcsec to pc
radii = (rad*arcsec_to_radian)*distance # in pc

# define HST pixel scale in pc for slope calculation
scale = (0.0456*arcsec_to_radian)*distance

#%%
# compute the inner power-law slope of the 2D SB data
#SB_slope, blah,blah,blah = u.get_density_and_slope(np.log10(radii), np.log10(SBs), scale, 10, 2)
SB_slope, blah = u.get_density_and_slope_simple(np.log10(radii), np.log10(SBs))                            

#%%

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

#print('L_tot from SB data =  ', L_tot_SB, ' L_sol')
#print('L_tot from MGE data = ', L_tot_mge, ' L_sol')
#print('Ratio (SB/MGE) = ',L_tot_SB/L_tot_mge)

#%%
# let's measure the innermost slope of smallest sigma gaussian
A = mge.sol[0][2]/(np.sqrt(2*np.pi)*mge.sol[1][2])
sigma = mge.sol[1][2]
# =============================================================================
# def gauss_deriv(x,A,sigma):
#     return -x/(np.sqrt(2*np.pi)*sigma**3)*np.exp(-x**2/(2*sigma**2))
# =============================================================================
def gauss(height, sig, r):
    return height*np.exp((-0.5)*r**2/sig**2)
y = np.log10(gauss(A,sigma,rad[0:2]))
x = np.log10(rad[0:2])
s = (y[1]-y[0])/(x[1]-x[0])


print('slope @ minimum radius for smallest Gaussian: ',s)

#%%

# get power_law_index and density at 5pc
slope, inner_dens, blah, blah = u.get_density_and_slope(np.log10(radii),np.log10(densities), scale, 10, 2)

#%%  
# read data from Stone & Metzger (2016) who deprojected densities using nuker law 
# SB profiles (essentially double power-laws) and assuming spherical symmetry
stone_data =  u.get_stone_data(gal_name_nospace)   
if stone_data != 0:
    stone_rad = stone_data[1]
    stone_dens = stone_data[2]

#%%
# plot density profile
# =============================================================================
# plt.figure(dpi=500)
# plt.plot(np.log10(radii), np.log10(densities), label='This Work')
# if stone_data != 0:
#     plt.plot(np.log10(stone_rad), np.log10(stone_dens), alpha=0.7, label='Stone&Metzger2016')
# plt.xlabel('log(Radius [pc])')
# plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
# plt.legend()
# plt.xlim(np.min(np.log10(radii))-0.2,np.max(np.log10(radii))+0.2)
# plt.ylim(min([np.min(np.log10(densities)),np.min(np.log10(stone_dens))])-0.2,
#          max([np.max(np.log10(densities)),np.max(np.log10(stone_dens))])-1)
# x0,xmax = plt.xlim()
# y0,ymax = plt.ylim()
# width = np.abs(xmax-x0)
# height = np.abs(ymax-y0)
# plt.text(x0+width*.05, y0+height*0.05, gal_name+', M$_v$:'+str(gal_data['VMAG'][0])+
#          ', '+SB_data['Filt'][0]) 
# 
# =============================================================================
#plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_density_profile.pdf',
#            bbox_inches='tight', pad_inches=0.1, dpi=500)
#plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_density_profile.png',
#            bbox_inches='tight', pad_inches=0.1, dpi=500)

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
# plt.figure(dpi=500)
# plt.plot(rad, SBs, color='b', linestyle='', marker='o')
# for i in range(len(mge.sol[0])):
#     plt.plot(rad,gausses[i])
# 
# plt.plot(rad, summed_gauss, color='orange',)#, alpha=0.7)    
# 
# plt.yscale('log')
# plt.xscale('log')
# plt.ylim(min(SBs)-20,max(SBs)+1000)
# plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
# plt.xlabel('Radius ["]')
# =============================================================================


#plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_reconstructed_MGE_plot.pdf',
#            bbox_inches='tight', pad_inches=0.1, dpi=500)
#plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_reconstructed_MGE_plot.png',
#            bbox_inches='tight', pad_inches=0.1, dpi=500)

L_tot_gauss = u.integrate_SB(radii, summed_gauss)


# compute the inner power-law slope of summed MGE just using 3 innermost points
MGE_slope, blah = u.get_density_and_slope_simple(np.log10(radii), np.log10(summed_gauss))

print('True SB slope: ',SB_slope)
print('MGE SB slope:  ',MGE_slope) 

print()
print('L_tot from SB data =     ', L_tot_SB, ' L_sol')
print('L_tot from summed MGE  = ', L_tot_gauss, ' L_sol')
print()
print('Ratio (SB/MGE) = ',L_tot_SB/L_tot_gauss)


    
