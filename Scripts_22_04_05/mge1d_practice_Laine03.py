#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to perform 1D MGE fits to the SB profiles of galaxies in Lauer2005

Created on Tue Jul 20 13:02:27 2021

@author: christian
"""
import glob
import csv
from mgefit.mge_fit_1d import mge_fit_1d
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import mge1d_util as u
#%%

gal_name = 'a1020'


# modify gal name for getting M/Ls
if gal_name[0] == 'a':
    gal_name_mod = 'A'+gal_name[1:]+'-M1'
    gal_name_nospace = gal_name_mod
elif gal_name[0] == 'n':
    gal_name_mod = 'NGC '+gal_name[1:]
    gal_name_nospace = 'NGC'+gal_name[1:]
elif gal_name[0] == 'v':
    gal_name_mod = 'VCC '+gal_name[1:]
    gal_name_nospace = 'VCC'+gal_name[1:]

# =============================================================================
# number = '3377'
# 
# # check to see if a galaxy is in the Stone & Metzger (2016) paper sample
# if u.check_for_galaxy_in_stone('NGC'+number) == False:
#     print('Galaxy not in Stone&Metzger(2016).')
# =============================================================================


#%%

# =============================================================================
# # specify name of galaxy for fit
# gal_name = 'NGC '+number
# gal_name_nospace = u.format_gal_name(gal_name)
# if len(number) == 3:
#     gal_name_mod = 'NGC 0'+number
# else:
#     gal_name_mod = gal_name
# =============================================================================
# ========================= READ IN DATA ======================================
#%%
SB_files = glob.glob('/Users/christian/OneDrive/Desktop/TDE Code/Laine2003_BCG_data/bcg_profs/*.txt')
all_SB_data = []
abell_nums = []

for i in range(len(SB_files)):
    abell_nums.append(SB_files[i].split('/')[-1].split('.')[0][5:]) 
    file = open(SB_files[i], 'r')
    csvreader = csv.reader(file)
    SBs = []
    rads = []
    for row in csvreader:
        s = row[0].split()
        rads.append(s[0])
        SBs.append(s[1])
    all_SB_data.append([rads,SBs])
abell_nums = np.array(abell_nums)
    #pdb.set_trace()
#%%
# =============================================================================
# # read in fits table with SB data
# SB_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/SBs_Lauer_2005.fit'
# # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
# hdul = fits.open(SB_table_filename)  # open a FITS file
# data1 = hdul[1].data 
# hdul.close()
# pdb.set_trace()
# # read in fits table with general galaxy properties
# gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/gal_properties_Lauer_2005.fit'
# # fields for this file: ['Name','Morph','Dist','r_Dist','VMAG','sigma','N','Filt','PS96','Prof','Simbad','NED','_RA','_DE','recno']
# hdul = fits.open(gal_table_filename)  # open a FITS file
# data2 = hdul[1].data 
# hdul.close()
# =============================================================================
#%%
# read in fits table with general galaxy properties including R_eff from Lauer+2007
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2007_data/gal_names_and_notes_Lauer_2007.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMag','logR','r1','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data3 = hdul[1].data 
hdul.close()
abell_notes = np.core.defchararray.find(data3['Note'],'A')
cross_inds = np.where(abell_notes == 0)
lauer_names = data3['name'][cross_inds]
abell_names = data3['Note'][cross_inds]

cross_abell_nums = []
for i in range(len(abell_names)):
    cross_abell_nums.append(abell_names[i][1:5])
cross_abell_nums = np.array(cross_abell_nums)

# define the galaxy names array using names consistent with Lauer+07
names = []
for i in range(len(abell_nums)):
    ref_ind = np.where(cross_abell_nums == abell_nums[i])
    if len(ref_ind[0]) == 0:
        names.append('A'+abell_nums[i]+'-M1')
    else:
        names.append(lauer_names[ref_ind[0]][0])
names = np.array(names)

#%%
# store the data for specific galaxy
gal_ind = np.where(names == gal_name)[0][0]
    
    
SB_data = np.array(all_SB_data[gal_ind]).astype(np.float)
# store the data for specific galaxy
#SB_data = data1[np.where(data1['name']==gal_name)]
#gal_data = data2[np.where(data2['name']==gal_name)]
gal_data_w_reff = data3[np.where(data3['Name']==gal_name_mod)]

# ensure the galaxy name returned data
#if len(SB_data) == 0:
#    print('*** ERROR: No SB data found for '+gal_name+'. ***')
#if len(gal_data) == 0:
#    print('*** ERROR: No galaxy properties found for '+gal_name+' in the '+
#          'Lauer2005 table. ***')
#if len(gal_data_w_reff) == 0:
#    print('*** ERROR: No galaxy properties found for '+gal_name+' in the '+
#          'Lauer2007 table. ***')
    
# =============================================================================


###############################################################################


# ======================== REPARE DATA FOR FIT ================================

# ensure the radii are logarithmically spaced 
num_rad = 100
rad = np.geomspace(np.min(SB_data[0]), np.max(SB_data[0]), num_rad)
SBs_i = np.interp(rad, SB_data[0], SB_data[1])

#pdb.set_trace()

# can switch interpolaton to cubic spline
#f = interpolate.CubicSpline(SB_data['MajAxis'], SB_data['SuBr'])
#SBs = f(rad)

# plot interpolated data with the original
plt.figure(dpi=500)
#plt.plot(rad,SBs, marker='.')
plt.plot(rad,SBs_i, marker='.', label='Interpolated')
plt.plot(SB_data[0], SB_data[1], marker='+', alpha=0.7, label='Original')
plt.ylim(np.max(SB_data[1])+0.2, np.min(SB_data[1])-0.2)
plt.ylabel('$\mu$ [mag/arcsec$^2$]')
plt.xlabel('Radius [arcsec]')
plt.text(np.max(rad)-7.5, np.min(SBs_i)+0.5, gal_name+' (F555W)')
plt.legend()
#plt.show()
plt.clf()


# convert SBs from mag/arcsec^2 to L_sol/pc^2
#if SB_data['Filt'][0] == 'F555W':
M_sol = 4.84 # Vega system (Lauer95 is in F555W)
# =============================================================================
# elif SB_data['Filt'][0] == 'F814W':
#     M_sol = 4.12 # Vega system
# elif SB_data['Filt'][0] == 'F606W':
#     M_sol = 4.62 # Vega system
# =============================================================================
SBs = 10**((SBs_i - M_sol - 21.572)/(-2.5))

# plot converted SB profile
plt.figure(dpi=500)
plt.plot(rad, SBs)
plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
plt.xlabel('Radius [arcsec]')
plt.text(np.max(rad)-5.5, np.max(SBs)-0.5, gal_name+' (F555W)')
#plt.show()
plt.clf()

# =============================================================================


###############################################################################


# ======================== DEPROJECT SB PROFILE================================

mge = mge_fit_1d(rad, SBs, ngauss=20, plot=True)
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_mge_fit.pdf',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_mge_fit.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# store distance to convert sigma from arcsec to pc
distance = gal_data_w_reff['Dist'][0] * 10**6 # in pc
# conversion factor between radians and arcsec
arcsec_to_radian = 4.8481e-6



# compute M/L ratio for the galaxy
#M_L = 11 # temporary
L_v = 10**(0.4*(-gal_data_w_reff['VMag'][0] + 4.83)) # L_sol_v
# =============================================================================
# if len(gal_data_w_reff) > 0:
#     logRe = gal_data_w_reff['logR'][0] # log(pc)
# else:
#     logRe = 0.0
# disp = gal_data['sigma'][0] # km/s
# G = 4.301e-3 # pc (km/s)^2 M_sol^-1
# if logRe != 0.0:
#     Re = 10**logRe
#     M_L = (2*disp**2*Re)/(G*L_v)
# else:
#     M_L = 4.9*(L_v/10**10)**0.18
# =============================================================================
M_L = 4.9*(L_v/10**10)**0.18

#pdb.set_trace()
print('M/L = ', M_L)

# convert radii from arcsec to pc
radii = (rad*arcsec_to_radian)*distance # in pc

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

print('L_tot from SB data =  ', L_tot_SB, ' L_sol')
print('L_tot from MGE data = ', L_tot_mge, ' L_sol')
print('Ratio (SB/MGE) = ',L_tot_SB/L_tot_mge)

#%%  
# read data from Stone & Metzger (2016) who deprojected densities using nuker law 
# SB profiles (essentially double power-laws) and assuming spherical symmetry
stone_data =  u.get_stone_data(gal_name_nospace)   
stone_rad = stone_data[1]
stone_dens = stone_data[2]

#%%
# plot density profile
plt.figure(dpi=500)
plt.plot(np.log10(radii), np.log10(densities), label='This Work')
if stone_data != 0:
    plt.plot(np.log10(stone_rad), np.log10(stone_dens), alpha=0.7, label='Stone&Metzger2016')
plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.legend()
plt.xlim(np.min(np.log10(radii))-0.2,np.max(np.log10(radii))+0.2)
plt.ylim(min([np.min(np.log10(densities)),np.min(np.log10(stone_dens))])-0.2,
         max([np.max(np.log10(densities)),np.max(np.log10(stone_dens))])-1)
x0,xmax = plt.xlim()
y0,ymax = plt.ylim()
width = np.abs(xmax-x0)
height = np.abs(ymax-y0)
plt.text(x0+width*.05, y0+height*0.05, gal_name+', M$_v$:'+str(gal_data_w_reff['VMAG'][0])+
         ', '+SB_data['Filt'][0]) 

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_density_profile.pdf',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_density_profile.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

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

plt.figure(dpi=500)
plt.plot(rad, SBs, color='b', linestyle='', marker='o')
for i in range(len(mge.sol[0])):
    plt.plot(rad,gausses[i])

plt.plot(rad, summed_gauss, color='orange',)#, alpha=0.7)    

plt.yscale('log')
plt.xscale('log')
plt.ylim(min(SBs)-20,max(SBs)+1000)
plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
plt.xlabel('Radius ["]')


plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_reconstructed_MGE_plot.pdf',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_reconstructed_MGE_plot.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

L_tot_gauss = u.integrate_SB(radii, summed_gauss)

print()
print('L_tot from SB data =     ', L_tot_SB, ' L_sol')
print('L_tot from summed MGE  = ', L_tot_gauss, ' L_sol')


    
