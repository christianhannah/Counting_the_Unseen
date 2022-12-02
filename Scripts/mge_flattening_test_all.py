#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 23:08:39 2022

@author: christian
"""

from tqdm import tqdm
import os
from matplotlib.collections import LineCollection


#%%

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:08:34 2021

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
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

#%%

def get_density_profile(desc,number,rad_ind):

    # check to see if a galaxy is in the Stone & Metzger (2016) paper sample
    if u.check_for_galaxy_in_stone(desc+number) == False:
        print('Galaxy not in Stone&Metzger(2016).')


#%%

    # specify name of galaxy for fit
    gal_name = desc+' '+number
    gal_name_nospace = u.format_gal_name(gal_name)
    if len(number) == 3:
        gal_name_mod = 'NGC 0'+number # for Lauer 2007 Reff table
    else:
        gal_name_mod = gal_name

    # ========================= READ IN DATA ======================================

    # read in fits table with SB data
    SB_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/SBs_Lauer_2005.fit'
    # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
    hdul = fits.open(SB_table_filename)  # open a FITS file
    data1 = hdul[1].data 
    hdul.close()

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
        return
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
    min_rad = np.min(SB_data['MajAxis'][rad_ind:])
    rad = np.geomspace(np.min(SB_data['MajAxis']), np.max(SB_data['MajAxis']), num_rad)
    SBs_i = np.interp(rad, SB_data['MajAxis'], SB_data['SuBr'])
    
    # truncate the rad and SBs_i to remove innermost points
    rad = rad[rad_ind:]
    SBs_i = SBs_i[rad_ind:]
    # can switch interpolaton to cubic spline
    #f = interpolate.CubicSpline(SB_data['MajAxis'], SB_data['SuBr'])
    #SBs = f(rad)

# =============================================================================
#     # plot interpolated data with the original
#     plt.figure(dpi=500)
#     #plt.plot(rad,SBs, marker='.')
#     plt.plot(rad,SBs_i, marker='.', label='Interpolated')
#     plt.plot(SB_data['MajAxis'], SB_data['SuBr'], marker='+', alpha=0.7, label='Original')
#     plt.ylim(np.max(SB_data['SuBr'])+0.2, np.min(SB_data['SuBr'])-0.2)
#     plt.ylabel('$\mu$ [mag/arcsec$^2$]')
#     plt.xlabel('Radius [arcsec]')
#     plt.text(np.max(rad)-7.5, np.min(SBs_i)+0.5, gal_name+' ('+str(SB_data['Filt'][0])+')')
#     plt.legend()
#     #plt.show()
#     plt.clf()
# =============================================================================
    
    
    # convert SBs from mag/arcsec^2 to L_sol/pc^2
    if SB_data['Filt'][0] == 'F555W':
        M_sol = 4.84 # Vega system
    elif SB_data['Filt'][0] == 'F814W':
        M_sol = 4.12 # Vega system
    elif SB_data['Filt'][0] == 'F606W':
        M_sol = 4.62 # Vega system
    SBs = 10**((SBs_i - M_sol - 21.572)/(-2.5))

# =============================================================================
#     # plot converted SB profile
#     plt.figure(dpi=500)
#     plt.plot(rad, SBs)
#     plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
#     plt.xlabel('Radius [arcsec]')
#     plt.text(np.max(rad)-5.5, np.max(SBs)-0.5, gal_name+' ('+str(SB_data['Filt'][0])+')')
#     #plt.show()
#     plt.clf()
# =============================================================================

# =============================================================================


###############################################################################


# ======================== DEPROJECT SB PROFILE================================

    mge = mge_fit_1d(rad, SBs, ngauss=20, plot=True)
# =============================================================================
#     plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/'+gal_name_nospace+'_mge_fit.pdf',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#     plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/'+gal_name_nospace+'_mge_fit.png',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================
    plt.clf()

    # store distance to convert sigma from arcsec to pc
    distance = gal_data['Dist'][0] * 10**6 # in pc
    # conversion factor between radians and arcsec
    arcsec_to_radian = 4.8481e-6


    
    # compute M/L ratio for the galaxy
    #M_L = 11 # temporary
    vmag = gal_data['VMAG'][0]
    L_v = 10**(0.4*(-vmag + 4.83)) # L_sol_v
    if len(gal_data_w_reff) > 0:
        logRe = gal_data_w_reff['logR'][0] # log(pc)
    else:
        logRe = 0.0
    disp = gal_data['sigma'][0] # km/s
    G = 4.301e-3 # pc (km/s)^2 M_sol^-1
    if logRe != 0.0:
        Re = 10**logRe
        M_L = (2*disp**2*Re)/(G*L_v)
        ML_type = 0
    else:
        M_L = 4.9*(L_v/10**10)**0.18
        ML_type = 1
    M_L_1 = 4.9*(L_v/10**10)**0.18
    print('M/L = ', M_L)
    
    # convert radii from arcsec to pc
    radii = (rad*arcsec_to_radian)*distance # in pc
    min_rad = (min_rad*arcsec_to_radian)*distance # in pc
    
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
    if stone_data != 0:
        stone_rad = stone_data[1]
        stone_dens = stone_data[2]
    else:
        stone_rad = 0
        stone_dens = 0

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
# =============================================================================


# =============================================================================
#     plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pdfs/reconstructed_MGEs/'+gal_name_nospace+'_reconstructed_MGE_plot.pdf',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
#     plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/MGE_FITS_and_DENSITIES_pngs/reconstructed_MGEs/'+gal_name_nospace+'_reconstructed_MGE_plot.png',
#                 bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================

    L_tot_gauss = u.integrate_SB(radii, summed_gauss)

    print()
    print('L_tot from SB data =     ', L_tot_SB, ' L_sol')
    print('L_tot from summed MGE  = ', L_tot_gauss, ' L_sol')
    
    return densities, radii, vmag, M_L, ML_type, distance, min_rad




#%%
# =============================================================================
# ======================= Lauer et al. (2005) =================================
# =============================================================================

# read in fits table with general galaxy properties
gal_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/gal_properties_Lauer_2005.fit'
# fields for this file: ['Name','Morph','Dist','r_Dist','VMAG','sigma','N','Filt','PS96','Prof','Simbad','NED','_RA','_DE','recno']
hdul = fits.open(gal_table_filename)  # open a FITS file
data = hdul[1].data 
hdul.close()
#%%
old_stdout = sys.stdout
#%%
names = data['name']
distances = data['Dist'] * 10**6 # in pc

# remove the galaxies whose centers were too obscurred for analysis in Lauer05
names = names[np.where(names != 'NGC 2768')] # no data for this galaxy
names = names[np.where(names != 'NGC 3557')] # no data for this galaxy
names = names[np.where(names != 'NGC 4125')] # no data for this galaxy
names = names[np.where(names != 'NGC 4786')] # no data for this galaxy
names = names[np.where(names != 'NGC 4936')] # no data for this galaxy
names = names[np.where(names != 'NGC 5018')] # no data for this galaxy
names = names[np.where(names != 'NGC 5322')] # no data for this galaxy
names = names[np.where(names != 'NGC 5845')] # no data for this galaxy
names = names[np.where(names != 'NGC 6776')] # no data for this galaxy
names = names[np.where(names != 'NGC 7626')] # no data for this galaxy
names = names[np.where(names != 'IC 3370')] # no data for this galaxy
names = names[np.where(names != 'IC 4296')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 2768')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 3557')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 4125')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 4786')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 4936')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 5018')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 5322')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 5845')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 6776')] # no data for this galaxy
distances = distances[np.where(names != 'NGC 7626')] # no data for this galaxy
distances = distances[np.where(names != 'IC 3370')] # no data for this galaxy
distances = distances[np.where(names != 'IC 4296')] # no data for this galaxy

# read in fits table with SB data
SB_table_filename = '/Users/christian/OneDrive/Desktop/TDE Code/Lauer2005_data/SBs_Lauer_2005.fit'
# fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
hdul = fits.open(SB_table_filename)  # open a FITS file
data1 = hdul[1].data 
hdul.close()

# conversion factor between radians and arcsec
arcsec_to_radian = 4.8481e-6

names_after_5pc_cut = []
for i in range(len(names)):
    SB_dat = data1[np.where(data1['name']==names[i])]
    r = (SB_dat['MajAxis']*arcsec_to_radian)*distances[i] 
    if r[0] <= 5:
        names_after_5pc_cut.append(names[i])
names = np.array(names_after_5pc_cut)
#%%

min_rad_inds = np.array([4,3,2,1,0])
all_slopes = np.zeros((len(names),len(min_rad_inds)))
all_rads = np.zeros((len(names),len(min_rad_inds)))

for w in tqdm(range(len(names)), position=0, leave=True):

    dens_profiles = np.zeros((100,len(min_rad_inds)))
    radii = np.zeros((100,len(min_rad_inds)))
    lums = np.zeros(len(min_rad_inds))
    ML_types = np.zeros(len(min_rad_inds))
    MLs = np.zeros(len(min_rad_inds))
    dists = np.zeros(len(min_rad_inds))
    slopes = np.zeros(len(min_rad_inds))
    cen_dens = np.zeros(len(min_rad_inds))

    # =============================================================================
    # min_rads = np.zeros(len(min_rad_inds))
    # for i in range(len(min_rad_inds)):
    #     min_rads[i] = np.min(radii[min_rad_inds[i]:,i])
    # =============================================================================

    for i in range(len(min_rad_inds)):
        s = names[w].split(' ')
        dens_profiles[min_rad_inds[i]:,i], radii[min_rad_inds[i]:,i], lums[i], MLs[i], ML_types[i], dists[i], min_rad = get_density_profile(s[0], s[1], min_rad_inds[i])
        #min_rad = np.min(radii[min_rad_inds[0]:,0])
        #min_ind = u.find_nearest(radii[min_rad_inds[i]:,i], min_rad)
        min_ind = 4
        #slopes[i], cen_dens[i] = u.get_density_and_slope(np.log10(dens_profiles[min_rad_inds[i]+min_ind:min_rad_inds[i]+min_ind+3,i]), 
        #                                                 np.log10(radii[min_rad_inds[i]+min_ind:min_rad_inds[i]+min_ind+3,i]))
        slopes[i], cen_dens[i] = u.get_density_and_slope(np.log10(radii[4:7,i]),np.log10(dens_profiles[4:7,i]))
        #slopes[i], cen_dens[i] = u.get_density_and_slope(np.log10(radii[min_rad_inds[i]:,i]),np.log10(dens_profiles[min_rad_inds[i]:,i]))
        #pdb.set_trace()
    
    all_slopes[w,:] = slopes
    all_rads[w,:] = np.abs(radii[0,4] - radii[min_rad_inds,np.flip(min_rad_inds)])
    # =============================================================================
    # =============================================================================
    # =============================================================================
#%%    

    ls = []
    for i in range(len(min_rad_inds)):
        a = np.zeros((len(radii[:,i]),2))
        a[:,0] = np.log10(radii[:,i])
        a[:,1] = np.log10(dens_profiles[:,i])
        ls.append(a)

# =============================================================================
#     fig, ax = plt.subplots(dpi=500)
#     lines = LineCollection(ls, array=min_rad_inds, cmap='viridis', alpha=0.6)
#     ax.add_collection(lines)
#     fig.colorbar(lines, label='# of removed points')
#     ax.set_xlabel('log(Radius [pc])')
#     ax.set_ylabel('log(Density [M$_{\odot}$/pc$^3$])')
#     ax.autoscale()
#     #ax.set_xlim(0.1,0.25)
#     #ax.set_ylim(5.2,5.35)
#     #plt.text(0.83,2.255, names[name_ind]+',  M$_v$='+str(lums[0])) 
#     plt.show()
#     plt.close()
# =============================================================================


#%%
    min_rads = np.zeros(len(min_rad_inds))
    for i in range(len(min_rad_inds)):
        min_rads[i] = np.min(radii[min_rad_inds[i]:,i])

# =============================================================================
#     plt.figure(dpi=500)
#     plt.plot(min_rad_inds, slopes, marker='.', linestyle='', color='b')
#     plt.xlabel('# of removed points')#'Minimum Radius[pc]')
#     plt.ylabel('Power-law Slope')
#     #plt.text(8.3,-4.5, names[name_ind]+',  M$_v$='+str(lums[0])) 
# =============================================================================



    dens_at_min = np.zeros(len(min_rad_inds))
    for i in range(len(min_rad_inds)):
        min_ind = u.find_nearest(radii[:,i], min_rads[0])
        dens_at_min[i] = dens_profiles[min_ind,i]

    #dens_at_cen = np.zeros(len(min_rad_inds))
    #for i in range(len(min_rad_inds)):
    

# =============================================================================
#     plt.figure(dpi=500)
#     plt.plot(min_rad_inds,np.log10(dens_at_min))
#     plt.xlabel('# of removed points')
#     plt.ylabel('log(3D Density at Least Common Radius [M$_{\odot}$/pc$^3$])')
#     #plt.text(8.3,-4.5, names[name_ind]+',  M$_v$='+str(lums[0])) 
#     
# 
# 
#     plt.figure(dpi=500)
#     plt.plot(min_rad_inds,cen_dens)
#     plt.xlabel('# of removed points')
#     plt.ylabel('log(3D Density at Least Common Radius [M$_{\odot}$/pc$^3$])')
#     #plt.text(8.3,-4.5, names[name_ind]+',  M$_v$='+str(lums[0])) 
# 
# =============================================================================
#%%

# insert code to plot the distributions of and median changes in slopes per 
# removed data point


devs = np.zeros((len(names),len(min_rad_inds)-1))
rad_removed = np.zeros((len(names),len(min_rad_inds)-1))
for i in range(len(names)):
        devs[i,:] = all_slopes[i,-1] - all_slopes[i,:-1]
        rad_removed[i,:] = all_rads[i,:-1]

plt.figure(dpi=500)
plt.plot(rad_removed, devs, linestyle='', marker='.', color='b')
plt.xlabel('Radial Coverage Removed [pc]')
plt.ylabel('Deviation of Power-Law Slope')


plt.figure(dpi=500)
plt.plot(all_rads, all_slopes, linestyle='',marker='.',color='b')
plt.xlabel('Radial Coverage Removed [pc]')
plt.ylabel('Power-Law Slope')


med_devs = np.zeros(len(min_rad_inds)-1)
for i in range(len(med_devs)):
    med_devs[i] = np.median(devs[:,i])
    
    
plt.figure(dpi=500)
plt.plot(min_rad_inds[:-1], med_devs, linestyle='', marker='.', color='b')
plt.xlabel('# of Removed Data Points')
plt.ylabel('Median Deviation of Power-law Slope')

labs = ['4 Removed Data Points','3 Removed Data Points','2 Removed Data Points',
        '1 Removed Data Points']
for i in range(len(devs[0,:])):
    plt.figure(dpi=500)
    plt.hist(devs[:,i], bins=30)
    plt.xlabel('Deviation of Power-law Slope')
    plt.ylabel('Counts')
    plt.title(labs[i])
    
    






