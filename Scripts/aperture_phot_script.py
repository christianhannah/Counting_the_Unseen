#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 00:43:18 2022

@author: christian
"""

# Write code to conduct 0.5" - 1.0" aperture photometry for galaxies in our 
# sample without nuclear color information.

import glob
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
from photutils.aperture.photometry import aperture_photometry
from photutils.aperture.circle import (CircularAperture, SkyCircularAperture)
from astropy.io import fits
import pdb
from astropy.wcs import WCS
from dustmaps.sfd import SFDQuery
from photutils.centroids import centroid_quadratic
from photutils.centroids import centroid_2dg
from photutils.centroids import centroid_com
import matplotlib.pyplot as plt
from scipy import stats
import mge1d_util as ut
from acstools import acszpt
from tqdm import tqdm
import time
import matplotlib.colors as colors
import warnings
warnings.filterwarnings("ignore")

#%%


filenames_wfpc2_f555w = glob.glob('../Data_Sets/COLOR_IMAGES/WFPC2_F555W/*.fits')
filenames_wfpc2_f814w = glob.glob('../Data_Sets/COLOR_IMAGES/WFPC2_F814W/*.fits')
filenames_wfpc2_f606w = glob.glob('../Data_Sets/COLOR_IMAGES/WFPC2_F606W/*.fits')
filenames_wfpc2_f547m = glob.glob('../Data_Sets/COLOR_IMAGES/WFPC2_F547M/*.fits')
filenames_wfc3_f547m = glob.glob('../Data_Sets/COLOR_IMAGES/WFC3_F547M/*.fits')
filenames_wfc3_f475w = glob.glob('../Data_Sets/COLOR_IMAGES/WFC3_F475W/*.fits')
filenames_wfc3_f814w = glob.glob('../Data_Sets/COLOR_IMAGES/WFC3_F814W/*.fits')
filenames_acs_f850lp = glob.glob('../Data_Sets/COLOR_IMAGES/ACS_F850LP/*.fits')
filenames_acs_f475w = glob.glob('../Data_Sets/COLOR_IMAGES/ACS_F475W/*.fits')
filenames_acs_f814w = glob.glob('../Data_Sets/COLOR_IMAGES/ACS_F814W/*.fits')

#%%
print("======================================")
print("==== BEGINNING APERTURE PHOTMETRY ====")
print("======================================")

# =============================================================================
# =============================================================================
# =============================================================================
# =========================       WFPC2 F555W        ==========================
# =============================================================================
# =============================================================================
# ============================================================================= 
#%%

names_wfpc2_f555w = []
RAs = np.zeros(len(filenames_wfpc2_f555w))
DECs = np.zeros(len(filenames_wfpc2_f555w))
fluxes = np.zeros(len(filenames_wfpc2_f555w))
f555w_mags = np.zeros(len(filenames_wfpc2_f555w))
print()
print('WFPC2 F555W')
print()
time.sleep(1)
for i in tqdm(range(len(filenames_wfpc2_f555w)), position=0, leave=True):
    
    hdul = fits.open(filenames_wfpc2_f555w[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    #cutind = int(np.min((data.shape[0]*0.20, data.shape[1]*0.20)))
    cutind = 0
    #data = data[cutind:-cutind,cutind:-cutind]
    
    
    # determine galaxy name and specify center coordinates
    names_wfpc2_f555w.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-200):int(position[0]+200),int(position[1]-200):int(position[1]+200)]

# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    # create a centroid mask
    mask = np.zeros_like(data).astype(bool)
# =============================================================================
#     mask_inds = np.where(data > 10)
#     mask[mask_inds] = False
# =============================================================================
    mask_inds = np.where(data < 0)
    mask[mask_inds] = False

# =============================================================================
#     cen = np.zeros(2)
#     m = np.where(data == np.max(data))
#     cen[0] = m[0][0]
#     cen[1] = m[1][0]
# =============================================================================
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct

    cen = np.zeros(2)
    m = np.where(data[int(sector[1]/n_sec*n_x):int((sector[1]+1)/n_sec*n_x),int(sector[0]/n_sec*n_y):int((sector[0]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[1]/n_sec*n_x):int((sector[1]+1)/n_sec*n_x),int(sector[0]/n_sec*n_y):int((sector[0]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    #pdb.set_trace()
# =============================================================================
#     x_dist = np.zeros(data.shape[0])
#     for v in range(len(x_dist)):
#         x_dist[v] = np.sum(data[v,:])
#     
#     y_dist = np.zeros(data.shape[1])
#     for v in range(len(y_dist)):
#         y_dist[v] = np.sum(data[:,v])
#         
#     cen = np.zeros(2)
#     cen[0] = ut.find_nearest(x_dist,np.median(x_dist))
#     cen[1] = ut.find_nearest(y_dist,np.median(y_dist))
# =============================================================================

    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 3.483e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    zero_pt = -22.545	
    
    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass WFPC2 F555W 
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*2.755
    
     
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(RAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================

    
    
    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    #position = cen_zoom
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[0]-0.2*n_x):int(cen_tot[0]+0.2*n_x),int(cen_tot[1]-0.2*n_y):int(cen_tot[1]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_wfpc2_f555w_image_w_aperture.png'.format(names_wfpc2_f555w[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data_all, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
#     #plt.imshow(data[int(cen_tot[0]-0.2*n_x):int(cen_tot[0]+0.2*n_x),int(cen_tot[1]-0.2*n_y):int(cen_tot[1]+0.2+n_y)],
#     #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
#     #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
#     circ = plt.Circle(cen_tot,aperture.r,alpha=0.3,color='r')
#     plt.gca().add_patch(circ)
#     #plt.plot(137,776,linestyle='',marker='.',color='r')
#     plt.colorbar()
#     #plt.show()
#     
#     plt.savefig('../Plots/aperture_photometry/whole_gal_w_aperture/{}_wfpc2_f555w_image_w_aperture.png'.format(names_wfpc2_f555w[i]),
#                      bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================
    
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['ATODGAIN']
    if gain != 7:
        if gain == 14:
            photflam += 1.987
            zero_pt += -2.5*np.log10(1.987)
    fluxes[i] = (phot_table['aperture_sum']/gain) # erg/s/cm^2/A
    f555w_mags[i] = - 2.5*np.log10(fluxes[i]) + zero_pt - A_v
    
    #pdb.set_trace()
    
#%%
plt.figure(dpi=500)
plt.hist(f555w_mags,bins=10)
plt.xlabel('m$_{F555W}$')
    
#pdb.set_trace()
    
    
# =============================================================================
# =============================================================================
# =============================================================================
# =========================       WFPC2 F814W        ==========================
# =============================================================================
# =============================================================================
# =============================================================================    
    
    
#%%

names_wfpc2_f814w = []
RAs = np.zeros(len(filenames_wfpc2_f814w))
DECs = np.zeros(len(filenames_wfpc2_f814w))
fluxes = np.zeros(len(filenames_wfpc2_f814w))
f814w_mags = np.zeros(len(filenames_wfpc2_f814w))
print()
print()
print('WFPC2 F814W')
time.sleep(1)

for i in tqdm(range(len(filenames_wfpc2_f814w)), position=0, leave=True):
    
    hdul = fits.open(filenames_wfpc2_f814w[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
# =============================================================================
    
    cutind = 0
    #data = data[cutind:-cutind,cutind:-cutind]
    
    # determine galaxy name and specify center coordinates
    names_wfpc2_f814w.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-200):int(position[0]+200),int(position[1]-200):int(position[1]+200)]

    
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    # create a centroid mask
    mask = np.zeros_like(data).astype(bool)
    mask_inds = np.where(data > 50)
    mask[mask_inds] = False
    mask_inds = np.where(data < -1)
    mask[mask_inds] = False

    #cen = np.where(data == np.max(data))
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(np.log10(data.flatten()))
    med = np.median(np.log10(data.flatten()))
    sigfac = 5
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(np.log10(data_zoom))
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 2.508e-18	
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    zero_pt = -21.639	
    

    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass WFPC2 F814W 
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*1.549
    
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================


    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_wfpc2_f814w_image_w_aperture.png'.format(names_wfpc2_f814w[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['ATODGAIN']
    if gain != 7:
        if gain == 14:
            photflam += 1.987
            zero_pt += -2.5*np.log10(1.987)
    fluxes[i] = (phot_table['aperture_sum']/gain) # erg/s/cm^2/A
    f814w_mags[i] = - 2.5*np.log10(fluxes[i]) + zero_pt - A_v
    
   # pdb.set_trace()
    
    
    
    
# =============================================================================
# =============================================================================
# =============================================================================
# ========================       WFPC2 F606W        ===========================
# =============================================================================
# =============================================================================
# ============================================================================= 

#%%

names_wfpc2_f606w = []
RAs = np.zeros(len(filenames_wfpc2_f606w))
DECs = np.zeros(len(filenames_wfpc2_f606w))
fluxes = np.zeros(len(filenames_wfpc2_f606w))
f606w_mags = np.zeros(len(filenames_wfpc2_f606w))
print()
print()
print('WFPC2 F606W')
time.sleep(1)

for i in tqdm(range(len(filenames_wfpc2_f606w)), position=0, leave=True):
    
    hdul = fits.open(filenames_wfpc2_f606w[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    cutind = 0
    #data = data[cutind:-cutind,cutind:-cutind]
    
    # determine galaxy name and specify center coordinates
    names_wfpc2_f606w.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-200):int(position[0]+200),int(position[1]-200):int(position[1]+200)]

    
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    #cen = np.where(data == np.max(data))
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 1.900e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    zero_pt = -22.887		

    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass WFPC2 F606W 
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*2.415
    
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================


    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_wfpc2_f606w_image_w_aperture.png'.format(names_wfpc2_f606w[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['ATODGAIN']
    if gain != 7:
        if gain == 14:
            photflam += 1.987
            zero_pt += -2.5*np.log10(1.987)
    fluxes[i] = (phot_table['aperture_sum']/gain) # erg/s/cm^2/A
    f606w_mags[i] = - 2.5*np.log10(fluxes[i]) + zero_pt - A_v
    
    #pdb.set_trace()



# =============================================================================
# =============================================================================
# =============================================================================
# =====================       WFPC2 F547M        ==============================
# =============================================================================
# =============================================================================
# =============================================================================

#%%

names_wfpc2_f547m = []
RAs = np.zeros(len(filenames_wfpc2_f547m))
DECs = np.zeros(len(filenames_wfpc2_f547m))
fluxes = np.zeros(len(filenames_wfpc2_f547m))
f547m_mags = np.zeros(len(filenames_wfpc2_f547m))
print()
print()
print('WFPC2 F547M')
time.sleep(1)

for i in tqdm(range(len(filenames_wfpc2_f547m)), position=0, leave=True):
    
    hdul = fits.open(filenames_wfpc2_f547m[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    #cutind = 0
    #data = data[cutind:-cutind,cutind:-cutind]
   
    # determine galaxy name and specify center coordinates
    names_wfpc2_f547m.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-200):int(position[0]+200),int(position[1]-200):int(position[1]+200)]

   
   # pdb.set_trace()
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    #cen = np.where(data == np.max(data))
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 7.691e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    zero_pt = -21.662
    
    # determine galaxy name and specify center coordinates
    names_wfpc2_f547m.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass WFPC2 F606W 
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*2.755
    
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================
    
    

    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_wfpc2_f547m_image_w_aperture.png'.format(names_wfpc2_f547m[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['ATODGAIN']
    if gain != 7:
        if gain == 14:
            photflam += 1.987
            zero_pt += -2.5*np.log10(1.987)
    fluxes[i] = (phot_table['aperture_sum']/gain) # erg/s/cm^2/A
    f547m_mags[i] = - 2.5*np.log10(fluxes[i]) + zero_pt - A_v
    
    #pdb.set_trace()



# =============================================================================
# =============================================================================
# =============================================================================
# =====================        WFC3 F475W        ==============================
# =============================================================================
# =============================================================================
# =============================================================================
#%%

names_wfc3_f475w = []
RAs = np.zeros(len(filenames_wfc3_f475w))
DECs = np.zeros(len(filenames_wfc3_f475w))
fluxes = np.zeros(len(filenames_wfc3_f475w))
wfc3_f475w_mags = np.zeros(len(filenames_wfc3_f475w))
print()
print()
print('WFC3 F475W')
time.sleep(1)

for i in tqdm(range(len(filenames_wfc3_f475w)), position=0, leave=True):
    
    hdul = fits.open(filenames_wfc3_f475w[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    
    # let's trim some fat off the image
    data = hdul['SCI'].data
    cutind = 0
    ##data = data[cutind:-cutind,cutind:-cutind]
    
    # determine galaxy name and specify center coordinates
    names_wfc3_f475w.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-0.2*data.shape[0]):int(position[0]+0.2*data.shape[0]),int(position[1]-0.2*data.shape[1]):int(position[1]+0.2*data.shape[1])]

    
    #pdb.set_trace()
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    #cen = np.where(data == np.max(data))
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 3.483e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    zero_pt = -22.545	
    
    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass WFC3 F475W
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*3.248
       
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================


    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #                      norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_wfc3_f475w_image_w_aperture.png'.format(names_wfc3_f475w[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['CCDGAIN']
    EE = 0.92373
    photflam = 2.4984e-19
    zero_pt = 25.8094
    
    fluxes[i] = phot_table['aperture_sum']/EE # e-/s
    wfc3_f475w_mags[i] = - 2.5*np.log10(fluxes[i]) + zero_pt - A_v
    
    #pdb.set_trace()
    
# =============================================================================
# =============================================================================
# =============================================================================
# =====================        WFC3 F814W        ==============================
# =============================================================================
# =============================================================================
# =============================================================================
#%%

names_wfc3_f814w = []
RAs = np.zeros(len(filenames_wfc3_f814w))
DECs = np.zeros(len(filenames_wfc3_f814w))
fluxes = np.zeros(len(filenames_wfc3_f814w))
wfc3_f814w_mags = np.zeros(len(filenames_wfc3_f814w))
print()
print()
print('WFC3 F814W')
time.sleep(1)

for i in tqdm(range(len(filenames_wfc3_f814w)), position=0, leave=True):
    
    hdul = fits.open(filenames_wfc3_f814w[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    cutind = int(np.min((data.shape[0]*0.3, data.shape[1]*0.3)))
    #data = data[cutind:-cutind,cutind:-cutind]
    
    # determine galaxy name and specify center coordinates
    names_wfc3_f814w.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print(filenames_wfc3_f814w[i])
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-0.2*data.shape[0]):int(position[0]+0.2*data.shape[0]),int(position[1]-0.2*data.shape[1]):int(position[1]+0.2*data.shape[1])]

    
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    #cen = np.where(data == np.max(data))
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 3.483e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    zero_pt = -22.545	
    
    # determine galaxy name and specify center coordinates
    names_wfc3_f814w.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass WFC3 F814W
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*1.536
       
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================


    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_wfc3_f814w_image_w_aperture.png'.format(names_wfc3_f814w[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['CCDGAIN']
    EE = 0.9127
    photflam = 1.4980e-19
    zero_pt =  24.6985
    
    fluxes[i] = phot_table['aperture_sum']/EE # e-/s
    wfc3_f814w_mags[i] = - 2.5*np.log10(fluxes[i]) + zero_pt - A_v
    
    #pdb.set_trace()



# =============================================================================
# =============================================================================
# =============================================================================
# =====================        WFC3 F547M        ==============================
# =============================================================================
# =============================================================================
# =============================================================================
#%%

names_wfc3_f547m = []
RAs = np.zeros(len(filenames_wfc3_f547m))
DECs = np.zeros(len(filenames_wfc3_f547m))
fluxes = np.zeros(len(filenames_wfc3_f547m))
wfc3_f547m_mags = np.zeros(len(filenames_wfc3_f547m))
print()
print()
print('WFC3 F547M')
time.sleep(1)

for i in tqdm(range(len(filenames_wfc3_f547m)), position=0, leave=True):
    
    hdul = fits.open(filenames_wfc3_f547m[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    cutind = 0
    #data = data[cutind:-cutind,cutind:-cutind]
    
    # determine galaxy name and specify center coordinates
    names_wfc3_f547m.append(name_header['TARGNAME'])
    RAs[i] = name_header['RA_TARG']
    DECs[i] = name_header['DEC_TARG']
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-0.2*data.shape[0]):int(position[0]+0.2*data.shape[0]),int(position[1]-0.2*data.shape[1]):int(position[1]+0.2*data.shape[1])]

    
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    #cen = np.where(data == np.max(data))
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 3.483e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    zero_pt = -22.545	
    
    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass WFC3 F555W (F547M not available)
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*2.855
       
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================


    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_wfc3_f547m_image_w_aperture.png'.format(names_wfc3_f547m[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['CCDGAIN']
    
    EE = 0.92612
    photflam = 4.5959e-19
    zero_pt = 24.7583
    
    fluxes[i] = phot_table['aperture_sum']/EE # e-/s
    wfc3_f547m_mags[i] = - 2.5*np.log10(fluxes[i]) + zero_pt - A_v
    
    #pdb.set_trace()
    



# =============================================================================
# =============================================================================
# =============================================================================
# =====================         ACS F475W       ==============================
# =============================================================================
# =============================================================================
# =============================================================================
#%%

names_acs_f475w = []
RAs = np.zeros(len(filenames_acs_f475w))
DECs = np.zeros(len(filenames_acs_f475w))
fluxes = np.zeros(len(filenames_acs_f475w))
acs_f475w_mags = np.zeros(len(filenames_acs_f475w))

print()
print()
print('FETCHING ACS ZEROPOINT...')
time.sleep(1)
w = acszpt.Query(date='2017-01-01', detector='WFC', filt='F475W')
zero_pt = w.fetch()
print()
print('ACS F475W')
time.sleep(1)

for i in tqdm(range(len(filenames_acs_f475w)), position=0, leave=True):
    
    hdul = fits.open(filenames_acs_f475w[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    cutind = 0
    #data = data[cutind:-cutind,cutind:-cutind]
    
    # determine galaxy name and specify center coordinates
    names_acs_f475w.append(name_header['TARGNAME'])
    
    acs475_ra = np.array([184.87750, 190.697583, 45.39958, 186.902841, 187.7176935, 51.61812, 187.33908, 192.14875, 187.6693575])
    acs475_dec = np.array([14.8775607, 11.4424897, -14.83639, 8.1543168, 16.7588163, -21.35517, 8.15642, -5.80083, 9.0156651])
    
    RAs[i] = acs475_ra[i]
    DECs[i] = acs475_dec[i]
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-0.2*data.shape[0]):int(position[0]+0.2*data.shape[0]),int(position[1]-0.2*data.shape[1]):int(position[1]+0.2*data.shape[1])]

    
    #pdb.set_trace()
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    #cen = np.where(data == np.max(data))
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 3.483e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    #zero_pt = -22.545	
    
# =============================================================================
#     w = acszpt.Query(date='2017-01-01', detector='WFC', filt='F475W')
#     zero_pt = w.fetch()
# =============================================================================
    
    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass ACS F475W 
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*3.268
       
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================


    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_acs_f475w_image_w_aperture.png'.format(names_acs_f475w[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['CCDGAIN']
    correction_inf = 0.912 # 0.5" aperture
    
    flux = (phot_table['aperture_sum']/gain) # erg/s/cm^2/A
    flux_inf = flux / correction_inf

    # Now convert instrumental fluxes to physical fluxes and magnitudes.
    # F_lambda is the flux density in units of erg/sec/cm^2/Angstrom.
    F_lambda = flux_inf * zero_pt['PHOTFLAM']

    
    fluxes[i] = flux_inf # erg/s/cm^2/A
    acs_f475w_mags[i] = -2.5 * np.log10(flux_inf) + zero_pt['VEGAmag'][0].value - A_v
    
    #pdb.set_trace()

# =============================================================================
# =============================================================================
# =============================================================================
# =====================        ACS F850LP        ==============================
# =============================================================================
# =============================================================================
# =============================================================================
#%%

names_acs_f850lp = []
RAs = np.zeros(len(filenames_acs_f850lp))
DECs = np.zeros(len(filenames_acs_f850lp))
fluxes = np.zeros(len(filenames_acs_f850lp))
acs_f850lp_mags = np.zeros(len(filenames_acs_f850lp))

print()
print()
print('FETCHING ACS ZEROPOINT...')
time.sleep(1)
w = acszpt.Query(date='2017-01-01', detector='WFC', filt='F850LP')
zero_pt = w.fetch()
print()
print('ACS F850LP')
time.sleep(1)

for i in tqdm(range(len(filenames_acs_f850lp)), position=0, leave=True):
    
    hdul = fits.open(filenames_acs_f850lp[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    cutind = 0
    #data = data[cutind:-cutind,cutind:-cutind]
    
    # determine galaxy name and specify center coordinates
    names_acs_f850lp.append(name_header['TARGNAME'])
    
    acs850_ra = np.array([186.902841, 45.39958, 187.7176935, 190.697583, 51.61812, 184.877508, 192.14875, 187.6693575, 187.33908])
    acs850_dec = np.array([8.1543168, -14.83639, 16.7588163, 11.4424897, -21.35517, 14.8775607, -5.80083, 9.0156651, 8.15642])
    
    RAs[i] = acs850_ra[i]
    DECs[i] = acs850_dec[i]
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    data = data_all[int(position[0]-0.2*data.shape[0]):int(position[0]+0.2*data.shape[0]),int(position[1]-0.2*data.shape[1]):int(position[1]+0.2*data.shape[1])]

    
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    #cen = np.where(data == np.max(data))
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5 
    data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.1*n_x),int(cen[1]+0.1*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.1*n_y),int(cen[0]+0.1*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 3.483e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    #zero_pt = -22.545	
    
# =============================================================================
#     w = acszpt.Query(date='2017-01-01', detector='WFC', filt='F850LP')
#     zero_pt = w.fetch()
# =============================================================================

    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass ACS F850LP 
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*1.243
       
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================

    

    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_acs_f850lp_image_w_aperture.png'.format(names_acs_f850lp[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['CCDGAIN']
    correction_inf = 0.893 # 0.5" aperture
    
    flux = (phot_table['aperture_sum']/gain) # erg/s/cm^2/A
    flux_inf = flux / correction_inf

    # Now convert instrumental fluxes to physical fluxes and magnitudes.
    # F_lambda is the flux density in units of erg/sec/cm^2/Angstrom.
    F_lambda = flux_inf * zero_pt['PHOTFLAM']

    
    fluxes[i] = flux_inf # erg/s/cm^2/A
    acs_f850lp_mags[i] = -2.5 * np.log10(flux_inf) + zero_pt['VEGAmag'][0].value - A_v
    
# =============================================================================
#     if i == 4:
#         pdb.set_trace()
# =============================================================================



# =============================================================================
# =============================================================================
# =============================================================================
# =====================         ACS F814W        ==============================
# =============================================================================
# =============================================================================
# =============================================================================
#%%

names_acs_f814w = []
RAs = np.zeros(len(filenames_acs_f814w))
DECs = np.zeros(len(filenames_acs_f814w))
fluxes = np.zeros(len(filenames_acs_f814w))
acs_f814w_mags = np.zeros(len(filenames_acs_f814w))

print()
print()
print('FETCHING ACS ZEROPOINT...')
time.sleep(1)
w = acszpt.Query(date='2017-01-01', detector='WFC', filt='F814W')
zero_pt = w.fetch()
print()
print('ACS F814W')
time.sleep(1)

for i in tqdm(range(len(filenames_acs_f814w)), position=0, leave=True):
    
    hdul = fits.open(filenames_acs_f814w[i])
    name_header = hdul[0].header
    header = hdul['SCI'].header
    data_all = hdul['SCI'].data
    # let's trim some fat off the image
    data = hdul['SCI'].data
    cutind = 0
    #data = data[cutind:-cutind,cutind:-cutind]
    
    # determine galaxy name and specify center coordinates
    names_acs_f814w.append(name_header['TARGNAME'])
    
    acs814_ra = np.array([214.6075, 53.81917])
    acs814_dec = np.array([36.49331, -35.22611])
    
    RAs[i] = acs814_ra[i]
    DECs[i] = acs814_dec[i]
    
    # use the header to manually compute the pixel size for the image
    RA_ref = header['CRVAL1']
    DEC_ref = header['CRVAL2']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']
    CD1_1 = header['CD1_1'] # partial of the right ascension w.r.t. x
    CD1_2 = header['CD1_2'] # partial of the right ascension w.r.t. y
    CD2_1 = header['CD2_1'] # partial of the declination w.r.t. x    
    CD2_2 = header['CD2_2'] # partial of the declination w.r.t. y
    pixel_size = np.sqrt((CD1_1+CD2_1)**2) * 3600
    arcsec_to_pix = 1/pixel_size
    
    if pixel_size < 0.049 or pixel_size > 0.051:
        print()
        print('WARNING: Possible pixel size error.')
        print('Pixel Size: {:.3f}'.format(pixel_size))
        print()
    
    position = np.zeros(2)
    position[0] = int(x_ref + (RA_ref-RAs[i])/CD1_1)  
    position[1] = int(y_ref + (DEC_ref-DECs[i])/CD2_2)

    #pdb.set_trace()
    # manually edit position
    if i == 1:
        position = position[::-1]

    data = data_all[int(position[0]-0.2*data.shape[0]):int(position[0]+0.2*data.shape[0]),int(position[1]-0.2*data.shape[1]):int(position[1]+0.2*data.shape[1])]
    
# =============================================================================
#     plt.figure(dpi=500)
#     plt.imshow(data, cmap='gray', vmin=0, vmax=2)
#     plt.colorbar()
#     plt.show()
# =============================================================================
    
    #cen = np.argmax(data)
    
    n_pix = data.shape[0]*data.shape[1]
    n_x = data.shape[0]
    n_y = data.shape[1]
    n_sec = 8
    sums_25 = np.zeros((n_sec,n_sec))
    
    # first mask the outliers
    sig = np.std(data.flatten())
    med = np.median(data.flatten())
    sigfac = 5
    data_temp = np.zeros_like(data)
    #data_temp = np.zeros_like(data)
    data_temp = np.copy(data)
    data_temp[np.where(data_temp > med+sigfac*sig)] = med+sigfac*sig
    data_temp[np.where(data_temp < med-sigfac*sig)] = med-sigfac*sig
    for v in range(n_sec):
        for w in range(n_sec):
            sums_25[v,w] = np.sum(data_temp[int(v/n_sec*n_x):int((v+1)/n_sec*n_x),int(w/n_sec*n_y):int((w+1)/n_sec*n_y)])
    sector = np.where(sums_25 == np.max(sums_25))[::-1] # reverse array so x and y are correct
    
    cen = np.zeros(2)
    m = np.where(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)] 
                 == np.max(data[int(sector[0]/n_sec*n_x):int((sector[0]+1)/n_sec*n_x),int(sector[1]/n_sec*n_y):int((sector[1]+1)/n_sec*n_y)]))
    cen[0] = m[0][0]+int(sector[0]/n_sec*n_x)
    cen[1] = m[1][0]+int(sector[1]/n_sec*n_y)
    
    x_range = np.array([int(cen[1]-0.2*n_x),int(cen[1]+0.2*n_x)])
    if x_range[0] < 0:
        x_range[0] = 0
    if x_range[1] >= n_x:
        x_range[1] = n_x
    y_range = np.array([int(cen[0]-0.2*n_y),int(cen[0]+0.2*n_y)])
    if y_range[0] < 0:
        y_range[0] = 0
    if y_range[1] >= n_y:
        y_range[1] = n_y

    data_zoom = data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    cen_zoom = centroid_2dg(data_zoom)
    
    med_bkg = np.median(data_all[0:50,0:50].flatten())
    mode = stats.mode(data_all.flatten())[0]
    dist = data_all.flatten()-mode
    std = np.std(dist)
    if med_bkg <= mode-std:
        print('WARNING: Possible background subtraction issue...')
        pdb.set_trace()
    #pdb.set_trace()
# =============================================================================
#     
#     # store important information from the header
#     #zero_pts[i] = header['PHOTZPT']
#     x_ref = header['CRPIX1']
#     y_ref = header['CRPIX2']
#     RA_ref = header['CRVAL1'] # deg
#     DEC_ref = header['CRVAL2'] # deg
#     dRA_dx = header['CD1_1']
#     dRA_dy = header['CD1_2']
#     dDec_dx = header['CD2_1']
#     dDec_dy = header['CD2_2']
# =============================================================================
    #photflam = header['PHOTFLAM']
    photflam = 3.483e-18
    exptime = name_header['EXPTIME']
    #photplam = header['PHOTPLAM']
    #zero_pt = header['PHOTZPT']
    #zero_pt = -22.545	
    
# =============================================================================
#     w = acszpt.Query(date='2017-01-01', detector='WFC', filt='F814W')
#     zero_pt = w.fetch()
# =============================================================================
    
    coords = SkyCoord(ra=RAs[i]*u.degree, dec=DECs[i]*u.degree)
    sfd = SFDQuery()
    E_BV = sfd(coords)
    #A_v = 3.1*E_BV
    # correct the A_v value for the appropriate bandpass ACS F814W 
    # (value from Table 6 in Schlafly & Finkbeiner (2011))
    A_v = E_BV*1.526
       
# =============================================================================
#     # specify galaxy center in pixel coordinates
#     if dDec_dx == 0:
#         x_step = (RAs[i]-RA_ref)/dRA_dx
#     else:
#         x_step = (RAs[i]-RA_ref)/dRA_dx + (DECs[i]-DEC_ref)/dDec_dx
#     if dRA_dy == 0:
#         y_step = (DECs[i]-DEC_ref)/dDec_dy
#     else: 
#         y_step = (RAs[i]-RA_ref)/dRA_dy + (DECs[i]-DEC_ref)/dDec_dy
#     
#     x_cen = x_ref + x_step
#     y_cen = y_ref + y_step
#     #position = (x_cen,y_cen)
#     pdb.set_trace()
# =============================================================================

# =============================================================================
#     wcs = WCS(header)  
#     positions = SkyCoord("{}{}".format(RAs[i],DECs[i]), unit=(u.deg, u.deg), frame='icrs')
#     #positions = SkyCoord(cRAs[i], DECs[i], frame='galactic')  
#     aperture = SkyCircularAperture(positions, r=0.5 * u.arcsec)  
# =============================================================================


    cen_tot = (int(cen[0]-0.2*n_x)+cen_zoom[0],int(cen[1]-0.2*n_y)+cen_zoom[1])
    position = cen_zoom
    #arcsec_to_pix = 1/0.05
    rad_arcsec = 0.5
    rad_pix = arcsec_to_pix * rad_arcsec
    aperture = CircularAperture(position, r=rad_pix)
    
    plt.figure(dpi=500)
    plt.imshow(data_zoom, norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.imshow(data[int(cen_tot[1]-0.2*n_x):int(cen_tot[1]+0.2*n_x),int(cen_tot[0]-0.2*n_y):int(cen_tot[0]+0.2+n_y)],
    #           norm=colors.LogNorm(vmin=0.01,vmax=100),cmap='Greys')
    #plt.plot(cen[0],cen[1],linestyle='',marker='.',color='b')
    circ = plt.Circle(cen_zoom,aperture.r,alpha=0.3,color='r')
    plt.gca().add_patch(circ)
    #plt.plot(137,776,linestyle='',marker='.',color='r')
    plt.colorbar()
    #plt.show()
    
    plt.savefig('../Plots/aperture_photometry/{}_acs_f814w_image_w_aperture.png'.format(names_acs_f814w[i]),
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
    phot_table = aperture_photometry(data_zoom, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    
    # get data from e/s to DN/s to compute magnitude
    gain = name_header['CCDGAIN']
    correction_inf = 0.914 # 0.5" aperture
    
    flux = (phot_table['aperture_sum']/gain) # erg/s/cm^2/A
    flux_inf = flux / correction_inf

    # Now convert instrumental fluxes to physical fluxes and magnitudes.
    # F_lambda is the flux density in units of erg/sec/cm^2/Angstrom.
    F_lambda = flux_inf * zero_pt['PHOTFLAM']

    
    fluxes[i] = flux_inf # erg/s/cm^2/A
    acs_f814w_mags[i] = -2.5 * np.log10(flux_inf) + zero_pt['VEGAmag'][0].value - A_v
    
    #pdb.set_trace()


print()
print()
print("======================================")
print("==== APERTURE PHOTMETRY FINISHED =====")
print("======================================")

#%%

# Let's make a file to store all of this info
import csv

with open('../Result_Tables/nuclear_mags.csv', 'w') as f:
    # create the csv writer
    writer = csv.writer(f)
    # write a row to the csv file
    writer.writerow('# WFPC2 F555W')
    writer.writerow(names_wfpc2_f555w)
    writer.writerow(f555w_mags)

    writer.writerow('# WFPC2 F814W')
    writer.writerow(names_wfpc2_f814w)
    writer.writerow(f814w_mags)
    
    writer.writerow('# WFPC2 F606W')
    writer.writerow(names_wfpc2_f606w)
    writer.writerow(f606w_mags)
    
    writer.writerow('# WFPC2 F547M')
    writer.writerow(names_wfpc2_f547m)
    writer.writerow(f547m_mags)
    
    writer.writerow('# WFC3 F475W')
    writer.writerow(names_wfc3_f475w)
    writer.writerow(wfc3_f475w_mags)
    
    writer.writerow('# WFC3 F814W')
    writer.writerow(names_wfc3_f814w)
    writer.writerow(wfc3_f814w_mags)
    
    writer.writerow('# WFC3 F547M')
    writer.writerow(names_wfc3_f547m)
    writer.writerow(wfc3_f547m_mags)
    
    writer.writerow('# ACS F475W')
    writer.writerow(names_acs_f475w)
    writer.writerow(acs_f475w_mags)
    
    writer.writerow('# ACS F850LP')
    writer.writerow(names_acs_f850lp)
    writer.writerow(acs_f850lp_mags)
    
    writer.writerow('# ACS F814W')
    writer.writerow(names_acs_f814w)
    writer.writerow(acs_f814w_mags)


#%%

other_delete_names = np.array(['NGC 4239', 'NGC 2549', 'NGC 2699'])

names_95 = np.array(['NGC 2841', 'NGC 5845', 'NGC 2636', 'NGC 4742', 'NGC 4697',
           'VCC 1440', 'NGC 1172', 'NGC 4636', 'NGC 7332', 'NGC 4239',
           'NGC 0524', 'NGC 4464', 'VCC 1545', 'NGC 0720', 'NGC 4467',
           'VCC 1627', 'VCC 1199', 'NGC 3605', 'NGC 3599', 'NGC 1331',
           'NGC 1400'])
lauer_95_dists = np.array([14.6       , 31.5       , 29.33983225, 15.5       , 11.2       ,
       15.99558   , 21.6       , 16.5       , 22.2       , 16.5       ,
       27.1       , 15.848925  , 16.826734  , 23.7       , 16.5       ,
       15.631468  , 16.5       , 21.7       , 20.8       , 21.1       ,
       24.89      ])

inner_rad_95 = np.array([1.55720972, 3.3597333 , 3.1293337 , 1.6532021 , 1.19457184,
       1.70605977, 2.30381712, 1.7598603 , 2.36781204, 1.7598603 ,
       2.89043722, 1.69041781, 1.79470916, 2.52779934, 1.7598603 ,
       1.66722424, 1.7598603 , 2.31448294, 2.21849056, 2.25048802,
       2.6547226 ])
logmass_95 = np.array([10.83419286, 10.23663002,  9.63350373,  9.93199102, 10.62223233,
        8.83220232, 10.12668444, 10.75267158, 10.36737072,  9.35366164,
       11.0179657 ,  9.57993237,  8.84774733, 10.86377766,  9.04294519,
        8.59608566,  8.80372546,  9.793021  ,  9.98065749,  9.38221741,
       10.64268276])

nsc_comp_95 = np.array([-1,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
       -1, -1, -1, -1])



arcsec_to_radian = 4.8481e-6
p1arcs = 0.1*arcsec_to_radian*lauer_95_dists*10**6


p = np.where(p1arcs < 10)[0]

k = np.where(logmass_95[p] <= 10.5)[0]



plt.figure(dpi=500)
plt.title('21 Lauer+95 WFPC1 Galaxies')
a = plt.hist(p1arcs,bins=20)
plt.hist(inner_rad_95,bins=20, color='0.8')
#plt.plot([5,5],[0,np.max(a[0])], linestyle='--', color='k')
plt.plot([10,10],[0,np.max(a[0])], linestyle='--', color='r')
plt.xlabel('pc radius corresponding to 0.1"')


plt.figure(dpi=500)
plt.title('21 Lauer+95 WFPC1 Galaxies')
plt.hist(logmass_95,bins=30)
plt.xlabel('log(M$_{*}$ [M$_{\odot}$])')

#%%

# Lauer+95 galaxies that are low mass (< 10.5), have nuclear color info, 
# and 0.1" is within the 10pc limit
names_95_c = np.array(['VCC 1440', 'NGC 4464', 'VCC 1545',
       'NGC 4467', 'VCC 1627', 'VCC 1199'])
logmass_95_c = np.array([8.83220232, 9.57993237, 8.84774733,
       9.04294519, 8.59608566, 8.80372546])
p1arcs_95_c = np.array([7.75481714, 7.68371733, 8.15776891,
       7.999365  , 7.578292  , 7.999365  ])
inner_rad_95_c = np.array([1.70605977, 1.69041781, 1.79470916,
       1.7598603 , 1.66722424, 1.7598603 ])


plt.figure(dpi=500)
plt.title('21 Lauer+95 WFPC1 Galaxies (GOOD)')
a = plt.hist(p1arcs_95_c,bins=20)
plt.hist(inner_rad_95_c,bins=20, color='0.8')
#plt.plot([5,5],[0,np.max(a[0])], linestyle='--', color='k')
plt.plot([10,10],[0,np.max(a[0])], linestyle='--', color='r')
plt.xlabel('pc radius corresponding to 0.1"')


plt.figure(dpi=500)
plt.title('21 Lauer+95 WFPC1 Galaxies (GOOD)')
plt.hist(logmass_95_c,bins=30)
plt.xlabel('log(M$_{*}$ [M$_{\odot}$])')

