#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 19:15:56 2022

Code to search the Véron-Cetty & Véron (2006) AGN Catalog for our sample galaxies.

@author: christian
"""

import numpy as np
from astropy.io import fits
import astropy.units as uni
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.table import Table, join, vstack, Column, MaskedColumn
import matplotlib.colors as colors
from astropy.coordinates import Angle
import pdb
from astropy.table import Table

def HMS2deg(ra='', dec=''):
  RA, DEC, rs, ds = '', '', 1, 1
  if dec:
    D, M, S = [float(i) for i in dec.split()]
    if str(D)[0] == '-':
      ds, D = -1, abs(D)
    deg = D + (M/60) + (S/3600)
    DEC = '{0}'.format(deg*ds)
  
  if ra:
    H, M, S = [float(i) for i in ra.split()]
    if str(H)[0] == '-':
      rs, H = -1, abs(H)
    deg = (H*15) + (M/4) + (S/240)
    RA = '{0}'.format(deg*rs)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC



gal_file = './Data_Sets/bi20_xray_table.fit'
hdul = fits.open(gal_file)
head = hdul[0].data
data = hdul[1].data
hdul.close()
#pdb.set_trace()
ra3 = np.array(data['RAJ2000'])
dec3 = np.array(data['DEJ2000'])

for i in range(len(ra3)):
    ra3[i], dec3[i] = HMS2deg(ra=ra3[i],dec=dec3[i])
    
x_catalog = SkyCoord(ra=ra3, dec=dec3, unit=uni.degree)



gal_file = '../Data_Sets/veron_cetty_agn_catalog_2010.fit'
hdul = fits.open(gal_file)
head = hdul[0].data
dat = hdul[1].data
hdul.close()


"""
Spectral Types for Veron-Cetty Catalog:
Note (3)  : the classification of the object is Q2 : type-2 quasar S1 : Seyfert 1 spectrum

S1h: broad polarized Balmer lines detected
S1i: broad Paschen lines observed in the infrared
S1n: narrow-line Seyfert 1
S1.0, S1.2, S1.5, S1.8, S1.9: intermediate Seyfert galaxies
S2 : Seyfert 2 spectrum
S2?: probable Seyfert 2
S3 : Seyfert 3 or liner
S3b: Seyfert 3 or liner with broad Balmer lines
S3h: Seyfert 3 or liner with broad polarized Balmer lines detected
S : unclassified Seyfert
S? : possibly a Seyfert
H2 : nuclear HII region
HP : high optical polarization (>3%)
BL : confirmed BL Lac object
BL?: probable BL Lac object
? : questionable BL Lac object
"""


ra2 = np.array(dat['RAJ2000'])
dec2 = np.array(dat['DEJ2000'])

for i in range(len(ra2)):
    ra2[i], dec2[i] = HMS2deg(ra=ra2[i],dec=dec2[i])
    
catalog = SkyCoord(ra=ra2, dec=dec2, unit=uni.degree)

# read in fits table with model galaxy parameters
gal_filename = '../Result_Tables/final_gal_data.fits'
hdul = fits.open(gal_filename)
gal_data = hdul[1].data 
hdul.close()

names = gal_data['name']
RAs = gal_data['RA']
DECs = gal_data['DEC']
types = gal_data['type']
gammas = gal_data['slope']
log_rhos = gal_data['cen_dens']
gal_logmass = gal_data['logmass']
dists = gal_data['dist']

num_agn = 0 
num_xray = 0
num_bi = 0
#%%
partial_match_inds = []
for i in range(len(ra2)):
    if '1399' in dat['Name'][i]:
        partial_match_inds.append(i)

#%%
import pdb

spec_class = np.zeros_like(names).astype(str)
xray_lum = np.zeros_like(gammas)-99.
xray_names = np.zeros_like(names)
spec_names = np.zeros_like(names)

for i in range(len(names)):
    ra1 = RAs[i]
    dec1 = DECs[i]
    
    c = SkyCoord(ra=ra1*uni.degree, dec=dec1*uni.degree)
    idx, d2d, d3d = match_coordinates_sky(c,catalog)
    bidx, bd2d, bd3d = match_coordinates_sky(c,x_catalog)
    
    
# =============================================================================
#     if i == 28:
#         pdb.set_trace()
# =============================================================================

    max_sep_1 = (1./60.) * uni.degree
    max_sep_2 = (1./10.) * uni.degree    
    
    if d2d < max_sep_1: 
        print('Veron-Cetty+10 AGN found...')
        print(dat['Name'][idx])
        print(names[i])
        print('Catalog Index: ', idx)
        print()
        spec_class[i] = dat['Sp'][idx]
        spec_names[i] = dat['Name'][idx]
        num_agn += 1

    if bd2d < max_sep_2: 
        num_bi += 1
        print('Bi+20 match found...')
        print(data['Name'][bidx])
        print(names[i])
        print('Catalog Index: ', bidx)
        print('log(Lx) = {:.2f}'.format(data['logLX'][bidx])) #log(10^7 Watts)
        print()
        xray_lum[i] = data['logLX'][bidx]
        xray_names[i] = data['Name'][bidx]
        num_xray += 1
        #pdb.set_trace()



# sort the data
a = np.argsort(names)
spec_names = spec_names[a]
xray_names = xray_names[a]
names = names[a]
xray_lum = xray_lum[a]
spec_class = spec_class[a]
gal_logmass = gal_logmass[a]
RAs = RAs[a]
DECs = DECs[a]
types = types[a]
gammas = gammas[a]
log_rhos = log_rhos[a]
dists = dists[a]
            
c1 = fits.Column(name='name', array=names, format='10A')
c2 = fits.Column(name='agn_spec_class', array=spec_class, format='5A')
c3 = fits.Column(name='xray_lum', array=xray_lum, format='D', unit='erg/s')
c4 = fits.Column(name='logmass', array=gal_logmass, format='D', unit='M_sol')
c5 = fits.Column(name='RA', array=RAs, format='D', unit='deg')
c6 = fits.Column(name='Dec', array=DECs, format='D', unit='deg')
c7 = fits.Column(name='type', array=types, format='I', unit='0=et,1=lt')
c8 = fits.Column(name='slope', array=gammas, format='D')
c9 = fits.Column(name='cen_dens', array=10**log_rhos, format='D', unit='M_sol/pc^3')
c10 = fits.Column(name='dist', array=dists, format='D', unit='Mpc')


t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])
t.writeto('../Result_Tables/gal_agn_table.fits',overwrite=True)   

#%%
# write as a CSV file too
table = Table.read('./gal_agn_table.fits')
table.write('../Result_Tables/gal_agn_table.csv', overwrite=True)



