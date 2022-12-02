#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 03:15:28 2022

Code to search for colors of our galaxies in Georgiev+14 data

@author: christian
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import pdb



model_gal_filename = '/Users/christian/OneDrive/Desktop/TDE Code/final_gal_data.fits'
hdul = fits.open(model_gal_filename)
gal_data = hdul[1].data 
hdul.close()

gal_names = gal_data['name']

georgiev_filename = '/Users/christian/OneDrive/Desktop/TDE Code/georgiev_14_data.fit'
hdul = fits.open(georgiev_filename)
g14_data = hdul[1].data 
hdul.close()

g_names = g14_data['Name']

f814w_mags = g14_data['m2']
f450w_mags = g14_data['m3']
f555w_mags = g14_data['m4']



num = '1426'
n = 'NGC'+num
found = False
for i in range(len(g_names)):
    if n == g_names[i]:
        print('Match found!')
        print()
        print(g_names[i])
        print('F814W: ',f814w_mags[i])
        print('F450W: ',f450w_mags[i])
        print('F555W: ',f555w_mags[i])
        found = True
if not found:
    print('No luck...')