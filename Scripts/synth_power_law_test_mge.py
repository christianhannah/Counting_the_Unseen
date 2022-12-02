#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code serves to implement a power law slope flattening test from MGE using
a synthetic power-law SB profile. 

Created on Sun Jan 30 20:20:04 2022

@author: christian
"""

from mgefit.mge_fit_1d import mge_fit_1d
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import mge1d_util as u


x = np.linspace(0.02,50,20000)
y = x**(-1.5)



plt.figure(dpi=500)
plt.plot(x,y,marker='.')
plt.xlabel('Radius')
plt.ylabel('Synthetic Surface Brightness')
plt.xlim(0,1)

plt.figure(dpi=500)
plt.plot(np.log10(x),np.log10(y),marker='.')
plt.xlabel('log(Radius)')
plt.ylabel('log(Synthetic Surface Brightness)')

# ensure the radii are logarithmically spaced 
num_rad = 5000
rad = np.geomspace(x[0], x[-1], num_rad)
SBs = 10**np.interp(np.log10(rad), np.log10(x), np.log10(y))
#f1 = interpolate.interp1d(np.log10(x),np.log10(y),kind='cubic')
#SBs = 10**f1(np.log10(rad))

plt.figure(dpi=500)
plt.plot(x,y)
plt.plot(rad,SBs,marker='+',linestyle='')
plt.xlabel('Radius')
plt.ylabel('Synthetic Surface Brightness')
plt.xlim(0,1)

plt.figure(dpi=500)
plt.plot(np.log10(x),np.log10(y),marker='.')
plt.plot(np.log10(rad),np.log10(SBs))#,marker='+')#,linestyle='')
plt.xlabel('log(Radius)')
plt.ylabel('log(Synthetic Surface Brightness)')
plt.xlim(-1.73,-1.6)
plt.ylim(2.37,2.57)
plt.text(-1.72,2.4,str(num_rad)+' points')

print(u.get_density_and_slope(np.log10(x), np.log10(y)))
print(u.get_density_and_slope(np.log10(rad), np.log10(SBs)))


# perform MGE fit
mge = mge_fit_1d(rad, SBs, ngauss=20, plot=True)

# define array to store radial densities
densities = np.zeros_like(rad)

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
        sig = (sig_i)#*arcsec_to_radian)*distance # in pc
        L_tot = (height)*2*np.pi*(sig)**2 # total area under 2-d gaussian
        
        if j == 0:
            tots.append(L_tot_i)
            heights.append(height)
            sigs.append(sig_i)
            
        # convert total luminosity to total mass using M/L ratio
        M_tot = L_tot#*M_L
        # compute desity contribution from one gaussian component
        dens_comp = (M_tot/((2*np.pi)**(3/2)*sig**3))*np.exp(-(rad[j]**2/(2*sig**2))) # M_sol/pc^3
        # add density from each gaussian to total density
        #print(dens_comp, radii[j])
        densities[j] += dens_comp
    #pdb.set_trace()

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
plt.ylim(min(SBs),max(SBs)+500)
plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
plt.xlabel('Radius ["]')


real_slope, real_density, blah, blah = u.get_density_and_slope(np.log10(rad),np.log10(SBs))
mge_slope, mge_density, blah, blah = u.get_density_and_slope(np.log10(rad),np.log10(summed_gauss))

plt.text(0.62, 0.8, 'Data Slope = {:.2f}'.format(real_slope), 
         transform=plt.gcf().transFigure, size=10)
plt.text(0.62, 0.75, 'MGE Slope = {:.2f}'.format(mge_slope), 
         transform=plt.gcf().transFigure, size=10)



