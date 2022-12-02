#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 18:12:31 2021

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
import os
from matplotlib.collections import LineCollection

# =============================================================================
# read in our data for our galaxies

gal_file = '/Users/christian/OneDrive/Desktop/TDE Code/all_gal_data_100000pc.fits'
hdul = fits.open(gal_file)
head = hdul[0].data
dat = hdul[1].data
hdul.close()

names = dat['name']
vmags = dat['vmag']
dists = dat['dist'] # Mpc
MLs = dat['ml']
ML_types = dat['ml_type']
slopes = dat['slope']
cen_dens = dat['dens_at_5pc'] # M_sol/pc^3
lograds = dat['lograd'] # log(pc)
logdens = dat['logdens'] # log(M_sol/pc^3)

# =============================================================================
# =============================================================================


# =============================================================================
# plot all density profiles together by galaxy v-band magnitude within the range 
# of v-mags and inner radii allowed

ls = []
vmags_adj = []
slopes_adj = []
names_adj = []
for i in range(len(names)):
    a = np.zeros((len(lograds[0,:]),2))
    a[:,0] = lograds[i,:]
    a[:,1] = logdens[i,:]
    if vmags[i] > -29 and a[0,0] <= np.log10(5):
        vmags_adj.append(vmags[i])
        slopes_adj.append(slopes[i])
        names_adj.append(names[i])
        ls.append(a)
vmags_adj = np.array(vmags_adj)
slopes_adj = np.array(slopes_adj)
names_adj = np.array(names_adj)

fig, ax = plt.subplots(dpi=500)
lines = LineCollection(ls, array=vmags_adj, cmap='viridis', alpha=0.7)
ax.add_collection(lines)
fig.colorbar(lines, label='M$_{v}$ [mag]')
ax.set_xlabel('log(Radius [pc])')
ax.set_ylabel('log(Density [M$_{\odot}$/pc$^3$])')
ax.autoscale()
#ax.set_xlim(-0.2,1)
#ax.set_ylim(1,7)
plt.show()
plt.close()


# plot the distribution of slopes againts v-mag
plt.figure(dpi=500)
plt.plot(vmags_adj, slopes_adj, marker='.', linestyle='',color='b')
plt.xlabel('M$_{v}$ [mag]')
plt.ylabel('Inner Power-law Slope')

#%%

# plot the distribution based on innermost radii and v-mag with accepted cutoffs
plt.figure(dpi=500)
plt.plot(vmags,10**lograds[:,0],marker='.',linestyle='',color='b')
plt.plot([min(vmags),max(vmags)], [5,5], linestyle='--', color='r')
plt.plot([-21,-21], [0,35], linestyle='--', color='g')
plt.xlabel('M$_{v}$ [mag]')
plt.ylabel('Innermost Radius [pc]')


#%%

#colors = ['r',(.8,.1,.0),(.8,.5,.0),(.8,.7,.0),'y',(.5,.7,.0),(.2,.8,0.),(.0,.9,.0),(.0,.8,.4),
#              (.0,.8,.6),'c',(0,.5,.6),(.0,.3,.7),'b',(.2,.0,.8),(.5,0,.8),(.7,0,.7),(.9,.0,.8)]
colors = ['r',(.8,.5,.0),'y',(.2,.8,0.),
              'c',(.0,.3,.7),(.2,.0,.8),(.7,0,.7)]

# plot the distribution based with magnitude bins selected
plt.figure(dpi=500)
fig, ax = plt.subplots(dpi=500)
plt.plot(vmags,10**lograds[:,0],marker='.',linestyle='',color='b')
#plt.plot([min(vmags),max(vmags)], [5,5], linestyle='--', color='r')
ax.axvspan(-15,-17,alpha=0.4,color=colors[0])
ax.axvspan(-17,-18,alpha=0.4,color=colors[1])
ax.axvspan(-18,-19,alpha=0.4,color=colors[2])
ax.axvspan(-19,-20.3,alpha=0.4,color=colors[3])
ax.axvspan(-20.3,-21.2,alpha=0.4,color=colors[4])
ax.axvspan(-21.2,-22,alpha=0.4,color=colors[5])

# =============================================================================
# plt.plot([-15,-15], [0,35], linestyle='--', color='g')
# plt.plot([-17,-17], [0,35], linestyle='--', color='g')
# plt.plot([-18,-18], [0,35], linestyle='--', color='g')
# plt.plot([-19,-19], [0,35], linestyle='--', color='g')
# plt.plot([-20.3,-20.3], [0,35], linestyle='--', color='g')
# plt.plot([-21.2,-21.2], [0,35], linestyle='--', color='g')
# plt.plot([-22,-22], [0,35], linestyle='--', color='g')
# =============================================================================
plt.ylim(0,5)
plt.xlabel('M$_{v}$ [mag]')
plt.ylabel('Innermost Radius [pc]')

# =============================================================================
# =============================================================================


#%%
# =============================================================================
# let's bin the density profiles together based on luminosity
# AND LET'S EXCLUDE THE PROFILES WITH NO COVERAGE TO 5 PC



# define the bins
bin_1 = np.array([-17,-15])
bin_2 = np.array([-18,bin_1[0]+.00001])
bin_3 = np.array([-19,bin_2[0]+.00001])
bin_4 = np.array([-20.3,bin_3[0]+.00001])
bin_5 = np.array([-21.2,bin_4[0]+.00001])
bin_6 = np.array([-22,bin_5[0]+.00001])

# get the indices of glaxies in each bin
inds_1 = []
inds_2 = []
inds_3 = []
inds_4 = []
inds_5 = []
inds_6 = []
for i in range(len(names)):
    if vmags[i] > bin_1[0] and vmags[i] <= bin_1[1] and lograds[i,0] <= np.log10(5):
        inds_1.append(i)
    elif vmags[i] > bin_2[0] and vmags[i] <= bin_2[1] and lograds[i,0] <= np.log10(5):
        inds_2.append(i)
    elif vmags[i] > bin_3[0] and vmags[i] <= bin_3[1] and lograds[i,0] <= np.log10(5):
        inds_3.append(i)
    elif vmags[i] > bin_4[0] and vmags[i] <= bin_4[1] and lograds[i,0] <= np.log10(5):
        inds_4.append(i)
    elif vmags[i] > bin_5[0] and vmags[i] <= bin_5[1] and lograds[i,0] <= np.log10(5):
        inds_5.append(i)
    elif vmags[i] > bin_6[0] and vmags[i] <= bin_6[1] and lograds[i,0] <= np.log10(5):
        inds_6.append(i)
        
inds_1 = np.array(inds_1)
inds_2 = np.array(inds_2)
inds_3 = np.array(inds_3)
inds_4 = np.array(inds_4)
inds_5 = np.array(inds_5)
inds_6 = np.array(inds_6)
    
# store radii and densities for each bin
rads_1 = 10**lograds[inds_1]
dens_1 = 10**logdens[inds_1]
rads_2 = 10**lograds[inds_2]
dens_2 = 10**logdens[inds_2]
rads_3 = 10**lograds[inds_3]
dens_3 = 10**logdens[inds_3]
rads_4 = 10**lograds[inds_4]
dens_4 = 10**logdens[inds_4]
rads_5 = 10**lograds[inds_5]
dens_5 = 10**logdens[inds_5]
rads_6 = 10**lograds[inds_6]
dens_6 = 10**logdens[inds_6]

# let's find the min and max radius for each bin
min_rad_1 = np.min(rads_1)
max_rad_1 = np.max(rads_1)
min_rad_2 = np.min(rads_2)
max_rad_2 = np.max(rads_2)
min_rad_3 = np.min(rads_3)
max_rad_3 = np.max(rads_3)
min_rad_4 = np.min(rads_4)
max_rad_4 = np.max(rads_4)
min_rad_5 = np.min(rads_5)
max_rad_5 = np.max(rads_5)
min_rad_6 = np.min(rads_6)
max_rad_6 = np.max(rads_6)

# let's find the min and max COMMON radius for each bin
min_c_rad_1 = np.max(rads_1[:,0])
max_c_rad_1 = np.min(rads_1[:,-1])
min_c_rad_2 = np.max(rads_2[:,0])
max_c_rad_2 = np.min(rads_2[:,-1])
min_c_rad_3 = np.max(rads_3[:,0])
max_c_rad_3 = np.min(rads_3[:,-1])
min_c_rad_4 = np.max(rads_4[:,0])
max_c_rad_4 = np.min(rads_4[:,-1])
min_c_rad_5 = np.max(rads_5[:,0])
max_c_rad_5 = np.min(rads_5[:,-1])
min_c_rad_6 = np.max(rads_6[:,0])
max_c_rad_6 = np.min(rads_6[:,-1])

# let's store the individual minimum radii from each bin
indiv_min_1 = rads_1[:,0]
indiv_min_2 = rads_2[:,0]
indiv_min_3 = rads_3[:,0]
indiv_min_4 = rads_4[:,0]
indiv_min_5 = rads_5[:,0]
indiv_min_6 = rads_6[:,0]

# =============================================================================
# Step 1: lets resample profiles within the common radii for each bin

num_rad = 100
c_rad_1 = np.linspace(min_c_rad_1,max_c_rad_1,num_rad)
dens_1_res = np.zeros((len(dens_1[:,0]),num_rad))
for i in range(len(dens_1_res[:,0])):
    dens_1_res[i,:] = np.interp(c_rad_1, rads_1[i,:], dens_1[i,:])

c_rad_2 = np.linspace(min_c_rad_2,max_c_rad_2,num_rad)
dens_2_res = np.zeros((len(dens_2[:,0]),num_rad))
for i in range(len(dens_2_res[:,0])):
    dens_2_res[i,:] = np.interp(c_rad_2, rads_2[i,:], dens_2[i,:])
    
c_rad_3 = np.linspace(min_c_rad_3,max_c_rad_3,num_rad)
dens_3_res = np.zeros((len(dens_3[:,0]),num_rad))
for i in range(len(dens_3_res[:,0])):
    dens_3_res[i,:] = np.interp(c_rad_3, rads_3[i,:], dens_3[i,:])
    
c_rad_4 = np.linspace(min_c_rad_4,max_c_rad_4,num_rad)
dens_4_res = np.zeros((len(dens_4[:,0]),num_rad))
for i in range(len(dens_4_res[:,0])):
    dens_4_res[i,:] = np.interp(c_rad_4, rads_4[i,:], dens_4[i,:])
    
c_rad_5 = np.linspace(min_c_rad_5,max_c_rad_5,num_rad)
dens_5_res = np.zeros((len(dens_5[:,0]),num_rad))
for i in range(len(dens_5_res[:,0])):
    dens_5_res[i,:] = np.interp(c_rad_5, rads_5[i,:], dens_5[i,:])
    
c_rad_6 = np.linspace(min_c_rad_6,max_c_rad_6,num_rad)
dens_6_res = np.zeros((len(dens_6[:,0]),num_rad))
for i in range(len(dens_6_res[:,0])):
    dens_6_res[i,:] = np.interp(c_rad_6, rads_6[i,:], dens_6[i,:])
    
# =============================================================================
# =============================================================================


# =============================================================================
# Step 2: let's bin together the radii within the common radii of each bin

density_1 = np.zeros((num_rad))
for i in range(len(density_1)):
    density_1[i] = np.mean(dens_1_res[:,i])
  
density_2 = np.zeros((num_rad))
for i in range(len(density_2)):
    density_2[i] = np.mean(dens_2_res[:,i])

density_3 = np.zeros((num_rad))
for i in range(len(density_3)):
    density_3[i] = np.mean(dens_3_res[:,i])

density_4 = np.zeros((num_rad))
for i in range(len(density_4)):
    density_4[i] = np.mean(dens_4_res[:,i])

density_5 = np.zeros((num_rad))
for i in range(len(density_5)):
    density_5[i] = np.mean(dens_5_res[:,i])
   
density_6 = np.zeros((num_rad))
for i in range(len(density_6)):
    density_6[i] = np.mean(dens_6_res[:,i]) 
   
# =============================================================================
# =============================================================================


# =============================================================================
# Step 3: let's resample the profiles inside the least common radius

num_rad_inner = 30
c_rad_1_inner = np.linspace(np.min(indiv_min_1),min_c_rad_1,num_rad_inner)
dens_1_res_inner = np.zeros((len(dens_1[:,0]),num_rad_inner))
for i in range(len(dens_1_res_inner[:,0])):
    dens_1_res_inner[i,:] = np.interp(c_rad_1_inner, rads_1[i,:], dens_1[i,:],left=np.nan)

c_rad_2_inner = np.linspace(np.min(indiv_min_2),min_c_rad_2,num_rad_inner)
dens_2_res_inner = np.zeros((len(dens_2[:,0]),num_rad_inner))
for i in range(len(dens_2_res_inner[:,0])):
    dens_2_res_inner[i,:] = np.interp(c_rad_2_inner, rads_2[i,:], dens_2[i,:],left=np.nan)

c_rad_3_inner = np.linspace(np.min(indiv_min_3),min_c_rad_3,num_rad_inner)
dens_3_res_inner = np.zeros((len(dens_3[:,0]),num_rad_inner))
for i in range(len(dens_3_res_inner[:,0])):
    dens_3_res_inner[i,:] = np.interp(c_rad_3_inner, rads_3[i,:], dens_3[i,:],left=np.nan)

c_rad_4_inner = np.linspace(np.min(indiv_min_4),min_c_rad_4,num_rad_inner)
dens_4_res_inner = np.zeros((len(dens_4[:,0]),num_rad_inner))
for i in range(len(dens_4_res_inner[:,0])):
    dens_4_res_inner[i,:] = np.interp(c_rad_4_inner, rads_4[i,:], dens_4[i,:],left=np.nan)

c_rad_5_inner = np.linspace(np.min(indiv_min_5),min_c_rad_5,num_rad_inner)
dens_5_res_inner = np.zeros((len(dens_5[:,0]),num_rad_inner))
for i in range(len(dens_5_res_inner[:,0])):
    dens_5_res_inner[i,:] = np.interp(c_rad_5_inner, rads_5[i,:], dens_5[i,:],left=np.nan)

c_rad_6_inner = np.linspace(np.min(indiv_min_6),min_c_rad_6,num_rad_inner)
dens_6_res_inner = np.zeros((len(dens_6[:,0]),num_rad_inner))
for i in range(len(dens_6_res_inner[:,0])):
    dens_6_res_inner[i,:] = np.interp(c_rad_6_inner, rads_6[i,:], dens_6[i,:],left=np.nan)

# =============================================================================
# =============================================================================


# =============================================================================
# Step 4: let's bin together profiles inward of the minimum common radius based
# on their original extent

density_1_inner = np.zeros((num_rad_inner))
for i in range(len(density_1_inner)):
    density_1_inner[i] = np.nanmean(dens_1_res_inner[:,i])

density_2_inner = np.zeros((num_rad_inner))
for i in range(len(density_2_inner)):
    density_2_inner[i] = np.nanmean(dens_2_res_inner[:,i])

density_3_inner = np.zeros((num_rad_inner))
for i in range(len(density_3_inner)):
    density_3_inner[i] = np.nanmean(dens_3_res_inner[:,i])

density_4_inner = np.zeros((num_rad_inner))
for i in range(len(density_4_inner)):
    density_4_inner[i] = np.nanmean(dens_4_res_inner[:,i])

density_5_inner = np.zeros((num_rad_inner))
for i in range(len(density_5_inner)):
    density_5_inner[i] = np.nanmean(dens_5_res_inner[:,i])

density_6_inner = np.zeros((num_rad_inner))
for i in range(len(density_6_inner)):
    density_6_inner[i] = np.nanmean(dens_6_res_inner[:,i])

# =============================================================================
# =============================================================================


# =============================================================================
# Now let's plot the binned profiles

#colors = ['b','g','r','c','m']
lab = "{lim_1:.1f} < Mv < {lim_2:.1f}"
plt.figure(dpi=500)

# plotting common radii
plt.plot(np.log10(c_rad_1), np.log10(density_1), color=colors[0], 
         label=lab.format(lim_1=bin_1[0],lim_2=bin_1[1]))
plt.plot(np.log10(c_rad_2), np.log10(density_2), color=colors[1], 
         label=lab.format(lim_1=bin_2[0],lim_2=bin_2[1]))
plt.plot(np.log10(c_rad_3), np.log10(density_3), color=colors[2], 
         label=lab.format(lim_1=bin_3[0],lim_2=bin_3[1]))
plt.plot(np.log10(c_rad_4), np.log10(density_4), color=colors[3], 
         label=lab.format(lim_1=bin_4[0],lim_2=bin_4[1]))
plt.plot(np.log10(c_rad_5), np.log10(density_5), color=colors[4], 
         label=lab.format(lim_1=bin_5[0],lim_2=bin_5[1]))
plt.plot(np.log10(c_rad_6), np.log10(density_6), color=colors[5], 
         label=lab.format(lim_1=bin_6[0],lim_2=bin_6[1]))

# plotting inner radii
plt.plot(np.log10(c_rad_1_inner), np.log10(density_1_inner), linestyle='--', 
         color=colors[0])
plt.plot(np.log10(c_rad_2_inner), np.log10(density_2_inner), linestyle='--', 
         color=colors[1])
plt.plot(np.log10(c_rad_3_inner), np.log10(density_3_inner), linestyle='--', 
         color=colors[2])
plt.plot(np.log10(c_rad_4_inner), np.log10(density_4_inner), linestyle='--', 
         color=colors[3])
plt.plot(np.log10(c_rad_5_inner), np.log10(density_5_inner), linestyle='--', 
         color=colors[4])
plt.plot(np.log10(c_rad_6_inner), np.log10(density_6_inner), linestyle='--', 
         color=colors[5])

plt.xlabel('log(Radius [pc])')
plt.ylabel('log(Density [M$_\odot$/pc$^3$])')
plt.legend()

# =============================================================================
# =============================================================================


# =============================================================================
# let's try to plot the the binned profiles without binning the innermost data

# construct line collections for profiles inward of least common radii
ls_1 = []
for i in range(len(dens_1_res_inner[:,0])):
    a = np.zeros((len(c_rad_1_inner),2))
    a[:,0] = np.log10(c_rad_1_inner)
    a[:,1] = np.log10(dens_1_res_inner[i,:])
    ls_1.append(a)
ls_2 = []
for i in range(len(dens_2_res_inner[:,0])):
    a = np.zeros((len(c_rad_2_inner),2))
    a[:,0] = np.log10(c_rad_2_inner)
    a[:,1] = np.log10(dens_2_res_inner[i,:])
    ls_2.append(a)
ls_3 = []
for i in range(len(dens_3_res_inner[:,0])):
    a = np.zeros((len(c_rad_3_inner),2))
    a[:,0] = np.log10(c_rad_3_inner)
    a[:,1] = np.log10(dens_3_res_inner[i,:])
    ls_3.append(a)
ls_4 = []
for i in range(len(dens_4_res_inner[:,0])):
    a = np.zeros((len(c_rad_4_inner),2))
    a[:,0] = np.log10(c_rad_4_inner)
    a[:,1] = np.log10(dens_4_res_inner[i,:])
    ls_4.append(a)
ls_5 = []
for i in range(len(dens_5_res_inner[:,0])):
    a = np.zeros((len(c_rad_5_inner),2))
    a[:,0] = np.log10(c_rad_5_inner)
    a[:,1] = np.log10(dens_5_res_inner[i,:])
    ls_5.append(a)
ls_6 = []
for i in range(len(dens_6_res_inner[:,0])):
    a = np.zeros((len(c_rad_6_inner),2))
    a[:,0] = np.log10(c_rad_6_inner)
    a[:,1] = np.log10(dens_6_res_inner[i,:])
    ls_6.append(a)
    
fig, ax = plt.subplots(dpi=500)
ins_1 = LineCollection(ls_1, color=colors[0], linestyle='--', alpha=0.7)
ins_2 = LineCollection(ls_2, color=colors[1], linestyle='--', alpha=0.7)
ins_3 = LineCollection(ls_3, color=colors[2], linestyle='--', alpha=0.7)
ins_4 = LineCollection(ls_4, color=colors[3], linestyle='--', alpha=0.7)
ins_5 = LineCollection(ls_5, color=colors[4], linestyle='--', alpha=0.7)
ins_6 = LineCollection(ls_6, color=colors[5], linestyle='--', alpha=0.7)
ax.add_collection(ins_1)
ax.add_collection(ins_2)
ax.add_collection(ins_3)
ax.add_collection(ins_4)
ax.add_collection(ins_5)
ax.add_collection(ins_6)

# plotting common radii
ax.plot(np.log10(c_rad_1), np.log10(density_1), color=colors[0], 
         label=lab.format(lim_1=bin_1[0],lim_2=bin_1[1]))
ax.plot(np.log10(c_rad_2), np.log10(density_2), color=colors[1], 
         label=lab.format(lim_1=bin_2[0],lim_2=bin_2[1]))
ax.plot(np.log10(c_rad_3), np.log10(density_3), color=colors[2], 
         label=lab.format(lim_1=bin_3[0],lim_2=bin_3[1]))
ax.plot(np.log10(c_rad_4), np.log10(density_4), color=colors[3], 
         label=lab.format(lim_1=bin_4[0],lim_2=bin_4[1]))
ax.plot(np.log10(c_rad_5), np.log10(density_5), color=colors[4], 
         label=lab.format(lim_1=bin_5[0],lim_2=bin_5[1]))
ax.plot(np.log10(c_rad_6), np.log10(density_6), color=colors[5], 
         label=lab.format(lim_1=bin_6[0],lim_2=bin_6[1]))

ax.set_xlabel('log(Radius [pc])')
ax.set_ylabel('log(Density [M$_{\odot}$/pc$^3$])')
plt.legend()
ax.autoscale()
ax.set_xlim(-0.2,1.5)
#ax.set_ylim(1,7)
plt.show()
plt.close()

# =============================================================================
# =============================================================================


