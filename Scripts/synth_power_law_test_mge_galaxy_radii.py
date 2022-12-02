#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 20:39:56 2022

@author: christian
"""

from mgefit.mge_fit_1d import mge_fit_1d
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import mge1d_util as u
from astropy.stats import median_absolute_deviation as mad
from scipy.optimize import curve_fit
from tqdm import tqdm

samp_ext = 'normal'
num_rad = 100

# define the inner and outer slope parameters for mge
ngauss = 20
inner_slope = 4
outer_slope = 1
slope_ext = '_inslope_{}_outslope_{}'.format(inner_slope,outer_slope)

#pl_slopes = np.array([-2,-1.5,-1.0,-0.5,-0.2])

pl_slopes = np.array([-2,-1.75,-1.5,-1.25,-1.0,-0.75,-0.6,-0.5,-0.2])

# read in our data for the radial sampling
names_c, slopes_c, cen_dens_c, vmags_c, stone_slopes_c, \
    stone_dens_at_5pc_c, NSC_comp_c, lograds_c, logdens_c, dists_c = u.get_our_data('1.5x_pixel_scale','_or_20pc',True)
    

# conversion factor between radians and arcsec
arcsec_to_radian = 4.8481e-6
# HST pixel scales below in pc (for use in slope calculation)
scales = (arcsec_to_radian*0.0456)*dists_c*10**6
phys_scale = 10
scale_mult = 2

real_slopes = np.zeros((len(pl_slopes),len(names_c)))
mge_slopes = np.zeros((len(pl_slopes),len(names_c)))    
dens_slopes = np.zeros((len(pl_slopes),len(names_c)))
rads_all = np.zeros((len(pl_slopes),len(names_c),num_rad))
SBs_all = np.zeros((len(pl_slopes),len(names_c),num_rad))
dens_all = np.zeros((len(pl_slopes),len(names_c),num_rad))
summed_gauss_all = np.zeros((len(pl_slopes),len(names_c),num_rad))
 
true_slope = np.zeros((len(pl_slopes),len(names_c)))
true_int = np.zeros((len(pl_slopes),len(names_c)))
def power_law(x,a,b):
    return a*x+b 

for k in tqdm(range(len(pl_slopes)), position=0, leave=True):
#for k in range(len(pl_slopes)):    
    for l in range(len(names_c)):
        
        rad = 10**lograds_c[l,:] # pc
        #x = np.linspace(0.02,50,20000)
        SBs = rad**(pl_slopes[k])# + 500
        
        # ensure the radii are logarithmically spaced 
        rad1 = np.geomspace(rad[0], rad[-1], num_rad)
        SBs1 = 10**np.interp(np.log10(rad1), np.log10(rad), np.log10(SBs))      
        rad = rad1
        SBs = SBs1
        
        rads_all[k,l,:] = rad
        SBs_all[k,l,:] = SBs
        
# =============================================================================
#         # plot the resampled power-law profiles
#         plt.figure(dpi=500)
#         plt.plot(np.log10(rad),np.log10(SBs),marker='.')
#         plt.xlabel('Radius')
#         plt.ylabel('Synthetic Surface Brightness')
#         plt.show()
#         
#         plt.figure(dpi=500)
#         plt.plot(rad,SBs,marker='.')
#         plt.xlabel('Radius')
#         plt.ylabel('Synthetic Surface Brightness')
#         plt.show()
#         pdb.set_trace()
# =============================================================================

        # perform MGE fit
        mge = mge_fit_1d(rad, SBs, ngauss=ngauss, plot=False, quiet=True, 
                         inner_slope=inner_slope,outer_slope=outer_slope)
    
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
                sig = mge.sol[1][i] # in pc
                L_tot = (height)*2*np.pi*(sig)**2 # total area under 2-d gaussian
                
                if j == 0:
                    tots.append(L_tot_i)
                    heights.append(height)
                    sigs.append(sig)
                
                # convert total luminosity to total mass using M/L ratio
                M_tot = L_tot#*M_L
                # compute desity contribution from one gaussian component
                dens_comp = (M_tot/((2*np.pi)**(3/2)*sig**3))*np.exp(-(rad[j]**2/(2*sig**2))) # M_sol/pc^3
                # add density from each gaussian to total density
                #print(dens_comp, radii[j])
                densities[j] += dens_comp
                #pdb.set_trace()
            
        #plt.plot(np.log10(rad),np.log10(densities))
        
        #true_slope[k] = u.get_density_and_slope_simple(np.log10(rad[39:]),np.log10(densities[39:]))[0]

    # attempt to get "true" density slope
        #x1 = np.log10(rad[39:59]) # ~ 1.5 - 2.0 [log(pc)]
        #y1 = np.log10(densities[39:59])
        x1 = np.log10(rad[21:40]) # ~ 1.0 - 2.5 [log(pc)]
        y1 = np.log10(densities[21:40])
        pars, cov = curve_fit(f=power_law, xdata=x1, ydata=y1)
        true_slope[k,l] = pars[0]
        true_int[k,l] = pars[1]
        #pdb.set_trace()
        dens_all[k,l,:] = densities
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
        
        summed_gauss_all[k,l,:] = summed_gauss
# =============================================================================
#         plt.figure(dpi=500)
#         plt.plot(rad, SBs, color='b', linestyle='', marker='o')
#         for i in range(len(mge.sol[0])):
#             plt.plot(rad,gausses[i])
#  
#         plt.plot(rad, summed_gauss, color='orange',)#, alpha=0.7)    
#         
#         plt.yscale('log')
#         plt.xscale('log')
#         plt.ylim(min(SBs),max(SBs)+500)
#         plt.ylabel('Surface Brightness [L$_\odot$/pc$^2$]')
#         plt.xlabel('Radius ["]')
# =============================================================================

        real_slope, blah, blah, blah = u.get_density_and_slope(np.log10(rad),np.log10(SBs), scales[l],phys_scale,scale_mult)
        mge_slope, blah, blah, blah = u.get_density_and_slope(np.log10(rad),np.log10(summed_gauss), scales[l],phys_scale,scale_mult)
        dens_slope, dens_at_5pc, blah, blah = u.get_density_and_slope(np.log10(rad),np.log10(densities), scales[l],phys_scale,scale_mult)
    
        real_slopes[k,l] = real_slope
        mge_slopes[k,l] = mge_slope
        dens_slopes[k,l] = dens_slope
# =============================================================================
#         plt.text(0.62, 0.8, 'Data Slope = {:.2f}'.format(real_slope), 
#                  transform=plt.gcf().transFigure, size=10)
#         plt.text(0.62, 0.75, 'MGE Slope = {:.2f}'.format(mge_slope), 
#                  transform=plt.gcf().transFigure, size=10)
#         
# =============================================================================
#%%

medians = np.zeros(len(pl_slopes))
mads = np.zeros(len(pl_slopes))
perc_16 = np.zeros(len(pl_slopes))
perc_84 = np.zeros(len(pl_slopes))
medians_dens = np.zeros(len(pl_slopes))
mads_dens = np.zeros(len(pl_slopes))
perc_16_dens = np.zeros(len(pl_slopes))
perc_84_dens = np.zeros(len(pl_slopes))     
# =============================================================================
# plot the slopes 
for i in range(len(pl_slopes)):
    
    plt.figure(dpi=500)
    plt.title('{:.2f} Power Law Slope (SB Slope)'.format(pl_slopes[i]))
    #plt.hist(real_slopes[i,:], bins=10, label='True')
    a = plt.hist(mge_slopes[i,:], bins=10, label='MGE')
    plt.plot([pl_slopes[i],pl_slopes[i]],[0,np.max(a[0])],linestyle='--',color='r', label='True')
    plt.ylabel('# of Galaxies')
    plt.xlabel('Power-law Slope')
    med = np.median(mge_slopes[i,:])
    madev = mad(mge_slopes[i,:])
    per_16 = np.percentile(mge_slopes[i,:],16)
    per_84 = np.percentile(mge_slopes[i,:],84)
    medians[i] = med
    mads[i] = madev
    perc_16[i] = per_16
    perc_84[i] = per_84
    
    plt.text(0.7, 0.7, 'median = {:.4f}'.format(med),
             transform=plt.gcf().transFigure, size=10)
    plt.text(0.7, 0.67, 'mad = {:.4f}'.format(madev),
             transform=plt.gcf().transFigure, size=10)
    plt.legend()
    #plt.show()
    
    #plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Synth_Power_Law_Test/Realistic_Radii/slope_{:.2f}_dist_'.format(pl_slopes[i])+'{}x_pixel_scale'.format(scale_mult)+'_or_{}pc_'.format(phys_scale)+samp_ext+'.png',
    #        bbox_inches='tight', pad_inches=0.1, dpi=500)
    #pdb.set_trace()

    plt.figure(dpi=500)
    plt.title('{:.2f} Power Law Slope (Density Slope)'.format(pl_slopes[i]))
    #plt.hist(real_slopes[i,:], bins=10, label='True')
    a = plt.hist(dens_slopes[i,:], bins=10, label='MGE')
    #plt.plot([pl_slopes[i],pl_slopes[i]],[0,np.max(a[0])],linestyle='--',color='r', label='True')
    plt.ylabel('# of Galaxies')
    plt.xlabel('Power-law Slope')
    med = np.median(dens_slopes[i,:])
    madev = mad(dens_slopes[i,:])
    per_16 = np.percentile(dens_slopes[i,:],16)
    per_84 = np.percentile(dens_slopes[i,:],84)
    medians_dens[i] = med
    mads_dens[i] = madev
    perc_16_dens[i] = per_16
    perc_84_dens[i] = per_84
    
    plt.text(0.7, 0.7, 'median = {:.4f}'.format(med),
             transform=plt.gcf().transFigure, size=10)
    plt.text(0.7, 0.67, 'mad = {:.4f}'.format(madev),
             transform=plt.gcf().transFigure, size=10)
    plt.legend()
    #plt.show()
    
    #plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Synth_Power_Law_Test/Realistic_Radii/slope_{:.2f}_dist_'.format(pl_slopes[i])+samp_ext+'.png',
    #        bbox_inches='tight', pad_inches=0.1, dpi=500)
    #pdb.set_trace()


# =============================================================================
#%%

plt.figure(dpi=500)
plt.plot(pl_slopes, medians-pl_slopes, linestyle='', marker='.',color='b',label='SB profile')
plt.errorbar(pl_slopes, medians-pl_slopes,mads,linestyle='',capsize=4,color='b')
#plt.plot(pl_slopes, medians_dens-(pl_slopes-1), linestyle='', marker='.',color='c',label='Dens slope median')
plt.plot([np.min(pl_slopes),np.max(pl_slopes)],[0,0],color='r',linestyle='--')

#plt.plot(pl_slopes, mads, linestyle='', marker='.',color='r',label='SB slope mad')
#plt.plot(pl_slopes, mads_dens, linestyle='', marker='.',color='m',label='Dens slope mad')
#plt.title('SB Profile')
#plt.xlabel('Synthetic Power-Law Slope')
#plt.ylabel('Median Slope - True Slope ')
plt.xlabel('2D Synthetic Power-Law Slope')
plt.ylabel('Median Slope - True Slope')
#plt.legend()

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Synth_Power_Law_Test/Realistic_Radii/median_slopes_w_mads_just_2D_2x_pixel_scale_or_{}pc_'.format(phys_scale)+samp_ext+slope_ext+'.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# =============================================================================
# true_slopes = np.zeros_like(pl_slopes)
# for i in range(len(pl_slopes)):
#     true_slopes[i] = np.median(true_slope[i,:])
# =============================================================================
true_slopes = pl_slopes-1

plt.figure(dpi=500)
plt.plot(pl_slopes-1, medians_dens-true_slopes, linestyle='', marker='.',color='m',label='Density profile')
plt.errorbar(pl_slopes-1, medians_dens-true_slopes, mads_dens,linestyle='',capsize=4,color='m')
#plt.plot(pl_slopes, medians_dens-(pl_slopes-1), linestyle='', marker='.',color='c',label='Dens slope median')

#plt.plot(pl_slopes, mads, linestyle='', marker='.',color='r',label='SB slope mad')
#plt.plot(pl_slopes, mads_dens, linestyle='', marker='.',color='m',label='Dens slope mad')
#plt.title('3D Luminosity Density Profile')
plt.xlabel('3D Synthetic Power-Law Slope')
plt.ylabel('Median Slope - True Slope')
#plt.legend()

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Synth_Power_Law_Test/Realistic_Radii/median_slopes_w_mads_just_3D_2x_pixel_scale_or_{}pc_'.format(phys_scale)+samp_ext+slope_ext+'.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


#%%

# plot all density profiles
for i in range(len(pl_slopes)):
    plt.figure(dpi=500)
    plt.xlabel('log(Radius [pc])')
    plt.ylabel('3D Luminosity Density')
    plt.title('{:.2f} Power-Law'.format(pl_slopes[i]))
    for j in range(len(names_c)):
        plt.plot(np.log10(rads_all[i,j,:]), np.log10(dens_all[i,j,:]), color='b',alpha=0.5)
        plt.plot(np.log10(rads_all[i,j,:]), 
                 np.log10(rads_all[i,j,:])*(np.median(true_slope[i,:]))+np.median(true_int[i,:]),
                 color='r',alpha=0.5)
    plt.xlim(-0.3,1.5)
    plt.ylim(-6,1)
    plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Synth_Power_Law_Test/Realistic_Radii/dens_profiles_{:.2f}_pl_slope_'.format(pl_slopes[i])+'{}x_pixel_scale'.format(scale_mult)+'_or_{}pc_'.format(phys_scale)+samp_ext+slope_ext+'.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
            
#%%        
# plot all SB profiles
for i in range(len(pl_slopes)):
    plt.figure(dpi=500)
    plt.xlabel('log(Radius [pc])')
    plt.ylabel('Synthetic Surface Brightness')
    plt.title('{:.2f} Power-Law'.format(pl_slopes[i]))
    for j in range(len(names_c)):
        plt.plot(np.log10(rads_all[i,j,:]), np.log10(SBs_all[i,j,:]), color='b',alpha=0.5,label='Synthetic Data')
        plt.plot(np.log10(rads_all[i,j,:]), np.log10(summed_gauss_all[i,j,:]), color='m',alpha=0.5,label='MGE fit')
        if j == 0:
            plt.legend()
    plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Synth_Power_Law_Test/Realistic_Radii/SB_and_MGE_profiles_{:.2f}_pl_slope_'.format(pl_slopes[i])+'{}x_pixel_scale'.format(scale_mult)+'_or_{}pc_'.format(phys_scale)+samp_ext+slope_ext+'.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
#%%
# =============================================================================


# =============================================================================
# plot the distributions of 3D slopes for -2 & -3 (so -1, -2 2D SB)
# =============================================================================

plt.figure(dpi=500)
w = 0.001
#bs = np.arange(-.13,.1+2*w,w)
bs = np.arange(-.05,.1+2*w,w)
#plt.hist(dens_slopes[0,:]+3,color='b',bins=bs,alpha=0.6,label='-3.0')
#plt.hist(dens_slopes[1,:]+2.5,color='g',bins=bs,alpha=0.6,label='-2.5')
#plt.hist(dens_slopes[2,:]+2,color='r',bins=bs,alpha=0.6,label='-2.0')
#plt.hist(dens_slopes[3,:]+1.5,color='c',bins=bs,alpha=0.6,label='-1.5')
#plt.hist(dens_slopes[4,:]+1.2,color='y',bins=bs,alpha=0.6,label='-1.2')

plt.hist(dens_slopes[4,:]+1.2,color='y',bins=bs,alpha=0.6,label='-1.2')
plt.hist(dens_slopes[3,:]+1.5,color='c',bins=bs,alpha=0.6,label='-1.5')
plt.hist(dens_slopes[2,:]+2,color='r',bins=bs,alpha=0.6,label='-2.0')
plt.hist(dens_slopes[1,:]+2.5,color='g',bins=bs,alpha=0.6,label='-2.5')
plt.hist(dens_slopes[0,:]+3,color='b',bins=bs,alpha=0.6,label='-3.0')
#plt.hist(dens_slopes[5,:]+1,color='c',bins=bs,alpha=0.6,label='-1.0')
plt.text(0.17, 0.84, 'inner_slope = {}'.format(inner_slope),
             transform=plt.gcf().transFigure, size=10)
plt.text(0.17, 0.8, 'outer_slope = {}'.format(outer_slope),
             transform=plt.gcf().transFigure, size=10)
plt.text(0.17, 0.76, 'ngauss = {}'.format(ngauss),
             transform=plt.gcf().transFigure, size=10)

plt.xlabel('Measured - True 3D Density Slope')
#plt.xlim(-.13,.13)
plt.xlim(-0.02,0.02)
plt.legend()

plt.savefig('/Users/christian/OneDrive/Desktop/TDE Code/Plots/Synth_Power_Law_Test/Realistic_Radii/All_3D_slope_dists_'+'{}x_pixel_scale'.format(scale_mult)+'_or_{}pc_'.format(phys_scale)+samp_ext+slope_ext+'.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)
# =============================================================================
# samp = [0.5,1.0,2]
# sigs_05 = [.0036,.0124,.0097]
# sigs_10 = [.0036,.0037,.0030]
# sigs_15 = [.0116,.0041,.0139]
# sigs_20 = [.0319,.0145,.0131]
# plt.figure(dpi=500)
# plt.plot(samp,sigs_05,linestyle='',marker='.',label='-0.5',alpha=0.6)
# plt.plot(samp,sigs_10,linestyle='',marker='.',label='-1.0',alpha=0.6)
# plt.plot(samp,sigs_15,linestyle='',marker='.',label='-1.5',alpha=0.6)
# plt.plot(samp,sigs_20,linestyle='',marker='.',label='-2.0',alpha=0.6)
# plt.xlabel("Sampling (compared to our original)")
# plt.ylabel("$\sigma$ of Slopes")
# plt.legend()
# 
# =============================================================================
# =============================================================================
#%%

meds_dens_41 = np.array([-2.99806009, -2.49897555, -1.99944993, -1.49941478, -1.19617891])
mads_dens_41 = np.array([0.00084932, 0.00044831, 0.00020317, 0.00016298, 0.00095972])
meds_dens_31 = np.array([-2.9919602 , -2.49558816, -1.9975731 , -1.49855155, -1.19589061])
mads_dens_31 = np.array([0.00183997, 0.00066387, 0.00021866, 0.00041747, 0.00092318])
meds_dens_21 = np.array([-2.97093293, -2.48444385, -1.99181388, -1.49542487, -1.19324988])
mads_dens_21 = np.array([0.00940537, 0.00478209, 0.00248083, 0.00119652, 0.00129874])
meds_dens_11 = np.array([-2.82178603, -2.38996315, -1.93650322, -1.46634728, -1.17648845])
mads_dens_11 = np.array([0.02470376, 0.0086878 , 0.00335665, 0.00180582, 0.00206359])

meds_dens_23 = np.array([-2.9711997 , -2.4835316 , -1.99162827, -1.49518651, -1.1954123 ])
mads_dens_23 = np.array([0.00914894, 0.00526386, 0.0026382 , 0.00142316, 0.00085191])

colors = ['r',(.8,.5,.0),'y',(.2,.8,0.),'c',(.0,.3,.7),(.2,.0,.8),(.7,0,.7)]

plt.figure(dpi=500)
plt.plot(true_slopes, meds_dens_41-true_slopes, linestyle='', marker='.',color=colors[0],label='inner = 4, outer = 1')
plt.errorbar(true_slopes, meds_dens_41-true_slopes, mads_dens_41,linestyle='',capsize=4,color=colors[0])

plt.plot(true_slopes, meds_dens_31-true_slopes, linestyle='', marker='.',color=colors[1],label='inner = 3, outer = 1')
plt.errorbar(true_slopes, meds_dens_31-true_slopes, mads_dens_31,linestyle='',capsize=4,color=colors[1])

#plt.plot(true_slopes, meds_dens_21-true_slopes, linestyle='', marker='.',color=colors[2],label='inner = 2, outer = 1')
#plt.errorbar(true_slopes, meds_dens_21-true_slopes, mads_dens_21,linestyle='',capsize=4,color=colors[2])

#plt.plot(true_slopes, meds_dens_11-true_slopes, linestyle='', marker='.',color=colors[3],label='inner = 1, outer = 1')
#plt.errorbar(true_slopes, meds_dens_11-true_slopes, mads_dens_11,linestyle='',capsize=4,color=colors[3])

plt.plot(true_slopes, meds_dens_23-true_slopes, linestyle='', marker='.',color=colors[4],label='inner = 2, outer = 3; Default')
plt.errorbar(true_slopes, meds_dens_23-true_slopes, mads_dens_23,linestyle='',capsize=4,color=colors[4])

plt.plot(true_slopes, np.zeros_like(true_slopes),linestyle='--',color='k')

plt.xlabel('3D Synthetic Power-Law Slope')
plt.ylabel('Median Slope - True Slope')

plt.ylim(-0.0005,0.0105)

plt.legend()


meds_41 = np.array([-1.9999747 , -1.49999514, -0.99999918, -0.50000108, -0.19999644])
mads_41 = np.array([2.15712048e-05, 8.55393929e-06, 5.77342555e-06, 5.71551916e-06,
       3.88635634e-05])
meds_31 = np.array([-1.99978871, -1.49993878, -0.99998012, -0.49999028, -0.20000222])
mads_31 = np.array([1.17586896e-04, 4.42990356e-05, 1.87096985e-05, 2.81115753e-05,
       1.57887758e-05])
meds_21 = np.array([-1.99802997, -1.4994753 , -0.99990807, -0.50002302, -0.20001139])
mads_21 = np.array([4.77383699e-04, 2.42061058e-04, 4.00177862e-05, 5.30230531e-05,
       4.40022924e-05])
meds_11 = np.array([-1.95568365, -1.4838768 , -0.99498679, -0.49880515, -0.19971803])
mads_11 = np.array([0.00066075, 0.00199205, 0.0013984 , 0.0005113 , 0.00015408])

meds_23 = np.array([-1.99811873, -1.49942098, -0.99993963, -0.50002295, -0.20001026])
mads_23 = np.array([4.79834459e-04, 2.33987328e-04, 1.20847474e-04, 9.88655345e-05,
       7.95748366e-05])

colors = ['r',(.8,.5,.0),'y',(.2,.8,0.),'c',(.0,.3,.7),(.2,.0,.8),(.7,0,.7)]

plt.figure(dpi=500)
plt.plot(pl_slopes, meds_41-pl_slopes, linestyle='', marker='.',color=colors[0],label='inner = 4, outer = 1')
plt.errorbar(pl_slopes, meds_41-pl_slopes, mads_41,linestyle='',capsize=4,color=colors[0])

plt.plot(pl_slopes, meds_31-pl_slopes, linestyle='', marker='.',color=colors[1],label='inner = 3, outer = 1')
plt.errorbar(pl_slopes, meds_31-pl_slopes, mads_31,linestyle='',capsize=4,color=colors[1])

#plt.plot(pl_slopes, meds_21-pl_slopes, linestyle='', marker='.',color=colors[2],label='inner = 2, outer = 1')
#plt.errorbar(pl_slopes, meds_21-pl_slopes, mads_21,linestyle='',capsize=4,color=colors[2])

#plt.plot(pl_slopes, meds_11-pl_slopes, linestyle='', marker='.',color=colors[3],label='inner = 1, outer = 1')
#plt.errorbar(pl_slopes, meds_11-pl_slopes, mads_11,linestyle='',capsize=4,color=colors[3])

plt.plot(pl_slopes, meds_23-pl_slopes, linestyle='', marker='.',color=colors[4],label='inner = 2, outer = 3; Default')
plt.errorbar(pl_slopes, meds_23-pl_slopes, mads_23,linestyle='',capsize=4,color=colors[4])

plt.plot(pl_slopes, np.zeros_like(pl_slopes),linestyle='--',color='k')

plt.xlabel('2D Synthetic Power-Law Slope')
plt.ylabel('Median Slope - True Slope')

plt.ylim(-0.00005,0.0005)

plt.legend()

#%%

meds_dens_41_20 = np.array([-2.99806009, -2.49897555, -1.99944993, -1.49941478, -1.19617891])
mads_dens_41_20 = np.array([0.00084932, 0.00044831, 0.00020317, 0.00016298, 0.00095972])
meds_dens_41_30 = np.array([-2.99795618, -2.49887746, -1.99938218, -1.49941984, -1.19648132])
mads_dens_41_30 = np.array([0.00074696, 0.00039869, 0.00015178, 0.00039177, 0.00051461])
meds_dens_41_40 = np.array([-2.99801929, -2.49888604, -1.99925443, -1.4993547 , -1.19644168])
mads_dens_41_40 = np.array([0.00065413, 0.00038197, 0.0002031 , 0.00014994, 0.00075191])
meds_dens_41_50 = np.array([-2.99796548, -2.49890598, -1.99941915, -1.49939334, -1.19641766])
mads_dens_41_50 = np.array([6.90028838e-04, 3.55309535e-04, 1.69552266e-04, 8.84509549e-05,
       7.13392072e-04])

colors = ['r',(.8,.5,.0),'y',(.2,.8,0.),'c',(.0,.3,.7),(.2,.0,.8),(.7,0,.7)]

plt.figure(dpi=500)
plt.title('inner_slope = 4, outer_slope = 1')
plt.plot(true_slopes, meds_dens_41_20-true_slopes, linestyle='', marker='.',color=colors[0],label='ngauss = 20')
plt.errorbar(true_slopes, meds_dens_41_20-true_slopes, mads_dens_41_20,linestyle='',capsize=4,color=colors[0])

plt.plot(true_slopes, meds_dens_41_30-true_slopes, linestyle='', marker='.',color=colors[1],label='ngauss = 30')
plt.errorbar(true_slopes, meds_dens_41_30-true_slopes, mads_dens_41_30,linestyle='',capsize=4,color=colors[1])

plt.plot(true_slopes, meds_dens_41_40-true_slopes, linestyle='', marker='.',color=colors[2],label='ngauss = 40')
plt.errorbar(true_slopes, meds_dens_41_40-true_slopes, mads_dens_41_40,linestyle='',capsize=4,color=colors[2])

plt.plot(true_slopes, meds_dens_41_50-true_slopes, linestyle='', marker='.',color=colors[3],label='ngauss = 50')
plt.errorbar(true_slopes, meds_dens_41_50-true_slopes, mads_dens_41_50,linestyle='',capsize=4,color=colors[3])

plt.plot(true_slopes, np.zeros_like(true_slopes),linestyle='--',color='k')

plt.xlabel('3D Synthetic Power-Law Slope')
plt.ylabel('Median Slope - True Slope')

plt.ylim(-0.0005,0.006)

plt.legend()


plt.figure(dpi=500)
plt.plot(true_slopes, mads_dens_41_20, linestyle='',marker='.',color=colors[0],label='ngauss = 20',alpha=0.7)
plt.plot(true_slopes, mads_dens_41_30, linestyle='',marker='.',color=colors[1],label='ngauss = 30',alpha=0.7)
plt.plot(true_slopes, mads_dens_41_40, linestyle='',marker='.',color=colors[2],label='ngauss = 40',alpha=0.7)
plt.plot(true_slopes, mads_dens_41_50, linestyle='',marker='.',color=colors[3],label='ngauss = 50',alpha=0.7)

plt.ylabel('MAD of Median Slope - True Slope')
plt.xlabel('3D Synthetic Power-Law Slope')
plt.legend(loc='upper center')


#%%

mads_dens_41_20_extended = np.array([8.49321618e-04, 6.42463575e-04, 4.48306979e-04, 2.76315273e-04,
       2.03169639e-04, 1.47947929e-04, 9.26351140e-05, 1.62979737e-04,
       9.59715573e-04])


plt.figure(dpi=500)
plt.plot(true_slopes, mads_dens_41_20_extended, linestyle='', marker='.', color=colors[6])
#pars = np.polyfit(true_slopes, mads_dens_41_20, deg=3)
slope_array = np.arange(-3,-1.15,0.01)
#plt.plot(slope_array, pars[0]*slope_array**3+pars[1]*slope_array**2+pars[2]*slope_array+pars[3],linestyle='--',color='r')

f = interpolate.CubicSpline(true_slopes, mads_dens_41_20_extended)
y = f(slope_array)
plt.plot(slope_array, y,linestyle='--',color='r')

plt.text(0.17, 0.84, 'inner_slope = {}'.format(inner_slope),
             transform=plt.gcf().transFigure, size=10)
plt.text(0.17, 0.8, 'outer_slope = {}'.format(outer_slope),
             transform=plt.gcf().transFigure, size=10)
plt.text(0.17, 0.76, 'ngauss = {}'.format(ngauss),
             transform=plt.gcf().transFigure, size=10)

med_mad = np.mean(mads_dens_41_20_extended)
plt.plot(slope_array,np.ones_like(slope_array)*med_mad,color='c')

plt.ylabel('MAD of Median Slope - True Slope')
plt.xlabel('3D Synthetic Power-Law Slope')










