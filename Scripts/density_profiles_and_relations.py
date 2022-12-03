#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 10:20:02 2021

This code serves the following purposes:
    
    - Read in Pechetti+20 density and MGE data (and reconstruct density profiles)
    - Concatenate our data/results with those of Pechetti+20
    - Plot color-mass distribution of our sample
    - Compare types between David and sample sources
    - Compare results to Stone&Metzger (2016)
    - Fit for (and plot) the density relation by type using linmix
    - Plot the mass CDF and nucleation fraction of our sample

@author: christian
"""

import matplotlib.pyplot as plt
import numpy as np
import mge1d_util as u
from astropy.io import fits
from astropy.io import ascii
import csv
import pdb
from linmix import linmix
import multiprocessing
from matplotlib.colors import LogNorm
import random
import re
import warnings
warnings.filterwarnings("ignore")


#%%

multiprocessing.set_start_method("fork", force=True)

data = ascii.read('../Data_Sets/Pechetti_20_tables/Pechetti_Table_A3.txt', format='latex') 
data_1 = ascii.read('../Data_Sets/Pechetti_20_tables/Pechetti_Table_1.txt', format='latex') 
data_2 = ascii.read('../Data_Sets/Pechetti_20_tables/Pechetti_Table_A1.txt', format='latex') 

#%%

# Read in table ffor Pechetti+20 galaxies with RA and Dec
all_names = []
all_RAs = []
all_DECs = []
pechetti_filename = '../Data_Sets/Pechetti_20_tables/pechetti_20_ra_dec_NED.csv'
# fields: name	instrument	channel	filter	pa	l_pa	u_pa	ell	l_ell	u_ell	n	l_n	u_n	reff	l_reff	
#         u_reff	ieff	l_ieff	u_ieff	m	l_m	u_m	v	u_v	i	u_i	logm	u_logm
with open(pechetti_filename, 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        #pdb.set_trace()
        if row[0] != 'No.': # marks the start of data row
            all_names.append(re.split('>|<',row[1])[2])
            all_RAs.append(float(row[2]))
            all_DECs.append(float(row[3]))
       
all_names = np.array(all_names)
all_RAs = np.array(all_RAs)
all_DECs = np.array(all_DECs)

# Adjust some of the names to match the names in the Pechetti+20 paper
all_names[np.where(all_names == 'NGC 3115 DW01')] = 'NGC 3115B'
all_names[np.where(all_names == 'ESO 274- G 001')] = 'ESO 274-1'
all_names[np.where(all_names == 'MESSIER 051a')] = 'NGC 5194'
all_names[np.where(all_names == 'MESSIER 051b')] = 'NGC 5195'
all_names[np.where(all_names == 'MESSIER 101')] = 'NGC 5457'
all_names[np.where(all_names == 'MESSIER 063')] = 'NGC 5055'
all_names[np.where(all_names == 'MESSIER 083')] = 'NGC 5236'
all_names[np.where(all_names == 'Circinus Galaxy')] = 'Circinus'

# Remove the spaces in the names
for i in range(len(all_names)):
    s = all_names[i].split()
    if len(s) > 1:
        all_names[i] = s[0]+s[1]


# from Pechetti+20
names_r = ['Circinus','ESO274-1','IC5052','IC5332','NGC2784','NGC2787',
             'NGC2903','NGC3115','NGC3115B','NGC3184','NGC3274','NGC3344',
             'NGC3593','NGC4242','NGC4460','NGC4517','NGC4592','NGC4600',
             'NGC4605','NGC4941','NGC5055','NGC5068','NGC5194','NGC5195',
             'NGC5236','NGC5238','NGC5457','NGC6503','NGC7713']


# now let's construct RA and DEC arrays for the galaxies in the order they appear above
RAs = np.zeros(len(names_r))
DECs = np.zeros(len(names_r))
temp_names = np.zeros(len(names_r)).astype(str)
for i in range(len(names_r)):
    ind = np.where(all_names == names_r[i])[0]
    RAs[i] = all_RAs[ind][0]
    DECs[i] = all_DECs[ind][0]
    #temp_names[i] = all_names[ind][0]


# get galaxy masses and colors from david's table 
gal_mass = np.zeros(len(names_r))
gi_color = np.zeros(len(names_r))
vmags = np.zeros(len(names_r))
dists = np.zeros(len(names_r))
galtypes = np.zeros(len(names_r))
for i in range(len(names_r)):
    gal_info = u.get_David_gal_data(names_r[i], RAs[i], DECs[i])
    temp_names[i] = gal_info[0]
    gal_mass[i] = 10**gal_info[2]
    gi_color[i] = gal_info[4]
    vmags[i] = gal_info[3]
    dists[i] = gal_info[1]
    galtypes[i] = gal_info[5]


renuka_bmags_init = data_2['M$_B$'][1:]
renuka_ttypes_init = data_2['T-type'][1:]
renuka_bmags = np.array(renuka_bmags_init).astype(np.float)
renuka_ttypes = np.array(renuka_ttypes_init).astype(int)

#convert t-types to Early/Late (0/1)
morphs_r = np.zeros(len(renuka_ttypes))
for i in range(len(morphs_r)):
    if renuka_ttypes[i] <= 0:
        morphs_r[i] = 0
    else:
        morphs_r[i] = 1


#%%
# put in density/slope values manually as well as errors
log_dens_5pc_r = np.array([3.84,2.78,3.06,-99,-99,3.93,4.03,-99,3.21,-99,2.35,3.83,-99,
                         2.03,-99,2.61,2.75,3.22,3.65,3.75,3.60,2.68,-99,3.93,3.66,
                         2.04,3.47,3.42,1.83])
log_dens_5pc_r_err = np.array([.10,.14,.11,-99,-99,.11,.10,-99,.10,-99,.11,.11,-99,.10,
                             -99,.11,.12,.10,.12,.11,.12,.11,-99,.10,.12,.10,.10,.16,
                             .10])
gamma_r = -1*np.array([.92,2.16,2.06,-99,-99,1.42,1.56,-99,1.86,-99,2.63,1.43,-99,1.80,
                     -99,2.44,2.87,2.09,2.75,1.80,1.87,1.68,-99,1.23,2.48,2.76,1.82,
                     2.75,3.22])
gamma_r_err = np.array([.12,.07,.04,-99,-99,.03,.08,-99,.04,-99,.04,.11,-99,.12,-99,
                     .03,.02,.04,.02,.12,.06,.26,-99,.24,.04,.03,.02,.15,.01])



# modify arrays to only include gals with the appropriate measurments
a = np.where(log_dens_5pc_r != -99)
names_r = np.array(names_r)
names_r_corrected = names_r[a]
gi_color = gi_color[a]
dists = dists[a]
RAs = RAs[a]
DECs = DECs[a]
gal_mass = gal_mass[a]
bmags = renuka_bmags[a]
vmags = vmags[a]
galtypes = galtypes[a]
morphs_r = morphs_r[a]
slopes_r = gamma_r[a]
slopes_r_e = gamma_r_err[a]
cen_dens = log_dens_5pc_r[a]
cen_dens_e = log_dens_5pc_r_err[a]
slope_err = np.median(gamma_r_err[a])
cen_dens_err = np.median(log_dens_5pc_r_err[a])


# =============================================================================
# =================== Get Pechetti+20 Density Data ============================
# =============================================================================

# let's get the density profiles for PEchetti+20 galaxies using the MGE data
nsc_ml = np.array([5.48,2.10,4.53,5.48,-99,3.43,5.48,-99,1.36,-99,1.26,0.68,5.48,
                   1.05,-99,0.81,1.43,1.58,5.48,1.55,0.44,1.2,1.2,5.48,4.20,0.50,
                   1.80,5.48,0.32])
scales = np.array([.05,.04,.05,.04,.045,.05,.05,.04,.05,.045,.04,.04,.04,.045,
                   .05,.045,.08,.05,.045,.045])
outer_rad_arcsec = np.array([0.93138189,0.23741352,0.76973705,0.38235687,0.43015148,0.32644360,0.43449645,
                             0.42059258,0.47794609,0.22073570,0.50891705,0.13401371,0.17872269,0.47316665,
                             0.32644360,0.24809267,0.20137295,0.29379924,0.35674454,0.63614631,0.29379924,
                             1.63349997])

dens_20 = np.zeros((len(names_r),1000))
rad_20 = np.zeros((len(names_r),1000))
outer_pcs = np.zeros(len(names_r))

mge_data = np.loadtxt('../Data_Sets/Pechetti_20_tables/Pechetti_MGE_Table_2.txt', dtype=str)

for i in range(0,len(names_r)):
	
    pc = dists[i]*np.pi/0.648
    radius = scales[i]*pc
    outer_pc = outer_rad_arcsec[i]*pc
    outer_pcs[i] = outer_pc
	
    gal_inds = np.where(mge_data[:,0] == names_r[i])[0]
    filt = mge_data[gal_inds[0],1]
    I = 10**mge_data[gal_inds,3]
    
    # convert luminosities to I-band if not already F814W for the I-Band M/L's
    if filt == 'F814W':
        M = I*nsc_ml[i]
    if filt == 'F110W':
        I = 0
        M = I*nsc_ml[i]
    if filt == 'F160W':
        I = 0
        M = I*nsc_ml[i]
    if filt == 'F547M':
        I = 0
        M = I*nsc_ml[i]
    if filt == 'F555W':
        I = 0
        M = I*nsc_ml[i]
    if filt == 'F606W':
        I = 0
        M = I*nsc_ml[i]
    if filt == 'F625W':
        I = 0
        M = I*nsc_ml[i]
    
    
    sigma = 10**mge_data[gal_inds,4]*pc
    q_prime = mge_data[gal_inds,5]
    inc = 60*np.pi/180
    q = (q_prime**2 - (np.cos(inc)**2))/(np.sin(inc)**2)
	
    r = np.linspace(radius,100,1000)
    rho = np.zeros(len(r))
    for k in range(0,len(r)):
        for j in range(0,len(gal_inds)):
            rho[k]+= M[j]*np.exp(-(r[k]**2)/(2*sigma[j]**2))/(q[j]*(np.sqrt(2*np.pi)*sigma[j])**3)

    dens_20[i,:] = rho 
    rad_20[i,:] = r

# =============================================================================
# =============================================================================
# =============================================================================


#%%


# =============================================================================
# Read in our data
slope_ext = '2x_pixel_scale'
phys_ext = '_or_10pc_extinction_corr_nsa_ml_w_vi_color'
names_c, slopes_c, cen_dens_c, vmags_c, stone_slopes_c, \
    stone_dens_at_5pc_c, NSC_comp_c, lograds_c, logdens_c, dists_c, \
    all_MLs_c, all_ML_types_c, MLs_c, ML_types_c, dist_flags_c = u.get_our_data(slope_ext,phys_ext,True)

names_c = np.array(names_c, dtype='<U100')

# alter names to reflect David's catalog names and correct truncation issue
names_c[np.where(names_c == 'PGC 02888')[0]] = 'PGC 028887'
names_c[np.where(names_c == 'PGC 05039')[0]] = 'PGC 050395'
names_c[np.where(names_c == 'bts76')[0]] = 'SDSSJ115844.10+273506.0'
names_c[np.where(names_c == 'ddo084')[0]] = 'UGC05829'
names_c[np.where(names_c == 'kk2000-03')[0]] = 'PGC009140'
names_c[np.where(names_c == 'kk2000-53')[0]] = '[KK2000] 53'
names_c[np.where(names_c == 'kk96')[0]] = 'KK96'
names_c[np.where(names_c == 'leg09')[0]] = 'LeG09'
names_c[np.where(names_c == 'lvj1217+47')[0]] = 'LV J1217+4703'
names_c[np.where(names_c == 'ngc5011c')[0]] = 'ESO269-068'
names_c[np.where(names_c == 'pgc4310323')[0]] = 'SDSSJ120531.04+310434.1'
names_c[np.where(names_c == 'ugc07242')[0]] = 'UGC07242'


# get the log(masses) from David for the galaxies he has in his table
logmass_c = []
gi_color_c = []
vmags_c = []
RAs_c = []
Decs_c = []
galtype_c = []
for i in range(len(names_c)):
    n = names_c[i].split(' ')
    if len(n) == 1:
        logmass_c.append(u.get_galaxy_logmass(n[0]))
        RA, DEC = u.get_galaxy_RA_and_Dec(n[0])
        RAs_c.append(RA)
        Decs_c.append(DEC)
        gi_color_c.append(u.get_David_gal_data(n[0],RA,DEC)[4])
        vmags_c.append(u.get_David_gal_data(n[0],RA,DEC)[3])
        galtype_c.append(u.get_David_gal_data(n[0],RA,DEC)[5])
    elif len(names_c)-i < 10:
        logmass_c.append(u.get_galaxy_logmass(n[0]+' '+n[1]))
        RA, DEC = u.get_galaxy_RA_and_Dec(n[0]+' '+n[1])
        RAs_c.append(RA)
        Decs_c.append(DEC)
        gi_color_c.append(u.get_David_gal_data(n[0]+' '+n[1],RA,DEC)[4])
        vmags_c.append(u.get_David_gal_data(n[0]+' '+n[1],RA,DEC)[3])
        galtype_c.append(u.get_David_gal_data(n[0]+' '+n[1],RA,DEC)[5])
    else:
        logmass_c.append(u.get_galaxy_logmass(n[0]+n[1]))
        RA, DEC = u.get_galaxy_RA_and_Dec(n[0]+n[1])
        RAs_c.append(RA)
        Decs_c.append(DEC)
        gi_color_c.append(u.get_David_gal_data(n[0]+n[1],RA,DEC)[4])
        vmags_c.append(u.get_David_gal_data(n[0]+n[1],RA,DEC)[3])
        galtype_c.append(u.get_David_gal_data(n[0]+n[1],RA,DEC)[5])
        
galtype_c = np.array(galtype_c)
logmass_c = np.array(logmass_c)
gi_color_c = np.array(gi_color_c)
vmags_c = np.array(vmags_c)
RAs_c = np.array(RAs_c)
Decs_c = np.array(Decs_c)

# order -> ['bts76', 'ddo084', 'kk2000-03', 'kk2000-53', 'kk96', 'leg09',
# 'lvj1217+4703', 'ngc5011c', 'pgc4310323', 'ugc07242']
# below updated with David's sample_4 masses for 'SDSSJ115844.10+273506.0', 'PGC009140', 
# '[KK2000] 53','LV J1217+4703'
logmass_22 = np.array([7.53, 8.52871022, 8.16, 6.85, 7.07330028,
                    6.99333377, 6.58, 7.51866251, 6.61198218, 7.75])

# manual assignment of Hoyer+22 masses
logmass_c[len(logmass_c)-10:] = logmass_22

# get index list for galaxies that had logmasses
mass_inds = []
for i in range(len(logmass_c)):
    if logmass_c[i] != 0:
        mass_inds.append(i)


# un-comment below to single out galaxies with masses and known NSCs
#q = np.where(NSC_comp_c == 1)
#mass_inds = np.intersect1d(mass_inds,q)

# remove galaxies without masses
names_corrected = names_c[mass_inds]
logmass_corrected = logmass_c[mass_inds]
lograds_corrected = lograds_c[mass_inds]
logdens_corrected = logdens_c[mass_inds]
slopes_corrected = slopes_c[mass_inds]
cen_dens_corrected = cen_dens_c[mass_inds]
morphs_corrected = galtype_c[mass_inds]
NSC_comp_corrected = NSC_comp_c[mass_inds]
gi_color_corrected = gi_color_c[mass_inds]
dists_corrected = dists_c[mass_inds]
vmags_corrected = vmags_c[mass_inds]
galtype_corrected = galtype_c[mass_inds]
RAs_corrected = RAs_c[mass_inds]
Decs_corrected = Decs_c[mass_inds]


#%%


# =============================================================================
# ================ NUCLEATION FRACTION AND MASS CDF OF OUR SAMPLE =============
# =============================================================================

# let's compute the nucleation fraction of our galaxies as a function of galmass
all_logm_init = np.concatenate((np.log10(gal_mass),logmass_corrected))
NSC_comp_r = np.zeros(len(gal_mass))+1
all_NSC_comp_init = np.concatenate((NSC_comp_r, NSC_comp_corrected))

# get rid of galaxies where we do not know if nucleated
all_logm = all_logm_init[np.where(all_NSC_comp_init != -1)]
all_NSC_comp = all_NSC_comp_init[np.where(all_NSC_comp_init != -1)]

blah, bins, blah = plt.hist(all_logm, bins=20)

nuc_frac = np.zeros(len(bins)-1)
# separate the galaxies by bin
for i in range(len(bins)-1):
    nsc_comps = []
    for j in range(len(all_logm)):
        if all_logm[j] >= bins[i] and all_logm[j] < bins[i+1]:
            nsc_comps.append(all_NSC_comp[j])
    nsc = np.array(nsc_comps)
    if len(nsc) != 0:
        nuc_frac[i] = len(np.where(nsc==1)[0])/len(nsc)
    else:
        nuc_frac[i] = -1.

# compute the mean of the bins
bin_means = bins[0:-1]+(bins[1:]-bins[0:-1])/2
bin_widths = (bins[1:]-bins[0:-1])/2


nuc_frac_corrected = nuc_frac[np.where(nuc_frac != -1.)]
bin_means_corrected = bin_means[np.where(nuc_frac != -1.)]
bin_widths_corrected = bin_widths[np.where(nuc_frac != -1.)]

#%%

# let's try to compute the CDF of nucleation over galaxy mass
ms = np.arange(6.5,11.5,0.01)
mass_cdf = np.zeros(len(ms))
tot_gals = len(all_logm_init)
for i in range(len(ms)):
    num_gals = 0
    for j in range(len(all_logm_init)):
        if all_logm_init[j] <= ms[i]:
            num_gals += 1
    mass_cdf[i] = num_gals/tot_gals

plt.figure(dpi=500)
plt.plot(ms, mass_cdf)

#%%
#let's do the same as above but for each source
inds_20 = (0,22)
inds_05 = (inds_20[1], inds_20[1]+45)
inds_95 = (inds_05[1], inds_05[1]+6)
inds_17 = (inds_95[1], inds_95[1]+25)
inds_22 = (inds_17[1], inds_17[1]+10)


mass_cdf_20 = np.zeros(len(ms))
logm_20 = all_logm_init[inds_20[0]:inds_20[1]]
tot_gals_20 = len(logm_20)
for i in range(len(ms)):
    num_gals = 0
    for j in range(len(logm_20)):
        if logm_20[j] <= ms[i]:
            num_gals += 1
    mass_cdf_20[i] = num_gals/tot_gals_20

mass_cdf_05 = np.zeros(len(ms))
logm_05 = all_logm_init[inds_05[0]:inds_05[1]]
tot_gals_05 = len(logm_05)
for i in range(len(ms)):
    num_gals = 0
    for j in range(len(logm_05)):
        if logm_05[j] <= ms[i]:
            num_gals += 1
    mass_cdf_05[i] = num_gals/tot_gals_05
    
mass_cdf_95 = np.zeros(len(ms))
logm_95 = all_logm_init[inds_95[0]:inds_95[1]]
tot_gals_95 = len(logm_95)
for i in range(len(ms)):
    num_gals = 0
    for j in range(len(logm_95)):
        if logm_95[j] <= ms[i]:
            num_gals += 1
    mass_cdf_95[i] = num_gals/tot_gals_95
    
mass_cdf_17 = np.zeros(len(ms))
logm_17 = all_logm_init[inds_17[0]:inds_17[1]]
tot_gals_17 = len(logm_17)
for i in range(len(ms)):
    num_gals = 0
    for j in range(len(logm_17)):
        if logm_17[j] <= ms[i]:
            num_gals += 1
    mass_cdf_17[i] = num_gals/tot_gals_17
    
mass_cdf_22 = np.zeros(len(ms))
logm_22 = all_logm_init[inds_22[0]:inds_22[1]]
tot_gals_22 = len(logm_22)
for i in range(len(ms)):
    num_gals = 0
    for j in range(len(logm_22)):
        if logm_22[j] <= ms[i]:
            num_gals += 1
    mass_cdf_22[i] = num_gals/tot_gals_22


plt.figure(dpi=500)
plt.title('CDF of Galaxy Stellar Masses (Our Sample)')
plt.plot(ms, mass_cdf, color='k', label='ALL')
plt.plot(ms, mass_cdf_20, label='Pechetti+20')
plt.plot(ms, mass_cdf_05, label='Lauer+05')
plt.plot(ms, mass_cdf_95, label='Lauer+95')
plt.plot(ms, mass_cdf_17, label='Pechetti+17')
plt.plot(ms, mass_cdf_22, label='Hoyer+22')
plt.xlabel('log(M$_\star$ [M$_\odot$])')
plt.legend(loc='upper left', fontsize=7)


#%%

fig, ax = plt.subplots(dpi=500)

ax.plot(ms, mass_cdf, color='k', label='ALL')
ax.plot(ms, mass_cdf_20, label='Pechetti+20')
ax.plot(ms, mass_cdf_05, label='Lauer+05')
ax.plot(ms, mass_cdf_95, label='Lauer+95')
ax.plot(ms, mass_cdf_17, label='Pechetti+17')
ax.plot(ms, mass_cdf_22, label='Hoyer+22')
ax.set_xlabel('log(M$_\star$ [M$_\odot$])')
ax.legend(loc='center left', fontsize=7)
ax2 = ax.twinx()
ax.errorbar(bin_means_corrected, nuc_frac_corrected, xerr=bin_widths_corrected,
             linestyle='', marker='.', color='b')
ax.set_xlabel('log(M$_\star$ [M$_\odot$])')
ax2.set_ylabel('Nucleation Fraction')

# =============================================================================
# ax2 = ax.twinx()
# ax2.hist(all_logm, bins=bins, color='0.3', alpha=0.5)
# ax2.set_ylabel('# of galaxies')
# 
# ax.text(0.2, 0.65, '{}/{} galaxies'.format(len(all_logm),len(all_logm_init)),color='k', transform=plt.gcf().transFigure, size=14)
# =============================================================================



# =============================================================================
# =============================================================================
# =============================================================================

#%%



# =============================================================================
# ======================== STORE DATA =============================
# =============================================================================


# Concatente renuka data and our data for storage
all_gal_names = np.concatenate((names_r_corrected, names_corrected))
all_gal_types = np.concatenate((morphs_r, morphs_corrected))
all_logmass = np.concatenate((np.log10(gal_mass), logmass_corrected)) 
all_slopes = np.concatenate((slopes_r, slopes_corrected)) 
all_cen_dens = np.concatenate((cen_dens, np.log10(cen_dens_corrected))) 
all_NSC_comp = np.concatenate((np.ones(len(names_r_corrected)), NSC_comp_corrected))
all_gi_colors = np.concatenate((gi_color,gi_color_corrected))
all_vmags = np.concatenate((vmags,vmags_corrected))
all_dists = np.concatenate((dists,dists_corrected))
all_RA = np.concatenate((RAs,RAs_corrected))
all_Dec = np.concatenate((DECs,Decs_corrected))
all_galtype = np.concatenate((galtypes,galtype_corrected))


c1 = fits.Column(name='name', array=all_gal_names, format='10A')
c7 = fits.Column(name='RA', array=all_RA, format='D', unit='deg')
c8 = fits.Column(name='Dec', array=all_Dec, format='D', unit='deg')
c9 = fits.Column(name='g_i', array=all_gi_colors, format='D', unit='mag')
c10 = fits.Column(name='vmag', array=all_vmags, format='D', unit='mag')
c11 = fits.Column(name='dist', array=all_dists, format='D', unit='Mpc')
c2 = fits.Column(name='type', array=all_gal_types, format='I', unit='0=et,1=lt')
c3 = fits.Column(name='logmass', array=all_logmass, format='D', unit='M_sol')
c4 = fits.Column(name='slope', array=all_slopes, format='D')
c5 = fits.Column(name='cen_dens', array=all_cen_dens, format='D', unit='M_sol/pc^3')
c6 = fits.Column(name='NSC_comp', array=all_NSC_comp, format='I')


t = fits.BinTableHDU.from_columns([c1,c7,c8,c9,c10,c11,c2,c3,c4,c5,c6])
t.writeto('../Result_Tables/final_gal_data.fits',overwrite=True)


# =============================================================================
# =============================================================================
# =============================================================================

#%%

# =============================================================================
# ================= Compare our galaxies with Davis's full sample =============
# =============================================================================

#################### COLOR-MASS PLOTS AND TYPE COMPARISON #####################

# read in data from David's Table
filename = '../Data_Sets/David_Tables/sample_4.fits'
# fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
hdul = fits.open(filename)  # open a FITS file
data = hdul[1].data
hdul.close()

every_mass = data['logmass']
every_gi_color = data['gi_color']
z = np.where(np.isnan(every_mass) == False)[0]
y = np.where(np.isnan(every_gi_color) == False)[0]
v = np.intersect1d(z,y)
every_mass = every_mass[v]
every_gi_color = every_gi_color[v]
early = np.where(data['best_type'][v]=='early')[0]
late = np.where(data['best_type'][v]=='late')[0]

### set number of bins & number of contour lines used in plot
bins=15
levs=6
    ### set [mass],[color] ranges used for contours
contour_range = [[6,12],[0,2]]
    ### create bins and counts for both 2d histograms
counts_early,xb_early,yb_early = np.histogram2d(every_mass[early],every_gi_color[early],
                                                bins=bins,range=contour_range)
counts_late,xb_late,yb_late = np.histogram2d(every_mass[late],every_gi_color[late],
                                             bins=[xb_early,yb_early])

t = np.where(all_gi_colors != -999)
g = np.where(np.isnan(all_gi_colors) == False)
f = np.where(np.isnan(all_gi_colors) == True)
h = np.intersect1d(t,g)

#%%
### TYPE COMAPRISON ###
# plot and compare TYPES between david and us
plt.figure(dpi=500)
plt.plot(all_gal_types-all_galtype, 'db')
#plt.xlim(-.2,1.2)
#plt.ylim(-.2,1.2)
plt.ylabel("Current Type - David's Type")

#%%

### COLOR-MASS PLOTS ###
nils_mass_inds = [98,100,101,104,107]
sat_inds = [60,76,82,68,91]

xlim = (6,11.5)
ylim = (0.55,1.4)
alp_contour=0.3

fig, ax = plt.subplots(dpi=500)
#ax = pu.density_scatter(every_mass, every_gi_color)
plt.contour(counts_early.T,extent=[xb_early.min(),xb_early.max(),yb_early.min(),yb_early.max()],levels=levs,
            linewidths=2,alpha=alp_contour,cmap='Reds')
plt.contour(counts_late.T,extent=[xb_late.min(),xb_late.max(),yb_late.min(),yb_late.max()],levels=levs,
            linewidths=2,alpha=alp_contour,cmap='Blues')
ax.plot(all_logmass[h], all_gi_colors[h],linestyle='',marker='*')

# highlight the galaxies with saturated ACS nuclear images
ax.plot(all_logmass[np.intersect1d(h,sat_inds)], all_gi_colors[np.intersect1d(h,sat_inds)],linestyle='',marker='*')

ax.set_xlabel('log(M$_{gal}$ [M$_\odot$])')
ax.set_ylabel('Galaxy g-i from David')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.savefig('../Plots/color_mass_of_sample/color_mass.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


fig, ax = plt.subplots(dpi=500)
#ax = pu.density_scatter(every_mass, every_gi_color)
plt.contour(counts_early.T,extent=[xb_early.min(),xb_early.max(),yb_early.min(),yb_early.max()],levels=levs,
            linewidths=2,alpha=alp_contour,cmap='Reds')
plt.contour(counts_late.T,extent=[xb_late.min(),xb_late.max(),yb_late.min(),yb_late.max()],levels=levs,
            linewidths=2,alpha=alp_contour,cmap='Blues')
w = ax.scatter(all_logmass[h], all_gi_colors[h],c=all_cen_dens[h],cmap='viridis',marker='*')

ax.set_xlabel('log(M$_{gal}$ [M$_\odot$])')
ax.set_ylabel('Galaxy g-i from David')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.colorbar(w,label='log($\\rho_{5pc}$ [M$_\odot/pc^3$])')
plt.savefig('../Plots/color_mass_of_sample/color_mass_w_dens.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

fig, ax = plt.subplots(dpi=500)
#ax = pu.density_scatter(every_mass, every_gi_color)
plt.contour(counts_early.T,extent=[xb_early.min(),xb_early.max(),yb_early.min(),yb_early.max()],levels=levs,
            linewidths=2,alpha=alp_contour,cmap='Reds')
plt.contour(counts_late.T,extent=[xb_late.min(),xb_late.max(),yb_late.min(),yb_late.max()],levels=levs,
            linewidths=2,alpha=alp_contour,cmap='Blues')
w = ax.scatter(all_logmass[h], all_gi_colors[h],c=all_slopes[h],cmap='viridis',marker='*')
ax.set_xlabel('log(M$_{gal}$ [M$_\odot$])')
ax.set_ylabel('Galaxy g-i from David')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.colorbar(w, label='$\gamma$')
plt.savefig('../Plots/color_mass_of_sample/color_mass_w_slope.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

fig, ax = plt.subplots(dpi=500)
#ax = pu.density_scatter(every_mass, every_gi_color)
plt.contour(counts_early.T,extent=[xb_early.min(),xb_early.max(),yb_early.min(),yb_early.max()],levels=levs,
            linewidths=2,alpha=alp_contour,cmap='Reds')
plt.contour(counts_late.T,extent=[xb_late.min(),xb_late.max(),yb_late.min(),yb_late.max()],levels=levs,
            linewidths=2,alpha=alp_contour,cmap='Blues')
w = ax.scatter(all_logmass[h], all_gi_colors[h],c=all_gal_types[h],cmap='bwr_r',marker='*')
ax.set_xlabel('log(M$_{gal}$ [M$_\odot$])')
ax.set_ylabel('Galaxy g-i from David')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.colorbar(w, label='Type: 0 = Early, 1 = Late')
plt.savefig('../Plots/color_mass_of_sample/color_mass_w_type.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# =============================================================================
# =============================================================================
# =============================================================================


#%%


# =============================================================================
# =========== SEPARATE DENSITY MEASURMENTS BASED ON TYPE ======================
# =============================================================================

# divide samples based on type for plotting
lt_inds_r = np.where(morphs_r == 1.)
et_inds_r = np.where(morphs_r == 0.)
gal_mass_lt = gal_mass[lt_inds_r]
gal_mass_et = gal_mass[et_inds_r]
slopes_lt = slopes_r[lt_inds_r]
slopes_err_lt = slopes_r_e[lt_inds_r]
slopes_et = slopes_r[et_inds_r]
slopes_err_et = slopes_r_e[et_inds_r]
cen_dens_lt = cen_dens[lt_inds_r]
cen_dens_err_lt = cen_dens_e[lt_inds_r]
cen_dens_et = cen_dens[et_inds_r]
cen_dens_err_et = cen_dens_e[et_inds_r]
names_r_et = names_r_corrected[et_inds_r]
names_r_lt = names_r_corrected[lt_inds_r]

lt_inds_c = np.where(morphs_corrected == 1.)
et_inds_c = np.where(morphs_corrected == 0.)
names_c_et = names_corrected[et_inds_c]
names_c_lt = names_corrected[lt_inds_c]
logmass_c_lt = logmass_corrected[lt_inds_c]
logmass_c_et = logmass_corrected[et_inds_c]
slopes_c_lt = slopes_corrected[lt_inds_c]
slopes_c_err_lt = np.ones_like(slopes_c_lt)*slope_err
slopes_c_et = slopes_corrected[et_inds_c]
slopes_c_err_et = np.ones_like(slopes_c_et)*slope_err
cen_dens_c_lt = cen_dens_corrected[lt_inds_c]
cen_dens_c_err_lt = np.ones_like(cen_dens_c_lt)*cen_dens_err
cen_dens_c_et = cen_dens_corrected[et_inds_c]
cen_dens_c_err_et = np.ones_like(cen_dens_c_et)*cen_dens_err

# isolate the nucleated galaxy masses
q = np.where(NSC_comp_corrected == 1)
logmass_c_ltnsc = logmass_corrected[np.intersect1d(lt_inds_c,q)]
logmass_c_etnsc = logmass_corrected[np.intersect1d(et_inds_c,q)]
all_ltnsc_mass = np.concatenate((np.log10(gal_mass_lt),logmass_c_ltnsc))
all_etnsc_mass = np.concatenate((np.log10(gal_mass_et),logmass_c_etnsc))

# power-law function
def power_law(x,a,b):
    return a*x+b

# concatenate Renuka's data with ours
all_masses = np.concatenate((np.log10(gal_mass),logmass_corrected))
all_lt_mass = np.concatenate((np.log10(gal_mass_lt),logmass_c_lt))
all_lt_slopes = np.concatenate((slopes_lt, slopes_c_lt))
all_lt_slope_errs = np.concatenate((slopes_err_lt,slopes_c_err_lt))
all_dens_lt = np.concatenate((cen_dens_lt, np.log10(cen_dens_c_lt)))
all_dens_errs_lt = np.concatenate((cen_dens_err_lt,cen_dens_err_lt))
all_et_mass = np.concatenate((np.log10(gal_mass_et),logmass_c_et))
all_et_slopes = np.concatenate((slopes_et, slopes_c_et))
all_et_slope_errs = np.concatenate((slopes_err_et,slopes_c_err_et))
all_dens_et = np.concatenate((cen_dens_et, np.log10(cen_dens_c_et)))
all_dens_errs_et = np.concatenate((cen_dens_err_et,cen_dens_c_err_et))
all_et_names = np.concatenate((names_r_et,names_c_et))
all_lt_names = np.concatenate((names_r_lt,names_c_lt))


# define mass offset for the y-intercept for relation fits
int_offset = 9.0

# sort the arrays based on mass
lt_sort = np.argsort(all_lt_mass)
et_sort = np.argsort(all_et_mass)
all_et_names_sort = all_et_names[et_sort]
ltm = all_lt_mass[lt_sort]-int_offset # adjust intercept fit
etm = all_et_mass[et_sort]-int_offset # adjust intercept fit
ltma = all_lt_mass[lt_sort] # for plotting
etma = all_et_mass[et_sort] # for plotting
lts = all_lt_slopes[lt_sort]
ets = all_et_slopes[et_sort]
ltd = all_dens_lt[lt_sort]
etd = all_dens_et[et_sort]
ltserr = all_lt_slope_errs[lt_sort]
etserr = all_et_slope_errs[et_sort]
ltderr = all_dens_errs_lt[lt_sort]
etderr = all_dens_errs_et[et_sort]

# assume 30% errors on the masses (symmetric so using larger of 2 directions)
ltmerr = ltma-np.log10(10**ltma - 10**ltma*.3)
etmerr = etma-np.log10(10**etma - 10**etma*.3)

# set the maximum galaxy masses to be considered in the fits
max_mass = 10.5
max_lt = u.find_nearest(ltma, max_mass) + 1
max_et = u.find_nearest(etma, max_mass) + 1

# get the indices of Lauer+95 galaxies to highlight them
names_95 = np.array(['VCC 1440', 'NGC 4464', 'VCC 1545', 'NGC 4467', 'VCC 1627', 'VCC 1199'])
six_inds = np.zeros((6)).astype(int)
for i in range(6):
    six_inds[i] = np.where(all_et_names_sort[0:max_et] == names_95[i])[0]


#%%

def histogram(array, bins, norm=False):
    h = np.histogram(array, bins = bins)
    x = (h[1][1:]+h[1][:-1])/2
    if norm: return x, x[1]-x[0], h[0] / sum(h[0])
    else: return x, x[1]-x[0], h[0]
    
def fit_linmix(x, y, sig_x, sig_y):
    lm = linmix.LinMix(x, y, xsig=sig_x, ysig=sig_y, nchains=8)
    lm.run_mcmc(silent=True)
    return lm.chain['alpha'], lm.chain['beta'], lm.chain['sigsqr'], lm.chain['corr']


#%%

######### MASS vs. SLOPE LT #########

Alts, Blts, Slts, Clts  = fit_linmix(ltm[0:max_lt],lts[0:max_lt],ltmerr[0:max_lt],ltserr[0:max_lt])

a_lts, ae_lts, b_lts, be_lts = Alts.mean(), Alts.std(), Blts.mean(), Blts.std()
p_linmix = [b_lts, be_lts, a_lts, ae_lts]
cor_lts, core_lts, cor_med_lts, cor_q25_lts, cor_q75_lts = Clts.mean(), Clts.std(), np.median(Clts), np.quantile(Clts, 0.25), np.quantile(Clts, 0.75)
cor_q25_lts, cor_q75_lts = cor_med_lts - cor_q25_lts, cor_q75_lts - cor_med_lts

kwargs = {"bins" : 100,
          "ec" : None}

legend_kw = {"frameon" : False,
             "handlelength" : 1.1,
             "fontsize" : 10.5}

# plot the MCMC results of the fit

fig, ax = plt.subplots(2,2, figsize=(8,7.5),dpi=500)
fig.suptitle("M$_*$ vs. $\gamma$ for Late-Types")
ax[1,0].hist2d(Alts, Blts, bins=(100,100), norm=LogNorm())
ax[1,0].set_xlabel("Intercept")
ax[1,0].set_ylabel("Slope")

ax[0,0].set_title("Itercept")
xA_lts, wA_lts, cA_lts = histogram(Alts, 1000, norm=True)
Amax_lts = max(cA_lts) * 1.1
ax[0,0].bar(xA_lts, cA_lts, wA_lts)
ax[0,0].plot([a_lts,a_lts],[0, Amax_lts], "--k", lw=1.5, label="$\mu$ = {0:.2f}".format(a_lts))
ax[0,0].plot([a_lts+ae_lts,a_lts+ae_lts],[0, Amax_lts], ":k", lw=1.5, label="$\sigma$ = $\pm${0:.2f}".format(ae_lts))
ax[0,0].plot([a_lts-ae_lts,a_lts-ae_lts],[0, Amax_lts], ":k", lw=1.5)
ax[0,0].set_ylim(0, Amax_lts)
#ax[0,0].set_xlim(-250,250)
ax[0,0].legend(**legend_kw, loc="upper left")

ax[1,1].set_title("Slope")
xB_lts, wB_lts, cB_lts = histogram(Blts, 1000, norm=True)
Bmax_lts = max(cB_lts) * 1.1
ax[1,1].bar(xB_lts, cB_lts, wB_lts)
ax[1,1].plot([b_lts,b_lts],[0, Bmax_lts], "--k", lw=1.5, label="$\mu$ = {0:.2f}".format(b_lts))
ax[1,1].plot([b_lts+be_lts,b_lts+be_lts],[0, Bmax_lts], ":k", lw=1.5, label="$\sigma$ = $\pm${0:.2f}".format(be_lts))
ax[1,1].plot([b_lts-be_lts,b_lts-be_lts],[0, Bmax_lts], ":k", lw=1.5)
ax[1,1].set_ylim(0, Bmax_lts)
#ax[1,1].set_xlim(-25,25)
ax[1,1].legend(**legend_kw,  loc="upper left");

ax[0,1].set_title("Correlation")
xcor_lts, wcor_lts, ccor_lts = histogram(Clts, 100, norm=True)
cormax_lts = max(ccor_lts) * 1.1
ax[0,1].bar(xcor_lts, ccor_lts, wcor_lts)
ax[0,1].plot([cor_lts,cor_lts],[0, cormax_lts], "--k", lw=1.5, label="mean = {0:.2f} $\pm$ {1:.2f}".format(cor_lts, core_lts))
ax[0,1].plot([cor_med_lts,cor_med_lts],[0, cormax_lts], "-.k", lw=1.5, label="median = {0:.2f}$_{{-{1:.2f}}}^{{+{2:.2f}}}$".format(cor_med_lts, cor_q25_lts, cor_q75_lts))
ax[0,1].set_ylim(0, cormax_lts)
ax[0,1].legend(**legend_kw, loc='upper left')

plt.savefig('../Plots/linmix_plots/lts_mcmc.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# plot the relation

color_lt = (0.9,0.3,0.)
color_et = (0.,0.,0.52)
plt.figure(dpi=500)
plt.title("M$_*$ vs. $\gamma$ for Late-Types")
plt.errorbar(ltma[0:max_lt], lts[0:max_lt], yerr=ltserr[0:max_lt], 
             xerr=ltmerr[0:max_lt], capsize=2,
             linestyle='',marker='.',color=color_lt)
plt.ylabel('$\gamma$')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.55)
plt.ylim(-4.8,0.2)
mass_array = np.arange(6,10.5,0.1)
plt.plot(mass_array, b_lts*(mass_array-int_offset)+a_lts,
         linestyle='-',color=color_et, label='linmix')
plt.legend()
plt.savefig('../Plots/linmix_plots/lts_trendlines.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


#%%

######### MASS vs. SLOPE ET #########

Aets, Bets, Sets, Cets = fit_linmix(etm[0:max_et],ets[0:max_et],etmerr[0:max_et],etserr[0:max_et])

a_ets, ae_ets, b_ets, be_ets = Aets.mean(), Aets.std(), Bets.mean(), Bets.std()
p_linmix = [b_ets, be_ets, a_ets, ae_ets]
cor_ets, core_ets, cor_med_ets, cor_q25_ets, cor_q75_ets = Cets.mean(), Cets.std(), np.median(Cets), np.quantile(Cets, 0.25), np.quantile(Cets, 0.75)
cor_q25_ets, cor_q75_ets = cor_med_ets - cor_q25_ets, cor_q75_ets - cor_med_ets

kwargs = {"bins" : 100,
          "ec" : None}

legend_kw = {"frameon" : False,
             "handlelength" : 1.1,
             "fontsize" : 10.5}

# plot the MCMC results of the fit

fig, ax = plt.subplots(2,2, figsize=(8,7.5),dpi=500)
fig.suptitle("M$_*$ vs. $\gamma$ for Early-Types")
ax[1,0].hist2d(Aets, Bets, bins=(100,100), norm=LogNorm())
ax[1,0].set_xlabel("Intercept")
ax[1,0].set_ylabel("Slope")

ax[0,0].set_title("Itercept")
xA_ets, wA_ets, cA_ets = histogram(Aets, 1000, norm=True)
Amax_ets = max(cA_ets) * 1.1
ax[0,0].bar(xA_ets, cA_ets, wA_ets)
ax[0,0].plot([a_ets,a_ets],[0, Amax_ets], "--k", lw=1.5, label="$\mu$ = {0:.2f}".format(a_ets))
ax[0,0].plot([a_ets+ae_ets,a_ets+ae_ets],[0, Amax_ets], ":k", lw=1.5, label="$\sigma$ = $\pm${0:.2f}".format(ae_ets))
ax[0,0].plot([a_ets-ae_ets,a_ets-ae_ets],[0, Amax_ets], ":k", lw=1.5)
ax[0,0].set_ylim(0, Amax_ets)
#ax[0,0].set_xlim(-250,250)
ax[0,0].legend(**legend_kw, loc="upper left")

ax[1,1].set_title("Slope")
xB_ets, wB_ets, cB_ets = histogram(Bets, 1000, norm=True)
Bmax_ets = max(cB_ets) * 1.1
ax[1,1].bar(xB_ets, cB_ets, wB_ets)
ax[1,1].plot([b_ets,b_ets],[0, Bmax_ets], "--k", lw=1.5, label="$\mu$ = {0:.2f}".format(b_ets))
ax[1,1].plot([b_ets+be_ets,b_ets+be_ets],[0, Bmax_ets], ":k", lw=1.5, label="$\sigma$ = $\pm${0:.2f}".format(be_ets))
ax[1,1].plot([b_ets-be_ets,b_ets-be_ets],[0, Bmax_ets], ":k", lw=1.5)
ax[1,1].set_ylim(0, Bmax_ets)
#ax[1,1].set_xlim(-25,25)
ax[1,1].legend(**legend_kw,  loc="upper left");

ax[0,1].set_title("Correlation")
xcor_ets, wcor_ets, ccor_ets = histogram(Cets, 100, norm=True)
cormax_ets = max(ccor_ets) * 1.1
ax[0,1].bar(xcor_ets, ccor_ets, wcor_ets)
ax[0,1].plot([cor_ets,cor_ets],[0, cormax_ets], "--k", lw=1.5, label="mean = {0:.2f} $\pm$ {1:.2f}".format(cor_ets, core_ets))
ax[0,1].plot([cor_med_ets,cor_med_ets],[0, cormax_ets], "-.k", lw=1.5, label="median = {0:.2f}$_{{-{1:.2f}}}^{{+{2:.2f}}}$".format(cor_med_ets, cor_q25_ets, cor_q75_ets))
ax[0,1].set_ylim(0, cormax_ets)
ax[0,1].legend(**legend_kw, loc='upper left')

plt.savefig('../Plots/linmix_plots/ets_mcmc.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# plot the relation

color_lt = (0.9,0.3,0.)
color_et = (0.,0.,0.52)
plt.figure(dpi=500)
plt.title("M$_*$ vs. $\gamma$ for Early-Types")
#plt.plot(ltm[0:max_lt],lts[0:max_lt],linestyle='',marker='o',color=color_lt)
plt.errorbar(etma[0:max_et], ets[0:max_et], yerr=etserr[0:max_et],
             xerr=etmerr[0:max_et], capsize=2,
             linestyle='',marker='.',color=color_lt)
plt.ylabel('$\gamma$')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.55)
plt.ylim(-4.8,0.2)
mass_array = np.arange(6,10.5,0.1)
plt.plot(mass_array, b_ets*(mass_array-int_offset)+a_ets,
         linestyle='-',color=color_et, label='linmix')
plt.legend()

plt.savefig('../Plots/linmix_plots/ets_trendline.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

#%%

######### MASS vs. rho_5pc LT #########

Altd, Bltd, Sltd, Cltd = fit_linmix(ltm[0:max_lt],ltd[0:max_lt],ltmerr[0:max_lt],ltderr[0:max_lt])

a_ltd, ae_ltd, b_ltd, be_ltd = Altd.mean(), Altd.std(), Bltd.mean(), Bltd.std()
p_linmix = [b_ltd, be_ltd, a_ltd, ae_ltd]
cor_ltd, core_ltd, cor_med_ltd, cor_q25_ltd, cor_q75_ltd = Cltd.mean(), Cltd.std(), np.median(Cltd), np.quantile(Cltd, 0.25), np.quantile(Cltd, 0.75)
cor_q25_ltd, cor_q75_ltd = cor_med_ltd - cor_q25_ltd, cor_q75_ltd - cor_med_ltd

kwargs = {"bins" : 100,
          "ec" : None}

legend_kw = {"frameon" : False,
             "handlelength" : 1.1,
             "fontsize" : 10.5}

# plot the MCMC results of the fit

fig, ax = plt.subplots(2,2, figsize=(8,7.5),dpi=500)
fig.suptitle("M$_*$ vs. log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) for Late-Types")
ax[1,0].hist2d(Altd, Bltd, bins=(100,100), norm=LogNorm())
ax[1,0].set_xlabel("Intercept")
ax[1,0].set_ylabel("Slope")

ax[0,0].set_title("Itercept")
xA_ltd, wA_ltd, cA_ltd = histogram(Altd, 1000, norm=True)
Amax_ltd = max(cA_ltd) * 1.1
ax[0,0].bar(xA_ltd, cA_ltd, wA_ltd)
ax[0,0].plot([a_ltd,a_ltd],[0, Amax_ltd], "--k", lw=1.5, label="$\mu$ = {0:.2f}".format(a_ltd))
ax[0,0].plot([a_ltd+ae_ltd,a_ltd+ae_ltd],[0, Amax_ltd], ":k", lw=1.5, label="$\sigma$ = $\pm${0:.2f}".format(ae_ltd))
ax[0,0].plot([a_ltd-ae_ltd,a_ltd-ae_ltd],[0, Amax_ltd], ":k", lw=1.5)
ax[0,0].set_ylim(0, Amax_ltd)
#ax[0,0].set_xlim(-250,250)
ax[0,0].legend(**legend_kw, loc="upper left")

ax[1,1].set_title("Slope")
xB_ltd, wB_ltd, cB_ltd = histogram(Bltd, 1000, norm=True)
Bmax_ltd = max(cB_ltd) * 1.1
ax[1,1].bar(xB_ltd, cB_ltd, wB_ltd)
ax[1,1].plot([b_ltd,b_ltd],[0, Bmax_ltd], "--k", lw=1.5, label="$\mu$ = {0:.2f}".format(b_ltd))
ax[1,1].plot([b_ltd+be_ltd,b_ltd+be_ltd],[0, Bmax_ltd], ":k", lw=1.5, label="$\sigma$ = $\pm${0:.2f}".format(be_ltd))
ax[1,1].plot([b_ltd-be_ltd,b_ltd-be_ltd],[0, Bmax_ltd], ":k", lw=1.5)
ax[1,1].set_ylim(0, Bmax_ltd)
#ax[1,1].set_xlim(-25,25)
ax[1,1].legend(**legend_kw,  loc="upper left");

ax[0,1].set_title("Correlation")
xcor_ltd, wcor_ltd, ccor_ltd = histogram(Cltd, 100, norm=True)
cormax_ltd = max(ccor_ltd) * 1.1
ax[0,1].bar(xcor_ltd, ccor_ltd, wcor_ltd)
ax[0,1].plot([cor_ltd,cor_ltd],[0, cormax_ltd], "--k", lw=1.5, label="mean = {0:.2f} $\pm$ {1:.2f}".format(cor_ltd, core_ltd))
ax[0,1].plot([cor_med_ltd,cor_med_ltd],[0, cormax_ltd], "-.k", lw=1.5, label="median = {0:.2f}$_{{-{1:.2f}}}^{{+{2:.2f}}}$".format(cor_med_ltd, cor_q25_ltd, cor_q75_ltd))
ax[0,1].set_ylim(0, cormax_ltd)
ax[0,1].legend(**legend_kw, loc='upper left')

plt.savefig('../Plots/linmix_plots/ltd_mcmc.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# plot the relation

color_lt = (0.9,0.3,0.)
color_et = (0.,0.,0.52)
plt.figure(dpi=500)
#plt.plot(ltm[0:max_lt],lts[0:max_lt],linestyle='',marker='o',color=color_lt)
plt.title("M$_*$ vs. log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) for Late-Types")
plt.errorbar(ltma[0:max_lt], ltd[0:max_lt], yerr=ltderr[0:max_lt], 
             xerr=ltmerr[0:max_lt], capsize=2,
             linestyle='',marker='.',color=color_lt)
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.55)
#plt.ylim(-4.8,0.2)
mass_array = np.arange(6,10.5,0.1)
plt.plot(mass_array, b_ltd*(mass_array-int_offset)+a_ltd,
         linestyle='-',color=color_et, label='linmix')
plt.legend()

plt.savefig('../Plots/linmix_plots/ltd_trendline.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

#%%

######### MASS vs. rho_5pc ET #########

Aetd, Betd, Setd, Cetd = fit_linmix(etm[0:max_et],etd[0:max_et],etmerr[0:max_et],etderr[0:max_et])

a_etd, ae_etd, b_etd, be_etd = Aetd.mean(), Aetd.std(), Betd.mean(), Betd.std()
p_linmix = [b_etd, be_etd, a_etd, ae_etd]
cor_etd, core_etd, cor_med_etd, cor_q25_etd, cor_q75_etd = Cetd.mean(), Cetd.std(), np.median(Cetd), np.quantile(Cetd, 0.25), np.quantile(Cetd, 0.75)
cor_q25_etd, cor_q75_etd = cor_med_etd - cor_q25_etd, cor_q75_etd - cor_med_etd

kwargs = {"bins" : 100,
          "ec" : None}

legend_kw = {"frameon" : False,
             "handlelength" : 1.1,
             "fontsize" : 10.5}

# plot the MCMC results of the fit

fig, ax = plt.subplots(2,2, figsize=(8,7.5),dpi=500)
fig.suptitle("M$_*$ vs. log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) for Early-Types")
ax[1,0].hist2d(Aetd, Betd, bins=(100,100), norm=LogNorm())
ax[1,0].set_xlabel("Intercept")
ax[1,0].set_ylabel("Slope")

ax[0,0].set_title("Itercept")
xA_etd, wA_etd, cA_etd = histogram(Aetd, 1000, norm=True)
Amax_etd = max(cA_etd) * 1.1
ax[0,0].bar(xA_etd, cA_etd, wA_etd)
ax[0,0].plot([a_etd,a_etd],[0, Amax_etd], "--k", lw=1.5, label="$\mu$ = {0:.2f}".format(a_etd))
ax[0,0].plot([a_etd+ae_etd,a_etd+ae_etd],[0, Amax_etd], ":k", lw=1.5, label="$\sigma$ = $\pm${0:.2f}".format(ae_etd))
ax[0,0].plot([a_etd-ae_etd,a_etd-ae_etd],[0, Amax_etd], ":k", lw=1.5)
ax[0,0].set_ylim(0, Amax_etd)
#ax[0,0].set_xlim(-250,250)
ax[0,0].legend(**legend_kw, loc="upper left")

ax[1,1].set_title("Slope")
xB_etd, wB_etd, cB_etd = histogram(Betd, 1000, norm=True)
Bmax_etd = max(cB_etd) * 1.1
ax[1,1].bar(xB_etd, cB_etd, wB_etd)
ax[1,1].plot([b_etd,b_etd],[0, Bmax_etd], "--k", lw=1.5, label="$\mu$ = {0:.2f}".format(b_etd))
ax[1,1].plot([b_etd+be_etd,b_etd+be_etd],[0, Bmax_etd], ":k", lw=1.5, label="$\sigma$ = $\pm${0:.2f}".format(be_etd))
ax[1,1].plot([b_etd-be_etd,b_etd-be_etd],[0, Bmax_etd], ":k", lw=1.5)
ax[1,1].set_ylim(0, Bmax_etd)
#ax[1,1].set_xlim(-25,25)
ax[1,1].legend(**legend_kw,  loc="upper left");

ax[0,1].set_title("Correlation")
xcor_etd, wcor_etd, ccor_etd = histogram(Cetd, 100, norm=True)
cormax_etd = max(ccor_etd) * 1.1
ax[0,1].bar(xcor_etd, ccor_etd, wcor_etd)
ax[0,1].plot([cor_etd,cor_etd],[0, cormax_etd], "--k", lw=1.5, label="mean = {0:.2f} $\pm$ {1:.2f}".format(cor_etd, core_etd))
ax[0,1].plot([cor_med_etd,cor_med_etd],[0, cormax_etd], "-.k", lw=1.5, label="median = {0:.2f}$_{{-{1:.2f}}}^{{+{2:.2f}}}$".format(cor_med_etd, cor_q25_etd, cor_q75_etd))
ax[0,1].set_ylim(0, cormax_etd)
ax[0,1].legend(**legend_kw, loc='upper left')

plt.savefig('../Plots/linmix_plots/etd_mcmc.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# plot the relation

color_lt = (0.9,0.3,0.)
color_et = (0.,0.,0.52)
plt.figure(dpi=500)
plt.title("M$_*$ vs. log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) for Early-Types")
#plt.plot(ltm[0:max_lt],lts[0:max_lt],linestyle='',marker='o',color=color_lt)
plt.errorbar(etma[0:max_et], etd[0:max_et], yerr=etderr[0:max_et],
             xerr=etmerr[0:max_et], capsize=2,
             linestyle='',marker='.',color=color_lt)
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.55)
#plt.ylim(-4.8,0.2)
mass_array = np.arange(6,10.5,0.1)
plt.plot(mass_array, b_etd*(mass_array-int_offset)+a_etd,
         linestyle='-',color=color_et, label='linmix')
plt.legend()

plt.savefig('../Plots/linmix_plots/etd_trendline.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

#%%
# =============================================================================
# =============================================================================
# Generate pretty plots of all slopes and densities vs logmass colored by 
# morphology with trendlines from linmix + scatter
# =============================================================================
# =============================================================================

# Define the scatter values for each relation
scat_lts = np.sqrt(np.mean(Slts))
scat_ets = np.sqrt(np.mean(Sets))
scat_ltd = np.sqrt(np.mean(Sltd))
scat_etd = np.sqrt(np.mean(Setd))

color_et = (0.9,0.3,0.)
color_lt = (0.,0.,0.52)

# density plot
plt.figure(dpi=500)
plt.errorbar(ltma[0:max_lt], ltd[0:max_lt], yerr=ltderr[0:max_lt], 
             xerr=ltmerr[0:max_lt], capsize=2,
             linestyle='',marker='.',color=color_lt)


plt.errorbar(etma[0:max_et], etd[0:max_et], yerr=etderr[0:max_et],
             xerr=etmerr[0:max_et], capsize=2,
             linestyle='',marker='.',color=color_et)

plt.errorbar(etma[0:max_et][six_inds], etd[0:max_et][six_inds], yerr=etderr[0:max_et][six_inds], 
             xerr=etmerr[0:max_et][six_inds], capsize=2,
             linestyle='',marker='*',color='g')

plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.7)
plt.ylim(1,5.25)
mass_array = np.arange(6,11,0.1)
plt.plot(mass_array, b_etd*(mass_array-int_offset)+a_etd,
         linestyle='-',color=color_et)
plt.plot(mass_array, b_ltd*(mass_array-int_offset)+a_ltd,
         linestyle='-',color=color_lt)
y1_etd = (b_etd*(mass_array-int_offset)+a_etd)+scat_etd
y2_etd = (b_etd*(mass_array-int_offset)+a_etd)-scat_etd
y1_ltd = (b_ltd*(mass_array-int_offset)+a_ltd)+scat_ltd
y2_ltd = (b_ltd*(mass_array-int_offset)+a_ltd)-scat_ltd
plt.fill_between(mass_array, y1_etd, y2_etd, color=color_et, alpha=0.4)
plt.fill_between(mass_array, y1_ltd, y2_ltd, color=color_lt, alpha=0.4)

plt.plot(6.7,5,linestyle='',marker="$Early$",color=color_et,markersize=20)
plt.plot(6.7,4.8,linestyle='',marker="$Late$",color=color_lt,markersize=20)

plt.savefig('../Plots/Pechetti_and_Our_Sample/densities_v_galmass_colored_by_morph.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# slope plot
plt.figure(dpi=500)
plt.errorbar(ltma[0:max_lt], lts[0:max_lt], yerr=ltserr[0:max_lt], 
             xerr=ltmerr[0:max_lt], capsize=2,
             linestyle='',marker='.',color=color_lt)


plt.errorbar(etma[0:max_et], ets[0:max_et], yerr=etserr[0:max_et],
             xerr=etmerr[0:max_et], capsize=2,
             linestyle='',marker='.',color=color_et)

plt.errorbar(etma[0:max_et][six_inds], ets[0:max_et][six_inds], yerr=etserr[0:max_et][six_inds], 
             xerr=etmerr[0:max_et][six_inds], capsize=2,
             linestyle='',marker='*',color='g')


plt.ylabel('$\gamma$')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.7)
plt.ylim(-3.3,0.2)
mass_array = np.arange(6,11,0.1)
plt.plot(mass_array, b_ets*(mass_array-int_offset)+a_ets,
         linestyle='-',color=color_et)
plt.plot(mass_array, b_lts*(mass_array-int_offset)+a_lts,
         linestyle='-',color=color_lt, label='linmix')
y1_ets = (b_ets*(mass_array-int_offset)+a_ets)+scat_ets
y2_ets = (b_ets*(mass_array-int_offset)+a_ets)-scat_ets
y1_lts = (b_lts*(mass_array-int_offset)+a_lts)+scat_lts
y2_lts = (b_lts*(mass_array-int_offset)+a_lts)-scat_lts
plt.fill_between(mass_array, y1_ets, y2_ets, color=color_et, alpha=0.4)
plt.fill_between(mass_array, y1_lts, y2_lts, color=color_lt, alpha=0.4)

plt.plot(6.7,0,linestyle='',marker="$Early$",color=color_et,markersize=20)
plt.plot(6.7,-0.15,linestyle='',marker="$Late$",color=color_lt,markersize=20)

plt.savefig('../Plots/Pechetti_and_Our_Sample/slopes_v_galmass_colored_by_morph.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


#%%
# =============================================================================
# generate the same plots as above but with mutiple trendlines instead of scatter

# draw 100 random slopes and intercepts from each relation
num_lines = 200
bs_ets = random.sample(list(Bets),num_lines)
as_ets = random.sample(list(Aets),num_lines)
bs_lts = random.sample(list(Blts),num_lines)
as_lts = random.sample(list(Alts),num_lines)
bs_etd = random.sample(list(Betd),num_lines)
as_etd = random.sample(list(Aetd),num_lines)
bs_ltd = random.sample(list(Bltd),num_lines)
as_ltd = random.sample(list(Altd),num_lines)

# lt density plot
plt.figure(dpi=500)
mass_array = np.arange(6,11,0.1)
for i in range(num_lines):
    plt.plot(mass_array, bs_ltd[i]*(mass_array-int_offset)+as_ltd[i],
             linestyle='-',color='0.8', alpha=0.5)
    
plt.plot(mass_array, b_ltd*(mass_array-int_offset)+a_ltd,
         linestyle='-',color='k')
plt.plot(mass_array, y1_ltd, linestyle='--', color='k')
plt.plot(mass_array, y2_ltd, linestyle='--', color='k')
plt.errorbar(ltma[0:max_lt], ltd[0:max_lt], yerr=ltderr[0:max_lt], 
             xerr=ltmerr[0:max_lt], capsize=2,
             linestyle='',marker='.',color=color_lt, zorder=3)
plt.plot(0,0,linestyle='',marker='',color='w', label='N$_{lines}$'+' = {}'.format(num_lines))
plt.legend(loc='upper left', frameon=False)
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.7)

plt.savefig('../Plots/Pechetti_and_Our_Sample/densities_v_galmass_colored_by_morph_mult_trendlines_lt.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)

# et density plot
plt.figure(dpi=500)
mass_array = np.arange(6,11,0.1)
for i in range(num_lines):
    plt.plot(mass_array, bs_etd[i]*(mass_array-int_offset)+as_etd[i],
             linestyle='-',color='0.8', alpha=0.5)

plt.plot(mass_array, b_etd*(mass_array-int_offset)+a_etd,
         linestyle='-',color='k')  
plt.plot(mass_array, y1_etd, linestyle='--', color='k')
plt.plot(mass_array, y2_etd, linestyle='--', color='k')
plt.errorbar(etma[0:max_et], etd[0:max_et], yerr=etderr[0:max_et],
             xerr=etmerr[0:max_et], capsize=2,
             linestyle='',marker='.',color=color_et, zorder=3)
plt.plot(0,0,linestyle='',marker='',color='w', label='N$_{lines}$'+' = {}'.format(num_lines))
plt.legend(loc='upper left', frameon=False)
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.7)

plt.savefig('../Plots/Pechetti_and_Our_Sample/densities_v_galmass_colored_by_morph_mult_trendlines_et.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# lt slope plot
plt.figure(dpi=500)
for i in range(num_lines):
    plt.plot(mass_array, bs_lts[i]*(mass_array-int_offset)+as_lts[i],
             linestyle='-',color='0.8', alpha=0.5)

plt.plot(mass_array, b_lts*(mass_array-int_offset)+a_lts,
         linestyle='-',color='k') 
plt.plot(mass_array, y1_lts, linestyle='--', color='k')
plt.plot(mass_array, y2_lts, linestyle='--', color='k')
plt.errorbar(ltma[0:max_lt], lts[0:max_lt], yerr=ltserr[0:max_lt], 
             xerr=ltmerr[0:max_lt], capsize=2,
             linestyle='',marker='.',color=color_lt, zorder=3)
plt.plot(0,0,linestyle='',marker='',color='w', label='N$_{lines}$'+' = {}'.format(num_lines))
plt.legend(loc='upper left', frameon=False)
plt.ylabel('$\gamma$')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.7)
plt.ylim(-4.8,0.2)
mass_array = np.arange(6,11,0.1)

plt.savefig('../Plots/Pechetti_and_Our_Sample/slopes_v_galmass_colored_by_morph_mult_trendlines_lt.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# et slope plot
plt.figure(dpi=500)
for i in range(num_lines):
    plt.plot(mass_array, bs_ets[i]*(mass_array-int_offset)+as_ets[i],
             linestyle='-',color='0.8', alpha=0.5)

plt.plot(mass_array, b_ets*(mass_array-int_offset)+a_ets,
         linestyle='-',color='k') 
plt.plot(mass_array, y1_ets, linestyle='--', color='k')
plt.plot(mass_array, y2_ets, linestyle='--', color='k')
plt.errorbar(etma[0:max_et], ets[0:max_et], yerr=etserr[0:max_et],
             xerr=etmerr[0:max_et], capsize=2,
             linestyle='',marker='.',color=color_et, zorder=3)
plt.plot(0,0,linestyle='',marker='',color='w', label='N$_{lines}$'+' = {}'.format(num_lines))
plt.legend(loc='upper left', frameon=False)
plt.ylabel('$\gamma$')
plt.xlabel('log(M$_*$ [M$_\odot$])')
plt.xlim(6.4,10.7)
plt.ylim(-4.8,0.2)
mass_array = np.arange(6,11,0.1)

plt.savefig('../Plots/Pechetti_and_Our_Sample/slopes_v_galmass_colored_by_morph_mult_trendlines_et.png',
            bbox_inches='tight', pad_inches=0.1, dpi=500)


# =============================================================================
# =============================================================================
# =============================================================================


#%%

# =============================================================================
# ================== COMPARE WITH STONE&METZGER (2016) ========================
# =============================================================================

# plot the slopes from our galaxies vs stone profile slopes
k = np.where(stone_slopes_c != -99)
p = np.where(NSC_comp_c == 0)
w = np.intersect1d(k,p)

stone_slopes_corrected = stone_slopes_c[k]
slopes_c_corrected = slopes_c[k]
stone_slopes_no_nsc = stone_slopes_c[w]
slopes_no_nsc = slopes_c[w]

plt.figure(dpi=500)
plt.plot(stone_slopes_corrected,slopes_c_corrected,linestyle='',marker='.',
         color='b',label='ALL')
plt.plot([-2,0],[-2,0],linestyle='--',color='r')
plt.xlabel('$\gamma$ from Stone&Metzger 2016')
plt.ylabel('$\gamma$ from This Work')
plt.ylim(-3.6,0.1)
plt.legend()

plt.figure(dpi=500)
plt.plot(stone_slopes_no_nsc,slopes_no_nsc,linestyle='',marker='.',
         color='b', label='No NSC')
plt.plot([-2,0],[-2,0],linestyle='--',color='r')
plt.xlabel('$\gamma$ from Stone&Metzger 2016')
plt.ylabel('$\gamma$ from This Work')
plt.legend()

# plot the densities from our galaxies vs stone profile densities
stone_dens_at_5pc_c_corrected = stone_dens_at_5pc_c[k]
cen_dens_c_corrected = cen_dens_c[k]
stone_dens_at_5pc_no_nsc = stone_dens_at_5pc_c[w]
cen_dens_no_nsc = cen_dens_c[w]

plt.figure(dpi=500)
plt.plot(np.log10(stone_dens_at_5pc_c_corrected),np.log10(cen_dens_c_corrected),
         linestyle='',marker='.',color='b',label='ALL')
plt.plot([2,6],[2,6],linestyle='--',color='r')
plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from Stone&Metzger 2016')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from This Work')
plt.legend()

plt.figure(dpi=500)
plt.plot(np.log10(stone_dens_at_5pc_no_nsc),np.log10(cen_dens_no_nsc),
         linestyle='',marker='.',color='b', label='No NSC')
plt.plot([2,6],[2,6],linestyle='--',color='r')
plt.xlabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from Stone&Metzger 2016')
plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$]) from This Work')
plt.legend()

# =============================================================================
# =============================================================================
# =============================================================================

#%%

# =============================================================================
# # =============================================================================
# # recreate Renuka's 2020 density and slope plots
# 
# fig1, ax1 = plt.subplots(nrows=1,ncols=1,dpi=500)
# 
# #ax1.scatter(np.log10(gal_mass),cen_dens,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
# ax1.plot(np.log10(gal_mass),cen_dens,linestyle='',marker='.',color='b',label='Pechetti+20')
# #plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
# ax1.plot(np.sort(np.log10(gal_mass)),0.61*np.log10(np.sort(gal_mass)/10**9)+2.78,linestyle='--',color='r')
# ax1.plot(np.log10(gal_mass[8]), cen_dens[8],marker='d',color='m')
# ax1.set_ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# plt.xlabel('log(M$_*$ [M$_\odot$]')
# ax1.set_xlim(8,11)
# ax1.set_ylim(1.5,4.2)
# ax1.set_aspect('equal', adjustable='box')
# 
# 
# fig2, ax2 = plt.subplots(nrows=1,ncols=1,dpi=500)
# 
# #ax2.scatter(np.log10(gal_mass),slopes,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
# ax2.plot(np.log10(gal_mass),slopes,linestyle='',marker='.',color='b',label='Pechetti+20')
# #plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
# ax2.plot(np.sort(np.log10(gal_mass)),0.42*np.log10(np.sort(gal_mass)/10**9)-2.39,linestyle='--',color='r')
# ax2.plot(np.log10(gal_mass[8]), slopes[8],marker='d',color='m')
# ax2.set_ylabel('Power-law ($\gamma$)')
# ax2.set_xlabel('log(M$_*$ [M$_\odot$]')
# ax2.set_xlim(8,11)
# ax2.set_ylim(-3.3,-0.9)
# ax2.set_aspect('equal', adjustable='box')
# 
# 
# #%%
# 
# plt.figure(dpi=500)
# plt.scatter(np.log10(gal_mass),cen_dens,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
# #plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
# plt.plot(np.log10(gal_mass),0.592*np.log10(gal_mass/10**9)+2.819,linestyle='--',color='r')
# plt.ylabel('log($\\rho_{5pc}$ [M$_\odot$/pc$^3$])')
# plt.xlabel('log(M$_*$ [M$_\odot$]')
# plt.xlim(8,11)
# plt.ylim(1.5,4.2)
# 
# plt.figure(dpi=500)
# plt.scatter(np.log10(gal_mass),slopes,c=morphs_r,marker='.',cmap='viridis',label='Pechetti+20')
# #plt.scatter(logmass_corrected, slopes_corrected,c=morphs_corrected,marker='d',cmap='viridis',label='This Work')
# plt.plot(np.log10(gal_mass),0.963*np.log10(gal_mass/10**9)-2.762,linestyle='--',color='r')
# plt.ylabel('$\gamma$')
# plt.xlabel('log(M$_*$ [M$_\odot$]')
# plt.xlim(8,11)
# plt.ylim(-3.1,-0.9)
# 
# # =============================================================================
# =============================================================================
