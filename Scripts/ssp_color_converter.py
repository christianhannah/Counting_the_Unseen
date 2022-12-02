#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 01:21:51 2022

@author: christian
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
import mge1d_util as u

# code to get color conversion relations needed from Padova SSP models


wfpc2_ssp_path = '../Data_Sets/Padova_SSPs/Padova_SSPs_Av_0_HST_WFPC2.dat.txt'
wfc3_wide_ssp_path = '../Data_Sets/Padova_SSPs/Padova_SSPs_Av_0_WFC3_UVIS_Wide.dat.txt'
acs_ssp_path = '../Data_Sets/Padova_SSPs/Padova_SSPs_Av_0_ACS_WFC.dat.txt'
johnson_ssp_path = '../Data_Sets/Padova_SSPs/Padova_SSPs_Av_0_UBVRIJHK.dat.txt'
wfc3_med_ssp_path = '../Data_Sets/Padova_SSPs/Padova_SSPs_Av_0_WFC3_UVIS_Medium.dat.txt'
sdss_ssp_path = '../Data_Sets/Padova_SSPs/Padova_SSPs_Av_0_SDSS.dat.txt'
wfc3_lp_ssp_path = '../Data_Sets/Padova_SSPs/Padova_SSPs_Av_0_WFC3_UVIS_LP.dat.txt'

#%%

# read in the SSP data for the various filters


### SDSS ###
with open(sdss_ssp_path) as f:
    all_lines = f.readlines()
    data_lines = []
    for i in range(len(all_lines)):
        if all_lines[i][0] != '#':
            data_lines.append(all_lines[i])

data_lines = np.array(data_lines)

# fields for model data
# age     Z   mbolmag  umag    gmag    rmag    imag    zmag

sdss_ages = np.zeros(len(data_lines))
sdss_zs = np.zeros(len(data_lines))
sdss_u = np.zeros(len(data_lines))
sdss_g = np.zeros(len(data_lines))
sdss_i = np.zeros(len(data_lines))
sdss_z = np.zeros(len(data_lines))
for i in range(len(data_lines)):
    s = data_lines[i].split()
    sdss_ages[i] = float(s[0])
    sdss_zs[i] = float(s[1])
    sdss_u[i] = float(s[3])
    sdss_g[i] = float(s[4])
    sdss_i[i] = float(s[6])
    sdss_z[i] = float(s[7])



### Johnson ###
with open(johnson_ssp_path) as f:
    all_lines = f.readlines()
    data_lines = []
    for i in range(len(all_lines)):
        if all_lines[i][0] != '#':
            data_lines.append(all_lines[i])

data_lines = np.array(data_lines)

# fields for model data
# age  Z  mbolmag  F170Wmag  F218Wmag  F255Wmag  F300Wmag  F336Wmag  F380Wmag  
# F439Wmag  F450Wmag  F555Wmag  F569Wmag  F606Wmag  F622Wmag  F675Wmag  F702Wmag  
# F791Wmag  F814Wmag  F850LPmag

johnson_ages = np.zeros(len(data_lines))
johnson_zs = np.zeros(len(data_lines))
johnson_B = np.zeros(len(data_lines))
johnson_V = np.zeros(len(data_lines))
johnson_R = np.zeros(len(data_lines))
johnson_I = np.zeros(len(data_lines))
for i in range(len(data_lines)):
    s = data_lines[i].split()
    johnson_ages[i] = float(s[0])
    johnson_zs[i] = float(s[1])
    johnson_B[i] = float(s[4])
    johnson_V[i] = float(s[5])
    johnson_R[i] = float(s[6])
    johnson_I[i] = float(s[7])


### WFPC2 ###
with open(wfpc2_ssp_path) as f:
    all_lines = f.readlines()
    data_lines = []
    for i in range(len(all_lines)):
        if all_lines[i][0] != '#':
            data_lines.append(all_lines[i])

data_lines = np.array(data_lines)

# fields for model data
# age  Z  mbolmag  F170Wmag  F218Wmag  F255Wmag  F300Wmag  F336Wmag  F380Wmag  
# F439Wmag  F450Wmag  F555Wmag  F569Wmag  F606Wmag  F622Wmag  F675Wmag  F702Wmag  
# F791Wmag  F814Wmag  F850LPmag

wfpc2_ages = np.zeros(len(data_lines))
wfpc2_zs = np.zeros(len(data_lines))
wfpc2_f555w = np.zeros(len(data_lines))
wfpc2_f450w = np.zeros(len(data_lines))
#wfpc2_f547m = np.zeros(len(data_lines))
wfpc2_f606w = np.zeros(len(data_lines))
wfpc2_f702w = np.zeros(len(data_lines))
wfpc2_f814w = np.zeros(len(data_lines))
for i in range(len(data_lines)):
    s = data_lines[i].split()
    wfpc2_ages[i] = float(s[0])
    wfpc2_zs[i] = float(s[1])
    wfpc2_f555w[i] = float(s[11])
    wfpc2_f450w[i] = float(s[10])
    wfpc2_f606w[i] = float(s[13])
    wfpc2_f702w[i] = float(s[16])
    wfpc2_f814w[i] = float(s[18])
    #wfpc2_f547m[i] = float(s[])

### ACS ###
with open(acs_ssp_path) as f:
    all_lines = f.readlines()
    data_lines = []
    for i in range(len(all_lines)):
        if all_lines[i][0] != '#':
            data_lines.append(all_lines[i])

data_lines = np.array(data_lines)

# fields for model data
# age   Z   mbolmag  F435Wmag  F475Wmag  F555Wmag  F606Wmag  F625Wmag  F775Wmag  
# F814Wmag

acs_ages = np.zeros(len(data_lines))
acs_zs = np.zeros(len(data_lines))
acs_f555w = np.zeros(len(data_lines))
acs_f475w = np.zeros(len(data_lines))
acs_f814w = np.zeros(len(data_lines))
acs_f814w = np.zeros(len(data_lines))
for i in range(len(data_lines)):
    s = data_lines[i].split()
    acs_ages[i] = float(s[0])
    acs_zs[i] = float(s[1])
    acs_f555w[i] = float(s[5])
    acs_f475w[i] = float(s[4])
    acs_f814w[i] = float(s[9])
    
    
    
### WFC3 Wide ###
with open(wfc3_wide_ssp_path) as f:
    all_lines = f.readlines()
    data_lines = []
    for i in range(len(all_lines)):
        if all_lines[i][0] != '#':
            data_lines.append(all_lines[i])

data_lines = np.array(data_lines)

# fields for model data
# age   Z   mbolmag  F218Wmag  F225Wmag  F275Wmag  F336Wmag  F390Wmag  F438Wmag  
# F475Wmag  F555Wmag  F606Wmag  F625Wmag  F775Wmag  F814Wmag  F105Wmag  F110Wmag  
# F125Wmag  F140Wmag  F160Wmag

wfc3w_ages = np.zeros(len(data_lines))
wfc3w_zs = np.zeros(len(data_lines))
wfc3w_f475w = np.zeros(len(data_lines))
wfc3w_f814w = np.zeros(len(data_lines))

for i in range(len(data_lines)):
    s = data_lines[i].split()
    wfc3w_ages[i] = float(s[0])
    wfc3w_zs[i] = float(s[1])
    wfc3w_f475w[i] = float(s[9])
    wfc3w_f814w[i] = float(s[13])
    
    
    
    
### WFC3 LP ###
with open(wfc3_lp_ssp_path) as f:
    all_lines = f.readlines()
    data_lines = []
    for i in range(len(all_lines)):
        if all_lines[i][0] != '#':
            data_lines.append(all_lines[i])

data_lines = np.array(data_lines)

# fields for model data
# age     Z  	 mbolmag  F200LPmag  F300Xmag  F350LPmag  F475Xmag  F600LPmag  F850LPmag

wfc3l_ages = np.zeros(len(data_lines))
wfc3l_zs = np.zeros(len(data_lines))
wfc3l_f850lp = np.zeros(len(data_lines))
for i in range(len(data_lines)):
    s = data_lines[i].split()
    wfc3l_ages[i] = float(s[0])
    wfc3l_zs[i] = float(s[1])
    wfc3l_f850lp[i] = float(s[8])
    
    
    
### WFC3 Medium ###
with open(wfc3_med_ssp_path) as f:
    all_lines = f.readlines()
    data_lines = []
    for i in range(len(all_lines)):
        if all_lines[i][0] != '#':
            data_lines.append(all_lines[i])

data_lines = np.array(data_lines)

# fields for model data
# age   Z   mbolmag  F390Mmag  F410Mmag  F467Mmag  F547Mmag  F621Mmag  F689Mmag  
# F763Mmag  F845Mmag  F098Mmag  F127Mmag  F139Mmag  F153Mmag

wfc3m_ages = np.zeros(len(data_lines))
wfc3m_zs = np.zeros(len(data_lines))
wfc3m_f547m = np.zeros(len(data_lines))
for i in range(len(data_lines)):
    s = data_lines[i].split()
    wfc3m_ages[i] = float(s[0])
    wfc3m_zs[i] = float(s[1])
    wfc3m_f547m[i] = float(s[6])
    

# compute colors for the filter combinations needed.
wfpc2_555_814 = wfpc2_f555w-wfpc2_f814w
wfc3_475_814 = wfc3w_f475w-wfc3w_f814w
wfpc2_606_814 = wfpc2_f606w-wfpc2_f814w
wfpc2_606_acs_814 = wfpc2_f606w-acs_f814w
wfpc2_555_acs_814 = wfpc2_f555w-acs_f814w
wfpc2_606_wfc3_814 = wfpc2_f606w-wfc3w_f814w
acs_475_wfc3_850lp = acs_f475w-wfc3l_f850lp
wfc3_547_814 = wfc3m_f547m-wfc3w_f814w



V_I = johnson_V-johnson_I
g_i = sdss_g-sdss_i

#%%

sort_inds = np.argsort(V_I)
V_I = V_I[sort_inds]
wfpc2_555_814_visort = wfpc2_555_814[sort_inds]
wfc3_475_814_visort = wfc3_475_814[sort_inds]
wfpc2_606_814_visort = wfpc2_606_814[sort_inds]
wfpc2_606_acs_814_visort = wfpc2_606_acs_814[sort_inds]
wfpc2_555_acs_814_visort = wfpc2_555_acs_814[sort_inds]
wfpc2_606_wfc3_814_visort = wfpc2_606_wfc3_814[sort_inds]
acs_475_wfc3_850lp_visort = acs_475_wfc3_850lp[sort_inds]
wfc3_547_814_visort = wfc3_547_814[sort_inds]

# plot the colors together that need to be converted and fit a trend line

#plt.figure(dpi=500)
#plt.plot(V_I, wfpc2_555_814_visort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F555W-F814W')
#plt.xlabel('V-I')

break_ind = 0

# fit inner relation
a_in, b_in = np.polyfit(V_I[0:break_ind+1]-0.6, wfpc2_555_814_visort[0:break_ind+1], 1)
#plt.plot(V_I[0:break_ind+1], a_in*V_I[0:break_ind+1]-0.6+b_in)

# fit outter relation
a_out, b_out = np.polyfit(V_I[break_ind:]-0.6, wfpc2_555_814_visort[break_ind:], 1)
#plt.plot(V_I[break_ind:], a_out*V_I[break_ind:]-0.6+b_out)

wfpc2_555_814_V_I_params = [a_out,b_out]


def wpc2_555_814_to_V_I(color):
    min_ind = u.find_nearest(wfpc2_555_814_visort, color-0.1)
    max_ind = u.find_nearest(wfpc2_555_814_visort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_555_814_visort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]


#%%

#plt.figure(dpi=500)
#plt.plot(V_I, wfc3_475_814_visort, linestyle='', marker='.')
#plt.ylabel('WFC3 F475W-F814W')
#plt.xlabel('V-I')

break_ind = 60


# fit inner relation
a_in, b_in = np.polyfit(V_I[0:break_ind+1], wfc3_475_814_visort[0:break_ind+1], 1)
#plt.plot(V_I[0:break_ind+1], a_in*V_I[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(V_I[break_ind:], wfc3_475_814_visort[break_ind:], 1)
#plt.plot(V_I[break_ind:], a_out*V_I[break_ind:]+b_out)

wfc3_475_814_V_I_break = wfc3_475_814_visort[break_ind]
wfc3_475_814_V_I_params = [a_in,b_in]
wfc3_475_814_V_I_params_out = [a_out,b_out]

def wfc3_475_814_to_V_I(color):
    min_ind = u.find_nearest(wfc3_475_814_visort, color-0.1)
    max_ind = u.find_nearest(wfc3_475_814_visort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfc3_475_814_visort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(V_I, wfpc2_606_814_visort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F606W-F814W')
#plt.xlabel('V-I')

break_ind = 180

# fit inner relation
a_in, b_in = np.polyfit(V_I[0:break_ind+1], wfpc2_606_814_visort[0:break_ind+1], 1)
#plt.plot(V_I[0:break_ind+1], a_in*V_I[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(V_I[break_ind:], wfpc2_606_814_visort[break_ind:], 1)
#plt.plot(V_I[break_ind:], a_out*V_I[break_ind:]+b_out)

wfpc2_606_814_V_I_break = wfpc2_606_814_visort[break_ind]
wfpc2_606_814_V_I_params = [a_in,b_in]
wfpc2_606_814_V_I_params_out = [a_out,b_out]

def wfpc2_606_814_to_V_I(color):
    min_ind = u.find_nearest(wfpc2_606_814_visort, color-0.1)
    max_ind = u.find_nearest(wfpc2_606_814_visort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_606_814_visort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(V_I, wfpc2_606_acs_814_visort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F606W - ACS F814W')
#plt.xlabel('V-I')

break_ind = 180

# fit inner relation
a_in, b_in = np.polyfit(V_I[0:break_ind+1], wfpc2_606_acs_814_visort[0:break_ind+1], 1)
#plt.plot(V_I[0:break_ind+1], a_in*V_I[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(V_I[break_ind:], wfpc2_606_acs_814_visort[break_ind:], 1)
#plt.plot(V_I[break_ind:], a_out*V_I[break_ind:]+b_out)

wfpc2_606_acs_814_V_I_break = wfpc2_606_acs_814_visort[break_ind]
wfpc2_606_acs_814_V_I_params = [a_in,b_in]
wfpc2_606_acs_814_V_I_params_out = [a_out,b_out]

def wfpc2_606_acs_814_to_V_I(color):
    min_ind = u.find_nearest(wfpc2_606_acs_814_visort, color-0.1)
    max_ind = u.find_nearest(wfpc2_606_acs_814_visort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_606_acs_814_visort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(V_I, wfpc2_555_acs_814_visort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F555W - ACS F814W')
#plt.xlabel('V-I')

break_ind = 0

# fit inner relation
a_in, b_in = np.polyfit(V_I[0:break_ind+1], wfpc2_555_acs_814_visort[0:break_ind+1], 1)
#plt.plot(V_I[0:break_ind+1], a_in*V_I[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(V_I[break_ind:], wfpc2_555_acs_814_visort[break_ind:], 1)
#plt.plot(V_I[break_ind:], a_out*V_I[break_ind:]+b_out)

wfpc2_555_acs_814_V_I_break = wfpc2_555_acs_814_visort[break_ind]
wfpc2_555_acs_814_V_I_params = [a_out,b_out]

def wfpc2_555_acs_814_to_V_I(color):
    min_ind = u.find_nearest(wfpc2_555_acs_814_visort, color-0.1)
    max_ind = u.find_nearest(wfpc2_555_acs_814_visort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_555_acs_814_visort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(V_I, wfpc2_606_wfc3_814_visort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F555W - WFC3 F814W')
#plt.xlabel('V-I')

break_ind = 200

# fit inner relation
a_in, b_in = np.polyfit(V_I[0:break_ind+1], wfpc2_606_wfc3_814_visort[0:break_ind+1], 1)
#plt.plot(V_I[0:break_ind+1], a_in*V_I[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(V_I[break_ind:], wfpc2_606_wfc3_814_visort[break_ind:], 1)
#plt.plot(V_I[break_ind:], a_out*V_I[break_ind:]+b_out)


wfpc2_606_wfc3_814_V_I_break = wfpc2_606_wfc3_814_visort[break_ind]
wfpc2_606_wfc3_814_V_I_params = [a_in,b_in]
wfpc2_606_wfc3_814_V_I_params_out = [a_out,b_out]

def wfpc2_606_wfc3_814_to_V_I(color):
    min_ind = u.find_nearest(wfpc2_606_wfc3_814_visort, color-0.1)
    max_ind = u.find_nearest(wfpc2_606_wfc3_814_visort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_606_wfc3_814_visort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(V_I, acs_475_wfc3_850lp_visort, linestyle='', marker='.')
#plt.ylabel('ACS F475W - WFC3 F850LP')
#plt.xlabel('V-I')

break_ind = 0

# fit inner relation
a_in, b_in = np.polyfit(V_I[0:break_ind+1], acs_475_wfc3_850lp_visort[0:break_ind+1], 1)
#plt.plot(V_I[0:break_ind+1], a_in*V_I[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(V_I[break_ind:], acs_475_wfc3_850lp_visort[break_ind:], 1)
#plt.plot(V_I[break_ind:], a_out*V_I[break_ind:]+b_out)

acs_475_wfc3_850lp_V_I_params = [a_out,b_out]

def acs_475_wfc3_850lp_to_V_I(color):
    min_ind = u.find_nearest(acs_475_wfc3_850lp_visort, color-0.1)
    max_ind = u.find_nearest(acs_475_wfc3_850lp_visort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], acs_475_wfc3_850lp_visort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(V_I, wfc3_547_814_visort, linestyle='', marker='.')
#plt.ylabel('WFC3 F547M-F814W')
#plt.xlabel('V-I')

break_ind = 0

# fit inner relation
a_in, b_in = np.polyfit(V_I[0:break_ind+1], wfc3_547_814_visort[0:break_ind+1], 1)
#plt.plot(V_I[0:break_ind+1], a_in*V_I[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(V_I[break_ind:], wfc3_547_814_visort[break_ind:], 1)
#plt.plot(V_I[break_ind:], a_out*V_I[break_ind:]+b_out)


wfc3_547_814_V_I_params = [a_out,b_out]

def wfc3_547_814_to_V_I(color):
    min_ind = u.find_nearest(wfc3_547_814_visort, color-0.1)
    max_ind = u.find_nearest(wfc3_547_814_visort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfc3_547_814_visort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

sort_inds = np.argsort(g_i)
g_i = g_i[sort_inds]
wfpc2_555_814_gisort = wfpc2_555_814[sort_inds]
wfc3_475_814_gisort = wfc3_475_814[sort_inds]
wfpc2_606_814_gisort = wfpc2_606_814[sort_inds]
wfpc2_606_acs_814_gisort = wfpc2_606_acs_814[sort_inds]
wfpc2_555_acs_814_gisort = wfpc2_555_acs_814[sort_inds]
wfpc2_606_wfc3_814_gisort = wfpc2_606_wfc3_814[sort_inds]
acs_475_wfc3_850lp_gisort = acs_475_wfc3_850lp[sort_inds]
wfc3_547_814_gisort = wfc3_547_814[sort_inds]



#%%
# plot the colors together that need to be converted and fit a trend line

#plt.figure(dpi=500)
#plt.plot(g_i, wfpc2_555_814_gisort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F555W-F814W')
#plt.xlabel('g_i')

break_ind = 85

# fit inner relation
a_in, b_in = np.polyfit(g_i[0:break_ind+1], wfpc2_555_814_gisort[0:break_ind+1], 1)
#plt.plot(g_i[0:break_ind+1], a_in*g_i[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(g_i[break_ind:], wfpc2_555_814_gisort[break_ind:], 1)
#plt.plot(g_i[break_ind:], a_out*g_i[break_ind:]+b_out)

wfpc2_555_814_g_i_break = wfpc2_555_814_gisort[break_ind]
wfpc2_555_814_g_i_params = [a_in,b_in]
wfpc2_555_814_g_i_params_out = [a_out,b_out]

def wpc2_555_814_to_g_i(color):
    min_ind = u.find_nearest(wfpc2_555_814_gisort, color-0.1)
    max_ind = u.find_nearest(wfpc2_555_814_gisort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_555_814_gisort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(g_i, wfc3_475_814_gisort, linestyle='', marker='.')
#plt.ylabel('WFC3 F475W-F814W')
#plt.xlabel('g_i')

break_ind = 200

# fit inner relation
a_in, b_in = np.polyfit(g_i[0:break_ind+1], wfc3_475_814_gisort[0:break_ind+1], 1)
#plt.plot(g_i[0:break_ind+1], a_in*g_i[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(g_i[break_ind:], wfc3_475_814_gisort[break_ind:], 1)
#plt.plot(g_i[break_ind:], a_out*g_i[break_ind:]+b_out)

wfc3_475_814_g_i_break = wfc3_475_814_gisort[break_ind]
wfc3_475_814_g_i_params = [a_in,b_in]
wfc3_475_814_g_i_params_out = [a_out,b_out]

def wfc3_475_814_to_g_i(color):
    min_ind = u.find_nearest(wfc3_475_814_gisort, color-0.1)
    max_ind = u.find_nearest(wfc3_475_814_gisort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfc3_475_814_gisort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]


#%%


#plt.figure(dpi=500)
#plt.plot(g_i, wfpc2_606_814_gisort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F606W-F814W')
#plt.xlabel('g_i')

break_ind = 60

# fit inner relation
a_in, b_in = np.polyfit(g_i[0:break_ind+1], wfpc2_606_814_gisort[0:break_ind+1], 1)
#plt.plot(g_i[0:break_ind+1], a_in*g_i[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(g_i[break_ind:], wfpc2_606_814_gisort[break_ind:], 1)
#plt.plot(g_i[break_ind:], a_out*g_i[break_ind:]+b_out)

wfpc2_606_814_g_i_break = wfpc2_606_814_gisort[break_ind]
wfpc2_606_814_g_i_params = [a_in,b_in]
wfpc2_606_814_g_i_params_out = [a_out,b_out]

def wfpc2_606_814_to_g_i(color):
    min_ind = u.find_nearest(wfpc2_606_814_gisort, color-0.1)
    max_ind = u.find_nearest(wfpc2_606_814_gisort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_606_814_gisort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(g_i, wfpc2_606_acs_814_gisort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F606W - ACS F814W')
#plt.xlabel('g_i')

break_ind = 60

# fit inner relation
a_in, b_in = np.polyfit(g_i[0:break_ind+1], wfpc2_606_acs_814_gisort[0:break_ind+1], 1)
#plt.plot(g_i[0:break_ind+1], a_in*g_i[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(g_i[break_ind:], wfpc2_606_acs_814_gisort[break_ind:], 1)
#plt.plot(g_i[break_ind:], a_out*g_i[break_ind:]+b_out)

wfpc2_606_acs_814_g_i_break = wfpc2_606_acs_814_gisort[break_ind]
wfpc2_606_acs_814_g_i_params = [a_in,b_in]
wfpc2_606_acs_814_g_i_params_out = [a_out,b_out]

def wfpc2_606_acs_814_to_g_i(color):
    min_ind = u.find_nearest(wfpc2_606_acs_814_gisort, color-0.1)
    max_ind = u.find_nearest(wfpc2_606_acs_814_gisort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_606_acs_814_gisort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]


#%%

#plt.figure(dpi=500)
#plt.plot(g_i, wfpc2_555_acs_814_gisort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F555W - ACS F814W')
#plt.xlabel('g_i')

break_ind = 70

# fit inner relation
a_in, b_in = np.polyfit(g_i[0:break_ind+1], wfpc2_555_acs_814_gisort[0:break_ind+1], 1)
#plt.plot(g_i[0:break_ind+1], a_in*g_i[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(g_i[break_ind:], wfpc2_555_acs_814_gisort[break_ind:], 1)
##plt.plot(g_i[break_ind:], a_out*g_i[break_ind:]+b_out)

wfpc2_555_acs_814_g_i_break = wfpc2_555_acs_814_gisort[break_ind]
wfpc2_555_acs_814_g_i_params = [a_in,b_in]
wfpc2_555_acs_814_g_i_params_out = [a_out,b_out]

def wfpc2_555_acs_814_to_g_i(color):
    min_ind = u.find_nearest(wfpc2_555_acs_814_gisort, color-0.1)
    max_ind = u.find_nearest(wfpc2_555_acs_814_gisort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_555_acs_814_gisort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]


#%%

##plt.figure(dpi=500)
#plt.plot(g_i, wfpc2_606_wfc3_814_gisort, linestyle='', marker='.')
#plt.ylabel('WFPC2 F606W - WFC3 F814W')
#plt.xlabel('g_i')

break_ind = 70

# fit inner relation
a_in, b_in = np.polyfit(g_i[0:break_ind+1], wfpc2_606_wfc3_814_gisort[0:break_ind+1], 1)
#plt.plot(g_i[0:break_ind+1], a_in*g_i[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(g_i[break_ind:], wfpc2_606_wfc3_814_gisort[break_ind:], 1)
#plt.plot(g_i[break_ind:], a_out*g_i[break_ind:]+b_out)

wfpc2_606_wfc3_814_g_i_break = wfpc2_606_wfc3_814_gisort[break_ind]
wfpc2_606_wfc3_814_g_i_params = [a_in,b_in]
wfpc2_606_wfc3_814_g_i_params_out = [a_out,b_out]

def wfpc2_606_wfc3_814_to_g_i(color):
    min_ind = u.find_nearest(wfpc2_606_wfc3_814_gisort, color-0.1)
    max_ind = u.find_nearest(wfpc2_606_wfc3_814_gisort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfpc2_606_wfc3_814_gisort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(g_i, acs_475_wfc3_850lp_gisort, linestyle='', marker='.')
#plt.ylabel('ACS F475W - WFC3 F850LP')
#plt.xlabel('g_i')

break_ind = 70

# fit inner relation
a_in, b_in = np.polyfit(g_i[0:break_ind+1], acs_475_wfc3_850lp_gisort[0:break_ind+1], 1)
#plt.plot(g_i[0:break_ind+1], a_in*g_i[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(g_i[break_ind:], acs_475_wfc3_850lp_gisort[break_ind:], 1)
#plt.plot(g_i[break_ind:], a_out*g_i[break_ind:]+b_out)

acs_475_wfc3_850lp_g_i_break = acs_475_wfc3_850lp_gisort[break_ind]
acs_475_wfc3_850lp_g_i_params = [a_in,b_in]
acs_475_wfc3_850lp_g_i_params_out = [a_out,b_out]

def acs_475_wfc3_850lp_to_g_i(color):
    min_ind = u.find_nearest(acs_475_wfc3_850lp_gisort, color-0.1)
    max_ind = u.find_nearest(acs_475_wfc3_850lp_gisort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], acs_475_wfc3_850lp_gisort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]

#%%

#plt.figure(dpi=500)
#plt.plot(g_i, wfc3_547_814_gisort, linestyle='', marker='.')
#plt.ylabel('WFC3 F547M-F814W')
#plt.xlabel('g_i')

break_ind = 80

# fit inner relation
a_in, b_in = np.polyfit(g_i[0:break_ind+1], wfc3_547_814_gisort[0:break_ind+1], 1)
#plt.plot(g_i[0:break_ind+1], a_in*g_i[0:break_ind+1]+b_in)

# fit outter relation
a_out, b_out = np.polyfit(g_i[break_ind:], wfc3_547_814_gisort[break_ind:], 1)
#plt.plot(g_i[break_ind:], a_out*g_i[break_ind:]+b_out)

wfc3_547_814_g_i_break = wfc3_547_814_gisort[break_ind]
wfc3_547_814_g_i_params = [a_in,b_in]
wfc3_547_814_g_i_params_out = [a_out,b_out]

def wfc3_547_814_to_g_i(color):
    min_ind = u.find_nearest(wfc3_547_814_gisort, color-0.1)
    max_ind = u.find_nearest(wfc3_547_814_gisort, color+0.1)
    a, b = np.polyfit(V_I[min_ind:max_ind+1], wfc3_547_814_gisort[min_ind:max_ind+1], 1)
    colors = np.arange(color-0.1,color+0.1,0.01)
    V_Is = a*colors + b
    color_ind = u.find_nearest(colors,color)
    return V_Is[color_ind]


#%%
# =============================================================================
# print()
# print('wfpc2_555_814_V_I_params = ',wfpc2_555_814_V_I_params)
# print('wfc3_475_814_V_I_params = ',wfc3_475_814_V_I_params)
# print('wfpc2_606_814_V_I_params = ',wfpc2_606_814_V_I_params)
# print('wfpc2_606_acs_814_V_I_params = ',wfpc2_606_acs_814_V_I_params)
# print('wfpc2_555_acs_814_V_I_params = ',wfpc2_555_acs_814_V_I_params)
# print('wfpc2_606_wfc3_814_V_I_params = ',wfpc2_606_wfc3_814_V_I_params)
# print('acs_475_wfc3_850lp_V_I_params = ',acs_475_wfc3_850lp_V_I_params)
# print('wfc3_547_814_V_I_params = ',wfc3_547_814_V_I_params)
# 
# 
# print('wfc3_475_814_V_I_params_out = ',wfc3_475_814_V_I_params_out)
# print('wfpc2_606_814_V_I_params_out = ',wfpc2_606_814_V_I_params_out)
# print('wfpc2_606_acs_814_V_I_params_out = ',wfpc2_606_acs_814_V_I_params_out)
# print('wfpc2_606_wfc3_814_V_I_params_out = ',wfpc2_606_wfc3_814_V_I_params_out)
# 
# print('wfc3_475_814_V_I_break = ',wfc3_475_814_V_I_break)
# print('wfpc2_606_814_V_I_break = ',wfpc2_606_814_V_I_break)
# print('wfpc2_606_acs_814_V_I_break = ',wfpc2_606_acs_814_V_I_break)
# print('wfpc2_606_wfc3_814_V_I_break = ',wfpc2_606_wfc3_814_V_I_break)
# 
# print('wfpc2_555_814_g_i_params = ',wfpc2_555_814_g_i_params)
# print('wfc3_475_814_g_i_params = ',wfc3_475_814_g_i_params)
# print('wfpc2_606_814_g_i_params = ',wfpc2_606_814_g_i_params)
# print('wfpc2_606_acs_814_g_i_params = ',wfpc2_606_acs_814_g_i_params)
# print('wfpc2_555_acs_814_g_i_params = ',wfpc2_555_acs_814_g_i_params)
# print('wfpc2_606_wfc3_814_g_i_params = ',wfpc2_606_wfc3_814_g_i_params)
# print('acs_475_wfc3_850lp_g_i_params = ',acs_475_wfc3_850lp_g_i_params)
# print('wfc3_547_814_g_i_params = ',wfc3_547_814_g_i_params)
# 
# print('wfpc2_555_814_g_i_params_out = ',wfpc2_555_814_g_i_params_out)
# print('wfc3_475_814_g_i_params_out = ',wfc3_475_814_g_i_params_out)
# print('wfpc2_606_814_g_i_params_out = ',wfpc2_606_814_g_i_params_out)
# print('wfpc2_606_acs_814_g_i_params_out = ',wfpc2_606_acs_814_g_i_params_out)
# print('wfpc2_555_acs_814_g_i_params_out = ',wfpc2_555_acs_814_g_i_params_out)
# print('wfpc2_606_wfc3_814_g_i_params_out = ',wfpc2_606_wfc3_814_g_i_params_out)
# print('acs_475_wfc3_850lp_g_i_params_out = ',acs_475_wfc3_850lp_g_i_params_out)
# print('wfc3_547_814_g_i_params_out = ',wfc3_547_814_g_i_params_out)
# 
# print('wfpc2_555_814_g_i_break = ',wfpc2_555_814_g_i_break)
# print('wfc3_475_814_g_i_break = ',wfc3_475_814_g_i_break)
# print('wfpc2_606_814_g_i_break = ',wfpc2_606_814_g_i_break)
# print('wfpc2_606_acs_814_g_i_break = ',wfpc2_606_acs_814_g_i_break)
# print('wfpc2_555_acs_814_g_i_break = ',wfpc2_555_acs_814_g_i_break)
# print('wfpc2_606_wfc3_814_g_i_break = ',wfpc2_606_wfc3_814_g_i_break)
# print('acs_475_wfc3_850lp_g_i_break = ',acs_475_wfc3_850lp_g_i_break)
# print('wfc3_547_814_g_i_break = ',wfc3_547_814_g_i_break)
# 
# print('wfc3_475_814_g_i_break = ',wfc3_475_814_g_i_break)
# print('wfpc2_606_814_g_i_break = ',wfpc2_606_814_g_i_break)
# print('wfpc2_606_acs_814_g_i_break = ',wfpc2_606_acs_814_g_i_break)
# print('wfpc2_606_wfc3_814_g_i_break = ',wfpc2_606_wfc3_814_g_i_break)
# =============================================================================
