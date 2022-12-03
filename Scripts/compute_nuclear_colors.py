#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 21:00:56 2022

@author: christian
"""

import csv 
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import pdb
import ssp_color_converter as ssp
#%%
all_names = []
all_mags = []
all_lines = []

# first let's read in the HST mag data
with open('../Result_Tables/nuclear_mags.csv') as f:
    reader = csv.reader(f, delimiter="\t")
    for line in reader:
        all_lines.append(line)

name_inds = np.array([1,4,7,10,13,16,19,22,25,28])
data_inds = name_inds+1

for i in range(len(name_inds)):
    all_names.append(all_lines[name_inds[i]][0].split(','))
    all_mags.append(all_lines[data_inds[i]][0].split(','))



# compute the colors for every galaxy
names = np.concatenate((np.array(all_names[0]),np.array(all_names[1]),
                              np.array(all_names[2]),np.array(all_names[3]),
                              np.array(all_names[4]),np.array(all_names[5]),
                              np.array(all_names[6]),np.array(all_names[7]),
                              np.array(all_names[8]),np.array(all_names[9])))

mags = np.concatenate((np.array(all_mags[0]),np.array(all_mags[1]),
                              np.array(all_mags[2]),np.array(all_mags[3]),
                              np.array(all_mags[4]),np.array(all_mags[5]),
                              np.array(all_mags[6]),np.array(all_mags[7]),
                              np.array(all_mags[8]),np.array(all_mags[9])))
# correct the first name to match other entry
names[0] = 'NGC1399'

names = np.concatenate((names[0:42],names[44:]))
mags = np.concatenate((mags[0:42],mags[44:])).astype(float)



wfpc2_555_lim = [0,len(all_names[0])]
wfpc2_814_lim = [wfpc2_555_lim[1],wfpc2_555_lim[1]+len(all_names[1])]
wfpc2_606_lim = [wfpc2_814_lim[1],wfpc2_814_lim[1]+len(all_names[2])-2]
wfpc2_547_lim = [wfpc2_606_lim[1],wfpc2_606_lim[1]+len(all_names[3])]
wfc3_475_lim = [wfpc2_547_lim[1],wfpc2_547_lim[1]+len(all_names[4])]
wfc3_814_lim = [wfc3_475_lim[1],wfc3_475_lim[1]+len(all_names[5])]
wfc3_547_lim = [wfc3_814_lim[1],wfc3_814_lim[1]+len(all_names[6])]
acs_475_lim = [wfc3_547_lim[1],wfc3_547_lim[1]+len(all_names[7])]
acs_850_lim = [acs_475_lim[1],acs_475_lim[1]+len(all_names[8])]
acs_814_lim = [acs_850_lim[1],acs_850_lim[1]+len(all_names[9])]

print(wfpc2_555_lim)
print(wfpc2_814_lim)
print(wfpc2_606_lim)
print(wfpc2_547_lim)
print(wfc3_475_lim)
print(wfc3_814_lim)
print(wfc3_547_lim)
print(acs_475_lim)
print(acs_850_lim)
print(acs_814_lim)

#%%
uni_names = np.unique(names)

uni_inds = np.zeros((2,len(uni_names))).astype(int)
for i in range(len(uni_names)):
    print(np.where(names == uni_names[i])[0])
    uni_inds[:,i] = np.where(names == uni_names[i])[0]

 
uni_filts = np.zeros_like(uni_inds).astype(str)
for i in range(len(uni_names)):
    for k in range(2):
        if uni_inds[k,i] >= wfpc2_555_lim[0] and uni_inds[k,i] < wfpc2_555_lim[1]:
            uni_filts[k,i] = 'WFPC2 F555W'
        elif uni_inds[k,i] >= wfpc2_814_lim[0] and uni_inds[k,i] < wfpc2_814_lim[1]:
            uni_filts[k,i] = 'WFPC2 F814W'
        elif uni_inds[k,i] >= wfpc2_606_lim[0] and uni_inds[k,i] < wfpc2_606_lim[1]:
            uni_filts[k,i] = 'WFPC2 F606W'
        elif uni_inds[k,i] >= wfpc2_547_lim[0] and uni_inds[k,i] < wfpc2_547_lim[1]:
            uni_filts[k,i] = 'WFPC2 F547M'
        elif uni_inds[k,i] >= wfc3_475_lim[0] and uni_inds[k,i] < wfc3_475_lim[1]:
            uni_filts[k,i] = 'WFC3 F475W'
        elif uni_inds[k,i] >= wfc3_814_lim[0] and uni_inds[k,i] < wfc3_814_lim[1]:
            uni_filts[k,i] = 'WFC3 F814W'
        elif uni_inds[k,i] >= wfc3_547_lim[0] and uni_inds[k,i] < wfc3_547_lim[1]:
            uni_filts[k,i] = 'WFC3 F547M'
        elif uni_inds[k,i] >= acs_475_lim[0] and uni_inds[k,i] < acs_475_lim[1]:
            uni_filts[k,i] = 'ACS F475W'
        elif uni_inds[k,i] >= acs_850_lim[0] and uni_inds[k,i] < acs_850_lim[1]:
            uni_filts[k,i] = 'ACS F850LP'
        elif uni_inds[k,i] >= acs_814_lim[0] and uni_inds[k,i] < acs_814_lim[1]:
            uni_filts[k,i] = 'ACS F814W'


wfpc2_555_814 = []
wfc3_475_814 = []
wfpc2_606_814 = []
wfpc2_606_acs_814 = []
wfpc2_555_acs_814 = []
wfpc2_606_wfc3_814 = []
acs_475_850 = []
wfc3_547_814= []
wfpc2_547_814 = []

names_wfpc2_555_814 = []
names_wfc3_475_814 = []
names_wfpc2_606_814 = []
names_wfpc2_606_acs_814 = []
names_wfpc2_555_acs_814 = []
names_wfpc2_606_wfc3_814 = []
names_acs_475_850 = []
names_wfc3_547_814= []
names_wfpc2_547_814 = []

for i in range(len(uni_names)):
    if len(np.where(uni_filts[:,i] == 'WFPC2 F555W')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'WFPC2 F814W')[0]) != 0:
            names_wfpc2_555_814.append(uni_names[i])
            ind555 = np.where(uni_filts[:,i] == 'WFPC2 F555W')[0]
            ind814 = np.where(uni_filts[:,i] == 'WFPC2 F814W')[0]
            wfpc2_555_814.append(mags[uni_inds[ind555,i]] - mags[uni_inds[ind814,i]])
            
    if len(np.where(uni_filts[:,i] == 'WFPC2 F606W')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'WFPC2 F814W')[0]) != 0:
            names_wfpc2_606_814.append(uni_names[i])
            ind606 = np.where(uni_filts[:,i] == 'WFPC2 F606W')[0]
            ind814 = np.where(uni_filts[:,i] == 'WFPC2 F814W')[0]
            wfpc2_606_814.append(mags[uni_inds[ind606,i]] - mags[uni_inds[ind814,i]])
            
    if len(np.where(uni_filts[:,i] == 'WFC3 F475W')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'WFC3 F814W')[0]) != 0:
            names_wfc3_475_814.append(uni_names[i])
            ind475 = np.where(uni_filts[:,i] == 'WFC3 F475W')[0]
            ind814 = np.where(uni_filts[:,i] == 'WFC3 F814W')[0]
            wfc3_475_814.append(mags[uni_inds[ind475,i]] - mags[uni_inds[ind814,i]])
            
    if len(np.where(uni_filts[:,i] == 'WFPC2 F606W')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'ACS F814W')[0]) != 0:
            names_wfpc2_606_acs_814.append(uni_names[i])
            ind606 = np.where(uni_filts[:,i] == 'WFPC2 F606W')[0]
            ind814 = np.where(uni_filts[:,i] == 'ACS F814W')[0]
            wfpc2_606_acs_814.append(mags[uni_inds[ind606,i]] - mags[uni_inds[ind814,i]])
            
    if len(np.where(uni_filts[:,i] == 'WFPC2 F555W')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'ACS F814W')[0]) != 0:
            names_wfpc2_555_acs_814.append(uni_names[i])
            ind555 = np.where(uni_filts[:,i] == 'WFPC2 F555W')[0]
            ind814 = np.where(uni_filts[:,i] == 'ACS F814W')[0]
            wfpc2_555_acs_814.append(mags[uni_inds[ind555,i]] - mags[uni_inds[ind814,i]])
            
    if len(np.where(uni_filts[:,i] == 'WFPC2 F606W')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'WFC3 F814W')[0]) != 0:
            names_wfpc2_606_wfc3_814.append(uni_names[i])
            ind606 = np.where(uni_filts[:,i] == 'WFPC2 F606W')[0]
            ind814 = np.where(uni_filts[:,i] == 'WFC3 F814W')[0]
            wfpc2_606_wfc3_814.append(mags[uni_inds[ind606,i]] - mags[uni_inds[ind814,i]])
            
    if len(np.where(uni_filts[:,i] == 'ACS F475W')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'ACS F850LP')[0]) != 0:
            names_acs_475_850.append(uni_names[i])
            ind475 = np.where(uni_filts[:,i] == 'ACS F475W')[0]
            ind850 = np.where(uni_filts[:,i] == 'ACS F850LP')[0]
            acs_475_850.append(mags[uni_inds[ind475,i]] - mags[uni_inds[ind850,i]])
            
    if len(np.where(uni_filts[:,i] == 'WFC3 F547M')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'WFC3 F814W')[0]) != 0:
            names_wfc3_547_814.append(uni_names[i])
            ind547 = np.where(uni_filts[:,i] == 'WFC3 F547M')[0]
            ind814 = np.where(uni_filts[:,i] == 'WFC3 F814W')[0]
            wfc3_547_814.append(mags[uni_inds[ind547,i]] - mags[uni_inds[ind814,i]])
            
    if len(np.where(uni_filts[:,i] == 'WFPC2 F547M')[0]) != 0 and \
        len(np.where(uni_filts[:,i] == 'WFPC2 F814W')[0]) != 0:
            names_wfpc2_547_814.append(uni_names[i])
            ind547 = np.where(uni_filts[:,i] == 'WFPC2 F547M')[0]
            ind814 = np.where(uni_filts[:,i] == 'WFPC2 F814W')[0]
            wfpc2_547_814.append(mags[uni_inds[ind547,i]] - mags[uni_inds[ind814,i]])

#%%

# let's add in the (1" aperture) V-I color measurments from Lauer+05





#%%
# =============================================================================
# # define the color converion linear relations.
# # copy and paste this below from output of compute_nuclear_colors.py
# wfpc2_555_814_V_I_params =  [0.9957296144724516, 0.5773626466270642]
# wfc3_475_814_V_I_params =  [1.2156149585203622, -0.054852525589466156]
# wfpc2_606_814_V_I_params =  [0.7716626702604199, -0.013397518998200667]
# wfpc2_606_acs_814_V_I_params =  [0.7853391345759677, 0.021896301353632058]
# wfpc2_555_acs_814_V_I_params =  [1.0068702782984291, 0.01588194989316296]
# wfpc2_606_wfc3_814_V_I_params =  [0.6711475600836552, 0.025780407046723387]
# acs_475_wfc3_850lp_V_I_params =  [1.73254832785054, -0.09951748944065845]
# wfc3_547_814_V_I_params =  [0.9089925292504973, -0.008732344771068278]
# wfc3_475_814_V_I_params_out =  [1.4574010187815944, -0.1518408592626391]
# wfpc2_606_814_V_I_params_out =  [0.666803737929715, 0.042096640433163815]
# wfpc2_606_acs_814_V_I_params_out =  [0.6762171747823783, 0.07951188698427143]
# wfpc2_606_wfc3_814_V_I_params_out =  [0.6009117379852231, 0.06747119723330842]
# wfc3_475_814_V_I_break =  0.242
# wfpc2_606_814_V_I_break =  0.36599999999999966
# wfpc2_606_acs_814_V_I_break =  0.4079999999999999
# wfpc2_606_wfc3_814_V_I_break =  0.3840000000000001
# wfpc2_555_814_g_i_params =  [0.8531851742184632, 0.43985245695834335]
# wfc3_475_814_g_i_params =  [1.0117789786909421, 0.49341285915865224]
# wfpc2_606_814_g_i_params =  [0.6534035277654185, 0.3407716949883519]
# wfpc2_606_acs_814_g_i_params =  [0.6655633245745755, 0.38269694190713294]
# wfpc2_555_acs_814_g_i_params =  [0.8628021567741856, 0.4801610269795961]
# wfpc2_606_wfc3_814_g_i_params =  [0.5708626543938582, 0.3356748493009951]
# acs_475_wfc3_850lp_g_i_params =  [1.4534153879003697, 0.6946532887143566]
# wfc3_547_814_g_i_params =  [0.7380059376543604, 0.40047673502288006]
# wfpc2_555_814_g_i_params_out =  [0.6455897117522612, 0.43241153012244]
# wfc3_475_814_g_i_params_out =  [0.9727609101305608, 0.49949791477822547]
# wfpc2_606_814_g_i_params_out =  [0.44303845411728315, 0.33478477936501505]
# wfpc2_606_acs_814_g_i_params_out =  [0.44901586566404195, 0.37635045093458647]
# wfpc2_555_acs_814_g_i_params_out =  [0.6543607193094688, 0.4719726302046884]
# wfpc2_606_wfc3_814_g_i_params_out =  [0.40081630912446786, 0.3306773673720151]
# acs_475_wfc3_850lp_g_i_params_out =  [1.1506006035685796, 0.6732915262722975]
# wfc3_547_814_g_i_params_out =  [0.6063126828625447, 0.3951730024592521]
# wfpc2_555_814_g_i_break =  0.3169999999999997
# wfc3_475_814_g_i_break =  0.609
# wfpc2_606_814_g_i_break =  0.18299999999999983
# wfpc2_606_acs_814_g_i_break =  0.22099999999999964
# wfpc2_555_acs_814_g_i_break =  0.33499999999999996
# wfpc2_606_wfc3_814_g_i_break =  0.24399999999999977
# acs_475_wfc3_850lp_g_i_break =  0.42700000000000005
# wfc3_547_814_g_i_break =  0.28600000000000003
# wfc3_475_814_g_i_break =  0.609
# wfpc2_606_814_g_i_break =  0.18299999999999983
# wfpc2_606_acs_814_g_i_break =  0.22099999999999964
# wfpc2_606_wfc3_814_g_i_break =  0.24399999999999977
# =============================================================================



# convert all measured colors to V-I
V_I_names = []
V_I = []

for i in range(len(wfpc2_555_814)):
    V_I_names.append(names_wfpc2_555_814[i])
    V_I.append(ssp.wfpc2_555_814_to_V_I(wfpc2_555_814[i]))
    #V_I.append(wfpc2_555_814_V_I_params[0]*wfpc2_555_814[i]+wfpc2_555_814_V_I_params[1])

for i in range(len(wfc3_475_814)):
    V_I_names.append(names_wfc3_475_814[i])
    V_I.append(ssp.wfc3_475_814_to_V_I(wfc3_475_814[i]))
# =============================================================================
#     if wfc3_475_814[i] <= wfc3_475_814_V_I_break:
#         V_I_names.append(names_wfc3_475_814[i])
#         V_I.append(wfc3_475_814_V_I_params[0]*wfc3_475_814[i]+wfc3_475_814_V_I_params[1])
#     else:
#         V_I_names.append(names_wfc3_475_814[i])
#         V_I.append(wfc3_475_814_V_I_params_out[0]*wfc3_475_814[i]+wfc3_475_814_V_I_params_out[1])
# =============================================================================
        
for i in range(len(wfpc2_606_814)):
    V_I_names.append(names_wfpc2_606_814[i])
    V_I.append(ssp.wfpc2_606_814_to_V_I(wfpc2_606_814[i]))
# =============================================================================
#     if wfpc2_606_814[i] <= wfpc2_606_814_V_I_break:
#         V_I_names.append(names_wfpc2_606_814[i])
#         V_I.append(wfpc2_606_814_V_I_params[0]*wfpc2_606_814[i]+wfpc2_606_814_V_I_params[1])
#     else:
#         V_I_names.append(names_wfpc2_606_814[i])
#         V_I.append(wfpc2_606_814_V_I_params_out[0]*wfpc2_606_814[i]+wfpc2_606_814_V_I_params_out[1])
# =============================================================================
        
for i in range(len(wfpc2_606_acs_814)):
    V_I_names.append(names_wfpc2_606_acs_814[i])
    V_I.append(ssp.wfpc2_606_acs_814_to_V_I(wfpc2_606_acs_814[i]))
# =============================================================================
#     if wfpc2_606_acs_814[i] <= wfpc2_606_acs_814_V_I_break:
#         V_I_names.append(names_wfpc2_606_acs_814[i])
#         V_I.append(wfpc2_606_acs_814_V_I_params[0]*wfpc2_606_acs_814[i]+wfpc2_606_acs_814_V_I_params[1])
#     else:
#         V_I_names.append(names_wfpc2_606_acs_814[i])
#         V_I.append(wfpc2_606_acs_814_V_I_params_out[0]*wfpc2_606_acs_814[i]+wfpc2_606_acs_814_V_I_params_out[1])
# =============================================================================
    
for i in range(len(wfpc2_555_acs_814)):
    V_I_names.append(names_wfpc2_555_acs_814[i])
    V_I.append(ssp.wfpc2_555_acs_814_to_V_I(wfpc2_555_acs_814[i]))
    #V_I.append(wfpc2_555_acs_814_V_I_params[0]*wfpc2_555_acs_814[i]+wfpc2_555_acs_814_V_I_params[1])
    
for i in range(len(wfpc2_606_wfc3_814)):
    V_I_names.append(names_wfpc2_606_wfc3_814[i])
    V_I.append(ssp.wfpc2_606_wfc3_814_to_V_I(wfpc2_606_wfc3_814[i]))
# =============================================================================
#     if wfpc2_606_wfc3_814[i] <= wfpc2_606_wfc3_814_V_I_break:
#         V_I_names.append(names_wfpc2_606_wfc3_814[i])
#         V_I.append(wfpc2_606_wfc3_814_V_I_params[0]*wfpc2_606_wfc3_814[i]+wfpc2_606_wfc3_814_V_I_params[1])
#     else:
#         V_I_names.append(names_wfpc2_606_wfc3_814[i])
#         V_I.append(wfpc2_606_wfc3_814_V_I_params_out[0]*wfpc2_606_wfc3_814[i]+wfpc2_606_wfc3_814_V_I_params_out[1])
# =============================================================================
    
for i in range(len(acs_475_850)):
    V_I_names.append(names_acs_475_850[i])
    V_I.append(ssp.acs_475_850_to_V_I(acs_475_850[i]))
    #V_I.append(acs_475_wfc3_850lp_V_I_params[0]*acs_475_850[i]+acs_475_wfc3_850lp_V_I_params[1])
    
for i in range(len(wfc3_547_814)):
    V_I_names.append(names_wfc3_547_814[i])
    V_I.append(ssp.wfc3_547_814_to_V_I(wfc3_547_814[i]))
    #V_I.append(wfc3_547_814_V_I_params[0]*wfc3_547_814[i]+wfc3_547_814_V_I_params[1])
    


# convert all measured colors to SDSS g-i
g_i_names = []
g_i = []

for i in range(len(wfpc2_555_814)):
    g_i_names.append(names_wfpc2_555_814[i])
    g_i.append(ssp.wfpc2_555_814_to_g_i(wfpc2_555_814[i]))
# =============================================================================
#     if wfpc2_555_814[i] <= wfpc2_555_814_g_i_break:
#         g_i_names.append(names_wfpc2_555_814[i])
#         g_i.append(wfpc2_555_814_g_i_params[0]*wfpc2_555_814[i]+wfpc2_555_814_g_i_params[1])
#     else:
#         g_i_names.append(names_wfpc2_555_814[i])
#         g_i.append(wfpc2_555_814_g_i_params_out[0]*wfpc2_555_814[i]+wfpc2_555_814_g_i_params_out[1])
# =============================================================================

for i in range(len(wfc3_475_814)):
    g_i_names.append(names_wfc3_475_814[i])
    g_i.append(ssp.wfc3_475_814_to_g_i(wfc3_475_814[i]))
# =============================================================================
#     if wfc3_475_814[i] <= wfc3_475_814_g_i_break:
#         g_i_names.append(names_wfc3_475_814[i])
#         g_i.append(wfc3_475_814_g_i_params[0]*wfc3_475_814[i]+wfc3_475_814_g_i_params[1])
#     else:
#         g_i_names.append(names_wfc3_475_814[i])
#         g_i.append(wfc3_475_814_g_i_params_out[0]*wfc3_475_814[i]+wfc3_475_814_g_i_params_out[1])
# =============================================================================
    
for i in range(len(wfpc2_606_814)):
    g_i_names.append(names_wfpc2_606_814[i])
    g_i.append(ssp.wfpc2_606_814_to_g_i(wfpc2_606_814[i]))
# =============================================================================
#     if wfpc2_606_814[i] <= wfpc2_606_814_g_i_break:
#         g_i_names.append(names_wfpc2_606_814[i])
#         g_i.append(wfpc2_606_814_g_i_params[0]*wfpc2_606_814[i]+wfpc2_606_814_g_i_params[1])
#     else:
#         g_i_names.append(names_wfpc2_606_814[i])
#         g_i.append(wfpc2_606_814_g_i_params_out[0]*wfpc2_606_814[i]+wfpc2_606_814_g_i_params_out[1])
# =============================================================================
    
for i in range(len(wfpc2_606_acs_814)):
    g_i_names.append(names_wfpc2_606_acs_814[i])
    g_i.append(ssp.wfpc2_606_acs_814_to_g_i(wfpc2_606_acs_814[i]))
# =============================================================================
#     if wfpc2_606_acs_814[i] <= wfpc2_606_acs_814_g_i_break:
#         g_i_names.append(names_wfpc2_606_acs_814[i])
#         g_i.append(wfpc2_606_acs_814_g_i_params[0]*wfpc2_606_acs_814[i]+wfpc2_606_acs_814_g_i_params[1])
#     else:
#         g_i_names.append(names_wfpc2_606_acs_814[i])
#         g_i.append(wfpc2_606_acs_814_g_i_params_out[0]*wfpc2_606_acs_814[i]+wfpc2_606_acs_814_g_i_params_out[1])
# =============================================================================
    
for i in range(len(wfpc2_555_acs_814)):
    g_i_names.append(names_wfpc2_555_acs_814[i])
    g_i.append(ssp.wfpc2_555_acs_814_to_g_i(wfpc2_555_acs_814[i]))
# =============================================================================
#     if wfpc2_555_acs_814[i] <= wfpc2_555_acs_814_g_i_break:
#         g_i_names.append(names_wfpc2_555_acs_814[i])
#         g_i.append(wfpc2_555_acs_814_g_i_params[0]*wfpc2_555_acs_814[i]+wfpc2_555_acs_814_g_i_params[1])
#     else:
#         g_i_names.append(names_wfpc2_555_acs_814[i])
#         g_i.append(wfpc2_555_acs_814_g_i_params_out[0]*wfpc2_555_acs_814[i]+wfpc2_555_acs_814_g_i_params_out[1])
# =============================================================================
    
for i in range(len(wfpc2_606_wfc3_814)):
    g_i_names.append(names_wfpc2_606_wfc3_814[i])
    g_i.append(ssp.wfpc2_606_wfc3_814_to_g_i(wfpc2_606_wfc3_814[i]))
# =============================================================================
#     if wfpc2_606_wfc3_814[i] <= wfpc2_606_wfc3_814_g_i_break:
#         g_i_names.append(names_wfpc2_606_wfc3_814[i])
#         g_i.append(wfpc2_606_wfc3_814_g_i_params[0]*wfpc2_606_wfc3_814[i]+wfpc2_606_wfc3_814_g_i_params[1])
#     else:
#         g_i_names.append(names_wfpc2_606_wfc3_814[i])
#         g_i.append(wfpc2_606_wfc3_814_g_i_params_out[0]*wfpc2_606_wfc3_814[i]+wfpc2_606_wfc3_814_g_i_params_out[1])
# =============================================================================
    
for i in range(len(acs_475_850)):
    g_i_names.append(names_acs_475_850[i])
    g_i.append(ssp.acs_475_850_to_g_i(acs_475_850[i]))
# =============================================================================
#     if acs_475_850[i] <= acs_475_wfc3_850lp_g_i_break:
#         g_i_names.append(names_acs_475_850[i])
#         g_i.append(acs_475_wfc3_850lp_g_i_params[0]*acs_475_850[i]+acs_475_wfc3_850lp_g_i_params[1])
#     else:
#         g_i_names.append(names_acs_475_850[i])
#         g_i.append(acs_475_wfc3_850lp_g_i_params_out[0]*acs_475_850[i]+acs_475_wfc3_850lp_g_i_params_out[1])
# =============================================================================
    
for i in range(len(wfc3_547_814)):
    g_i_names.append(names_wfc3_547_814[i])
    g_i.append(ssp.wfc3_547_814_to_g_i(wfc3_547_814[i]))
# =============================================================================
#     if wfc3_547_814[i] <= wfc3_547_814_g_i_break:
#         g_i_names.append(names_wfc3_547_814[i])
#         g_i.append(wfc3_547_814_g_i_params[0]*wfc3_547_814[i]+wfc3_547_814_g_i_params[1])
#     else:
#         g_i_names.append(names_wfc3_547_814[i])
#         g_i.append(wfc3_547_814_g_i_params_out[0]*wfc3_547_814[i]+wfc3_547_814_g_i_params_out[1])
# =============================================================================
    
for i in range(len(wfpc2_547_814)):
    g_i_names.append(names_wfpc2_547_814[i])
    g_i.append(ssp.wfpc2_547_814_to_g_i(wfpc2_547_814[i]))
# =============================================================================
#     if wfpc2_547_814[i] <= wfpc2_555_814_g_i_break:
#         g_i_names.append(names_wfpc2_547_814[i])
#         g_i.append(wfpc2_555_814_g_i_params[0]*wfpc2_547_814[i]+wfpc2_555_814_g_i_params[1])
#     else:
#         g_i_names.append(names_wfpc2_547_814[i])
#         g_i.append(wfpc2_555_814_g_i_params[0]*wfpc2_547_814[i]+wfpc2_555_814_g_i_params[1])
# =============================================================================


#%%
# let's add in the nuclear color measurements from Hoyer+22

# read in names for Hoyer+22 sample
table_b2_filename = '../Data_Sets/Hoyer22_data/tab_publication.ascii.csv'
# fields: ['galaxy', 'notes', 'filter', 'pa', 'l_pa', 'u_pa', 'ell', 'l_ell', 'u_ell', 'n', 'l_n', 'u_n', 'reff', 'u_reff', 'm', 'l_m', 'u_m', 'v', 'u_v', 'i', 'u_i', 'mlr', 'u_mlr', 'logm', 'u_logm']
all_names_22 = []
all_notes_22 = []
all_mags_22 = []
all_filts_22 = []

with open(table_b2_filename, 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row[0][0] != '#' and row[0][0] != 'g': # marks the start of data row
            all_names_22.append(row[0])
            all_notes_22.append(row[1])
            all_mags_22.append(row[14])
            all_filts_22.append(row[2])
          
all_names_22 = np.array(all_names_22)
all_notes_22 = np.array(all_notes_22)
all_mags_22 = np.array(all_mags_22)
all_filts_22 = np.array(all_filts_22) 
          
# remove galaxies where sersic indices were fixed or no fit was possible
keep_inds = []
for i in range(len(all_names_22)):
    if 'b' not in all_notes_22[i] and 'd' not in all_notes_22[i]:
        keep_inds.append(i)
        
names_22 = np.array(all_names_22[np.array(keep_inds).astype(int)])
notes_22 = np.array(all_notes_22[np.array(keep_inds).astype(int)])
mags_22 = np.array(all_mags_22[np.array(keep_inds).astype(int)])
filts_22 = np.array(all_filts_22[np.array(keep_inds).astype(int)])

# remove this galaxy
names_22 = names_22[np.where(names_22 != 'eso553-046')]
notes_22 = notes_22[np.where(names_22 != 'eso553-046')]
mags_22 = mags_22[np.where(names_22 != 'eso553-046')]
filts_22 = filts_22[np.where(names_22 != 'eso553-046')]
#pdb.set_trace()
names_22_uniq = np.unique(names_22)
mags_814_22 = np.zeros(len(names_22_uniq))
mags_606_22 = np.zeros(len(names_22_uniq))
h22_acs_606_814 = np.zeros(len(names_22_uniq))
for i in range(len(names_22_uniq)):
    gal_inds = np.where(names_22 == names_22_uniq[i])[0]
    for j in range(len(gal_inds)):
        if filts_22[gal_inds[j]] == 'f606w':
            mags_606_22[i] = mags_22[gal_inds[j]]
        elif filts_22[gal_inds[j]] == 'f814w':
            mags_814_22[i] = mags_22[gal_inds[j]]
    h22_acs_606_814[i] = mags_606_22[i] - mags_814_22[i]
    




# store galaxy names, mags, measured colors, V-I, and g-i in a fits file




