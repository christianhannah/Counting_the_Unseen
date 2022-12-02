#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 17:22:36 2022

@author: Decker French & Christian Hannah
"""

#astropy 5 env
import numpy as np
from astroquery.sdss import SDSS
import astropy.units as u
import requests
from urllib.error import HTTPError
import matplotlib.pyplot as plt
import mge1d_util as util
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.io import fits
import pdb
from matplotlib.backends.backend_pdf import PdfPages

# =============================================================================
# Function to parse David's Table for galaxy data based on coordinates
def get_David_gal_data(gal_name, ra1, dec1):
    
    # read in data from David's Table
    filename = '../Data_Sets/David_Tables/catalog.fits'
    # fields for this file: ['Name','Filt','MajAxis','SuBr','e_SuBr','PA','e_PA','Ell','e_Ell','recno']
    hdul = fits.open(filename)  # open a FITS file
    data = hdul[1].data
    hdul.close()
    
    ra2 = np.array(data['RA'])
    dec2 = np.array(data['Dec'])
    
    c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
    catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
    idx, d2d, d3d = match_coordinates_sky(c,catalog)

    max_sep = (1./60.) * u.degree
    if d2d > max_sep:
        print('*** ERROR ***: Max Separation Reached.')
        print(gal_name+" not found in David's table.")
        return '', -999,-999,-999, -999,-999
    
    #print('Match Success: Separation = {:.9f}'.format(d2d[0]))
    #pdb.set_trace()
    name = data['objname'][idx]
    dist = data['bestdist'][idx]
    logmass = data['logmass'][idx]
    ttype = data['t_type'][idx]
    gi_color = data['gi_color'][idx]
    if data['best_type'][idx] == 'early':
        galtype = 0
    elif data['best_type'][idx] == 'late':
        galtype = 1
    
    
    #pdb.set_trace()
    
    return name, dist, logmass, ttype, gi_color, galtype

# =============================================================================


from astropy.cosmology import FlatLambdaCDM, z_at_value
cosmo = FlatLambdaCDM(70,0.3)

zmax = z_at_value(cosmo.angular_diameter_distance, 50*u.Mpc,zmax=1.5)
zmax.value

# =============================================================================
# query = """select p.*, s.*, i.*, l.*, m.* from Galaxy g 
# join galSpecInfo s on s.specObjID = g.specObjID 
# join galSpecIndx i on i.specobjid = s.specobjid  
# join galSpecLine l on l.specobjid = i.specobjid  
# join galSpecExtra m on m.specobjid= l.specobjid  
# join specphotoall p on g.specobjid = p.specobjid
# where l.h_alpha_eqw > -3 
# and i.lick_hd_a - i.lick_hd_a_err > 4  
# and s.z<0.0118 
# and l.h_alpha_eqw_err>-1  
# and p.class = 'GALAXY' 
# and s.sn_median > 10 and l.h_alpha_eqw != 0"""
# =============================================================================


query = """select p.*, s.*, i.*, l.*, m.*, k.* from Galaxy g 
join galSpecInfo s on s.specObjID = g.specObjID 
join galSpecIndx i on i.specobjid = s.specobjid  
join galSpecLine l on l.specobjid = i.specobjid  
join galSpecExtra m on m.specobjid= l.specobjid  
join specphotoall p on g.specobjid = p.specobjid
join photoprofile k on p.objid = k.objid 
where l.h_alpha_eqw > -3 
and i.lick_hd_a - i.lick_hd_a_err > 4  
and s.z<0.0118 
and l.h_alpha_eqw_err>-1  
and p.class = 'GALAXY' 
and k.band = 4
and k.bin = 0
and s.sn_median > 10 and l.h_alpha_eqw != 0"""


res = SDSS.query_sql(query)
ra = res['ra']
dec = res['dec']

#cand = np.array([2, 3, 4, 9, 11, 12]); maybe = [1,7, 10]
ditch_inds = np.array([0,5,6,8,13])
ra = np.delete(ra, ditch_inds)
dec = np.delete(dec, ditch_inds)

#pdb.set_trace()
# get the spectra for each object
spectra = []
waves = []
mask = []
err = []
num_skipped = 0
count = 1
n_spec = len(res['specObjID'])
for i in range(n_spec):
    try:
        sp = SDSS.get_spectra(plate=res['plate'][i], mjd=res['mjd'][i], fiberID=res['fiberID'][i])[0]
        data = (sp[1].data)

        wave = (10**data.field('loglam'))
        flux = data.field('flux')
        error = data.field('ivar')
        masking = data.field('and_mask')


        mask.append(masking)
        spectra.append(flux*10**17) # convert as flux is in units of 10**-17 erg/cm^2/s/Ang
        err.append(error)
        waves.append(wave)
        count += 1
        
    except HTTPError:
        num_skipped += 1
        print("%i, %i, %i not found" % (res['plate'][i], res['mjd'][i], res['fiberID'][i]))
        continue
    except ValueError:
        num_skipped += 1
        print("%i, %i, %i ValueError" % (res['plate'][i], res['mjd'][i], res['fiberID'][i]))
        continue
    except TypeError:
        num_skipped += 1
        print("%i, %i, %i TypeError" % (res['plate'][i], res['mjd'][i], res['fiberID'][i]))
        continue

if num_skipped > 0:
    print("   %i spectra skipped" % num_skipped)


#%%

SB_nmy = res['profMean']
SB_vega = 22.5 - 2.5*np.log10(SB_nmy)
SB_AB = 0.54+SB_vega # from https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html


PSB_names = []
PSB_z_psfmags = []
PSB_zmags = []
PSB_dists = []
PSB_logmass = []
PSB_gicolor = []
PSB_type = []
PSB_ttype = []
PSB_cen_SB = []

for i in range(len(ra)):
    # returns: name, dist, logmass, vmag, gi_color, galtype
    info = get_David_gal_data('', ra[i], dec[i])     
    PSB_names.append(info[0])
    PSB_z_psfmags.append(res['psfMag_z'][i])
    PSB_zmags.append(res['fiberMag_z'][i])
    #PSB_imags.append(-2.5*np.log10(info[3])-4.19)
    PSB_dists.append(info[1])
    PSB_logmass.append(info[2])
    PSB_gicolor.append(info[4])
    PSB_type.append(info[5])
    PSB_ttype.append(info[3])
    PSB_cen_SB.append(SB_AB[i])    
    
    
# =============================================================================
# # manually add in name and dist for SDSS J080022.65+115023.2 as it wasn't found in David's table
# PSB_names[-2] = 'SDSS J080022.65+115023.2'
# PSB_dists[-2] = res['z'][-2]*299792/70
# PSB_logmass[-2] = 7
# =============================================================================
#w = np.where(np.array(PSB_names) != '')[0]
PSB_names = np.array(PSB_names)#[w]
PSB_z_psfmags = np.array(PSB_z_psfmags)#[w]
PSB_zmags = np.array(PSB_zmags)#[w]
PSB_dists = np.array(PSB_dists)#[w]
PSB_logmass = np.array(PSB_logmass)#[w]
PSB_gicolor = np.array(PSB_gicolor)#[w]
PSB_type = np.array(PSB_type)#[w]
PSB_ttype = np.array(PSB_ttype)#[w]
PSB_cen_SB = np.array(PSB_cen_SB)#[w] 
ra = np.array(ra)#[w]
dec = np.array(dec)#l[w]

for j in np.arange(len(ra)):
    # search for and save the SDSS image cutout of the galaxy
    url = 'https://www.legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra[j])+'&dec=' +str(dec[j])+'&layer=ls-dr9&pixscale=0.27&bands=grz'
    response = requests.get(url)
    if response.status_code == 200:
        newfile = '../Plots/JWST_Plots/cutouts_and_spectra/'+PSB_names[j]+'_SDSS_cutout_ra'+str(ra[j])+'dec'+str(dec[j])+'.jpg'
    with open(newfile, 'wb') as f:
        f.write(response.content)
       
    # plot the spectrum of each galaxy and save it as a jpeg
    plt.figure(dpi=500)
    plt.plot(waves[j],spectra[j])
    plt.xlabel('$\lambda$ [$\AA$]')
    plt.ylabel('Flux [erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]')
    plt.text(0.05, 0.96, PSB_names[j]+', log(M) = {:.2f}'.format(PSB_logmass[j])+
             ', Dist = {:.2f} Mpc'.format(PSB_dists[j])+', T-Type = {:.2f}'.format(PSB_ttype[j])+
             ', psfMag_z = {:.2f}'.format(PSB_z_psfmags[j])+
             ', SB(<0.22") = {:.2f}'.format(PSB_cen_SB[j]), transform=plt.gcf().transFigure, size=7)
    
    plt.savefig('../Plots/JWST_Plots/cutouts_and_spectra/'+PSB_names[j]+'_SDSS_spectrum_ra'+str(ra[j])+'dec'+str(dec[j])+'.jpg',
                     bbox_inches='tight', pad_inches=0.1, dpi=500)
    
        
        
#%%

# =============================================================================
# # search Chang+15 stellar mass catalog for SDSS J080022.65+115023.2 mass
# ### NEVERMIND, THE FIT WAS BAD FOR THIS GALAXY ###
# gal_file = '../Data_Sets/chang_2015_stellar_masses_table.fit'
# hdul = fits.open(gal_file)
# head = hdul[0].data
# data = hdul[1].data
# hdul.close()
# ra1 = np.array(data['RAJ2000'])
# dec1 = np.array(data['DEJ2000'])
#     
# mass_catalog = SkyCoord(ra=ra1, dec=dec1, unit=u.degree)
# 
# c = SkyCoord(ra=ra[-2]*u.degree, dec=dec[-2]*u.degree)
# idx, d2d, d3d = match_coordinates_sky(c,mass_catalog)
# 
# 
# 
# =============================================================================
#%%

arcsec_to_radian = 4.8481e-6
short_res = 0.031 # arcsec
long_res = 0.063 # arcsec
short_psf_fwhm = short_res*2
long_psf_fwhm = long_res*2
dist_array = np.arange(0.1, 50, 0.5)*10**6 # pc
short_pc_res = (short_res*arcsec_to_radian)*dist_array
long_pc_res = (long_res*arcsec_to_radian)*dist_array
short_psf_fwhm_pc = (short_psf_fwhm*arcsec_to_radian)*dist_array
long_psf_fwhm_pc = (long_psf_fwhm*arcsec_to_radian)*dist_array

# =============================================================================
# plt.figure(dpi=500)
# plt.plot(dist_array/10**6, short_pc_res, 'b', label='0.6-2.3 micron, pixel-size')
# plt.plot(dist_array/10**6, short_psf_fwhm_pc, '--b', label='0.6-2.3 micron, PSF FWHM')
# plt.plot(dist_array/10**6, long_pc_res, 'r', label='2.4-5.0 micron, pixel-size')
# plt.plot(dist_array/10**6, long_psf_fwhm_pc, '--r', label='2.4-5.0 micron, PSF FWHM')
# plt.xlabel('Distance [Mpc]')
# plt.ylabel('NIRCam Spatial Resolution [pc]')
# plt.legend()
# 
# =============================================================================

#e = np.where(np.array(PSB_dists) != -999)[0]
fig, ax = plt.subplots(dpi=500)
#ax.set_title('Post-Starburst Sample w/ D < 50 Mpc')
#ax2 = ax.twinx()
ax.plot(dist_array/10**6, short_pc_res, 'b', label='0.6-2.3 micron, 0.031"', zorder=10)
ax.plot(dist_array/10**6, short_psf_fwhm_pc, '--b', label='0.6-2.3 micron, PSF FWHM', zorder=10)
ax.plot(dist_array/10**6, long_pc_res, 'r', label='2.4-5.0 micron, 0.063"', zorder=10)
ax.plot(dist_array/10**6, long_psf_fwhm_pc, '--r', label='2.4-5.0 micron, PSF FWHM', zorder=10)
ax.set_ylabel('NIRCam Spatial Resolution [pc]')
#scat = ax2.scatter(PSB_dists, PSB_z_psfmags, s=15, c=PSB_logmass)
#ax2.hist(PSB_dists, bins=10, color='0.8',zorder=0)
ax.plot([PSB_dists[0],PSB_dists[0]],[(short_res*np.pi/0.648)*PSB_dists[0],(long_psf_fwhm*np.pi/0.648)*PSB_dists[0]],'k', label='Galaxy Distances')
for i in range(len(PSB_dists)):
    ax.plot([PSB_dists[i],PSB_dists[i]],[(short_res*np.pi/0.648)*PSB_dists[i],(long_psf_fwhm*np.pi/0.648)*PSB_dists[i]],'k')
ax.set_xlabel('Distance [Mpc]')
#ax2.set_ylabel('# of galaxies')
#fig.colorbar(scat, ax=ax2, location='top', label='log(M$_{*}$ [M$_\odot$])')
ax.set_xlim(0,50)
ax.plot([0,50],[5,5],'--k')
ax.legend(fontsize=7, frameon=False)
ax.set_ylim(0,27)
plt.savefig('../Plots/JWST_Plots/sample_overview.jpg',
                 bbox_inches='tight', pad_inches=0.1, dpi=500)


#%%
plt.figure(dpi=500)
plt.scatter(dec,ra,s=15,c=PSB_z_psfmags)
plt.xlabel('$\delta$ [deg]')
plt.ylabel('$\\alpha$ [deg]')
plt.colorbar(label='SDSS psfMag_z')
plt.savefig('../Plots/JWST_Plots/sky_distribution.jpg',
                 bbox_inches='tight', pad_inches=0.1, dpi=500)



#%%


# let's look at the number of JWST dithers needed for each filter
short_filts = np.array(['F070W','F090W','F115W','F140M','F150W','F162M','F164N',
                        'F150W2', 'F182M', 'F187N', 'F200W', 'F210M', 'F212N'])
short_pivots = np.array([.704,.901,1.154,1.404,1.501,1.626,1.644,1.671,
                         1.854,1.874,1.990,2.093,2.120])
long_filts = np.array(['F250M','F277W','F300M','F322W2','F323N','F335M','F356W',
                        'F360M', 'F405N', 'F410M', 'F430M', 'F444W', 'F460M',
                        'F466N', 'F470N', 'F480M'])
long_pivots = np.array([2.503,2.786,2.996,3.247,3.237,3.365,3.563,3.621,4.055,
                        4.092,4.280,4.421,4.624,4.654,4.707,4.834])

N_short = np.array((2/short_pivots)**2).astype(int)
N_long = np.array((4/long_pivots)**2).astype(int)

plt.figure(dpi=500)
plt.plot(short_pivots, N_short, '.b', label='Short')
plt.plot(long_pivots, N_long, '.r', label='Long')
plt.legend()
for i in range(len(short_filts[0:4])):
    plt.annotate(short_filts[i], (short_pivots[i], N_short[i]), textcoords='offset points',
                 xytext=(5,0))
plt.annotate(long_filts[0], (long_pivots[0], N_long[0]), textcoords='offset points',
             xytext=(-30,5))
plt.annotate(long_filts[1], (long_pivots[1], N_long[1]), textcoords='offset points',
             xytext=(0,5))

plt.xlabel('$\lambda$ [microns]')
plt.ylabel('N$_{dithers}$')



#%%


# let's plot the saturation time 
quer = """select p.* from ProfileDefs p"""
res1 = SDSS.query_sql(quer)

SB_nmy = res['profMean']
SB_vega = 22.5 - 2.5*np.log10(SB_nmy)
SB_AB = 0.54+SB_vega # from https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html



plt.figure(dpi=500)
#plt.plot(PSB_z_psfmags, PSB_cen_SB,'.')
plt.hist()




        