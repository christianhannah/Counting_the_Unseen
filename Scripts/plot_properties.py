import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from astropy.io import fits
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
import numpy.linalg as npl
from scipy.stats import linregress
import linmix
from astropy.table import Table
from astropy.table import join
from scipy.odr import *
# def lst_fit(X, Y):
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
# 	var = np.dot(np.transpose(X), X)
# 	inverse_el = np.ones((2, len(X))) 
# 	inverse_el[0, :] = X
# 	[a, b] = np.dot(npl.pinv(X), Y)
# 	return Z, 1 #, corrcoef(X,Y)[1,0]

rc('font',**{'family':'sans-serif','sans-serif':['Times New Roman'],'size':25,'weight':'normal'})
rc('axes',**{'linewidth':2})
plt.rcParams.update({'mathtext.default':'regular'})
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['xtick.major.size'] = 6
# mpl.rcParams['xtick.minor.width'] = 3
# mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['ytick.major.size'] = 6
# mpl.rcParams['ytick.minor.width'] = 3
# mpl.rcParams['ytick.minor.size'] = 3



#---------------------- mag_dist.pdf ------------------------------#

# a = np.loadtxt('../datafiles/10mpc_hst.dat',dtype=str)
# b = np.loadtxt('../datafiles/list.csv',dtype=str,delimiter=',')
# filt = a[:,-1].astype(int)
# Mb = a[:,3].astype(float)
# T = a[:,4].astype(float)
# D = a[:,-5].astype(float)

# Mb_n = b[:,9].astype(float)[29:]
# T_n = b[:,15].astype(float)[29:]
# D_n = b[:,3].astype(float)[29:]

# fig,(ax1,ax2) = plt.subplots(ncols=2, nrows=1, sharey=True,figsize=(6,8) )

# ax1.plot(T[filt==0],Mb[filt==0],'o',markersize=12,color='black',fillstyle='none')
# ax1.plot(T[filt==1],Mb[filt==1],'o',label='Selected Sample',color='black',markersize=12)
# ax2.plot(D[filt==0],Mb[filt==0],'o',markersize=12,color='black',fillstyle='none')
# ax2.plot(D[filt==1],Mb[filt==1],'o',label='Selected Sample',color='black',markersize=12)
# ax1.plot(T_n,Mb_n,'*',label='Nguyen(2018) + Carson(2015)',color='gray',markersize=13)
# ax2.plot(D_n,Mb_n,'*',label='Nguyen(2018) + Carson(2015)',color='gray',markersize=13)
# ax1.set_ylabel(r'$M_B (mag)$')
# ax1.set_xlabel('T-type')
# ax2.set_xlabel('Distance (Mpc)')
# ax1.legend(frameon=False,loc='upper left',fontsize=18)
# ax1.xaxis.set_minor_locator(MultipleLocator(5))
# ax1.tick_params(which='both', width=2)
# ax1.set_ylim(-14.5,-21.5)
# ax2.set_ylim(-14.5,-21.5)
# # ax1.yaxis.set_tick_params(which='both',width=4,length=8)
# # ax2.yaxis.set_tick_params(which='both',width=4,length=8)

# plt.show()



#------------------ galmass plot -----------------------------#
# f = lambda x: sorted(range(len(x)), key=lambda i: x[i]) #sort a column and sort others based on the index
# b = np.loadtxt('../datafiles/Table1.csv',skiprows=1, dtype=str,delimiter=',')
# c = np.loadtxt('../datafiles/table_modifiedsurf.csv',skiprows=2, dtype=str,delimiter=',')
# mass_gmi = c[:,5].astype(float)
# M_cook = c[:,7].astype(float)
# idx1 = np.where(np.isnan(M_cook))[0]
# idx = np.where(~np.isnan(M_cook))[0]
# mass = c[:,8].astype(float)
# distance = b[:,1].astype(float)
# ttype=b[:,2].astype(int)
# fig,ax = plt.subplots(figsize=(14,10))

# sc = ax.scatter(distance,np.log10(mass),c=ttype,s=0,cmap='jet_r',zorder=3,)
# err = np.repeat(np.log10(1.3),29)
# err1 = -np.repeat(np.log10(0.7),29)
# plt.show()
# color = sc.get_facecolor()[idx]
# color1 = sc.get_facecolor()[idx1]
# plt.plot(distance[idx][0],np.log10(mass[idx][0]),'*',color='black',markersize=18,markeredgewidth=2,mfc='none',label='Masses from Cook (2014)')
# plt.plot(distance[idx1][0],np.log10(mass[idx1][0]),'o',color='black',markersize=13,markeredgewidth=2,mfc='none',label='Masses from color-M/L relations')

# for j in range(0,len(color)):

# 	plt.plot(distance[idx][j],np.log10(mass[idx][j]),marker='*',markeredgecolor='black',markeredgewidth=1,markersize=26,color =color[j])
# for j in range(0,len(color1)):

# 	plt.plot(distance[idx1][j],np.log10(mass[idx1][j]),marker='o',markeredgecolor='black',markeredgewidth=1,markersize=19,color =color1[j])

# ax.errorbar(distance,np.log10(mass),yerr=[err1,err],fmt=None,marker=None,ecolor='black',elinewidth = 2,capsize=2,capthick=2,zorder=0)
# ax.set_ylabel(r'log (Galaxy stellar Mass)$[M_\odot]$')
# ax.set_xlabel('Distance (Mpc)')
# cbar = plt.colorbar(sc)
# cbar.set_label('T-type',fontsize=35)
# plt.text(11.3,7.6,'Early-Type',fontsize=30)
# plt.text(11.3,11.2,'Late-Type',fontsize=30)
# plt.xlim(2.5,11.5)
# #ax1.yaxis.set_tick_params(width=4,length=8)
# plt.tick_params(width=4,length=12)
# plt.legend(fontsize=20)

# # !!!! ---------------------------- nscmass color - galmass --------------#
# def fit_line(xs,A,B):
# 	return A*xs+B
# def linear_func(B,x):
# 	return B[0]*x+B[1]

# a = np.loadtxt('../datafiles/table_modifiedsurf.csv',dtype=str,delimiter=',',skiprows=2)
# b = np.loadtxt('../datafiles/asu.tsv',delimiter=';')
# c = np.loadtxt('../datafiles/list.csv',dtype=str,skiprows = 1,delimiter=',')
# d = np.loadtxt('../datafiles/asu_1.tsv',delimiter=';')
# s17 = Table.read('../datafiles/spengler17.dat',format='ascii')
# c06=Table.read('../datafiles/cote06_table1.dat',format='ascii')
# acsvcs=join(s17,c06,keys='VCC',join_type='inner')
# gal = a[:,0]
# qual = a[:,-2].astype(int)
# galmass_my = a[:,8].astype(float)
# nscmass_my = a[:,12].astype(float)
# reff_my = a[:,10].astype(float)
# ttype = c[:,15][0:29].astype(int)
# galmass_d = c[:,13].astype(float)
# galmass_d = galmass_d[29:]
# nscmass_d = c[:,16].astype(float)
# nscmass_d = nscmass_d[29:]
# ttype_d = c[:,15][29:].astype(int)
# reff_d = c[:,12][29:].astype(float)
# galmass_s = s17['M_star']
# nscmass_s = s17['M_nuc']

# galmass_g = b[:,6]*10**7
# nscmass_g = b[:,3]*10**4
# reff_g = b[:,0]
# ttype_g = d[:,0]

# sjgal=Table.read('../datafiles/sanchez-janssen19_tab4.tex')
# #print("sjgal Columns: ", sjgal.colnames)
# sjnuc=Table.read('../datafiles/sanchez-janssen19_tab5.tex')
# #print("sjnuc Columns: ", sjnuc.colnames)
# #Masses determined by Roediger relation fits to all available photometry
# sjgal['g-i']=sjgal['gmag']-sjgal['imag']
# sjgal['source']='NGVS'
# sj=join(sjnuc,sjgal,keys='id')
# sj['logmnuc']=sj['logmstar_1']
# sj['logmstargal']=sj['logmstar_2']
# sj['T']=-1

# fig,(ax1,ax2) = plt.subplots(ncols=2, nrows=1,figsize=(10,8))
# fig,ax = plt.subplots()

# err = np.repeat(0.15,24)

# err1 = 0.10*reff_my[qual==0][ttype[qual==0]<0]
# err2 = 0.10*reff_my[qual==0][ttype[qual==0]>0]
# err3 = np.repeat(0.15,len(galmass_my[qual==1]))

# err_nsc1 = 0.17*nscmass_my[qual==0][ttype[qual==0]<0]
# err_nsc1 = -np.log10(1-(err_nsc1/nscmass_my[qual==0][ttype[qual==0]<0]))
# err_nsc2 = 0.17*nscmass_my[qual==0][ttype[qual==0]<0]
# err_nsc2 = np.log10(1+(err_nsc2/nscmass_my[qual==0][ttype[qual==0]<0]))

# err_nsc3 = 0.17*nscmass_my[qual==0][ttype[qual==0]>0]
# err_nsc3 = -np.log10(1-(err_nsc3/nscmass_my[qual==0][ttype[qual==0]>0]))
# err_nsc4 = 0.17*nscmass_my[qual==0][ttype[qual==0]>0]
# err_nsc4 = np.log10(1+(err_nsc4/nscmass_my[qual==0][ttype[qual==0]>0]))
# galmass_my_q = galmass_my[qual==1]
# nscmass_my_q = nscmass_my[qual==1]

# # --------- plot galmass - nscmass ----------------#
# ax1.plot(np.log10(galmass_d),np.log10(nscmass_d),'*',markersize=16,color='black',mfc='none',label='Dynamical measurements (Literature)')
# ax1.plot(np.log10(galmass_d[ttype_d<0]),np.log10(nscmass_d[ttype_d<0]),'*',markersize=16,color='red',alpha=0.5)
# ax1.plot(np.log10(galmass_d[ttype_d>0]),np.log10(nscmass_d[ttype_d>0]),'*',markersize=16,color='blue',alpha=0.5)

# ax1.plot(np.log10(galmass_g),np.log10(nscmass_g),'o',markersize=8,color='blue',alpha=0.3,label='Georgiev 2016')
# ax1.plot(galmass_s,nscmass_s,'o',markersize=8,color='red',alpha=0.3,label='Spengler 2017')
# ax1.plot(sj['logmstargal'][sj['T']<0],sj['logmnuc'][sj['T']<0],'^',alpha=0.5,color='red',label='Sanchez-Jansen 2019')
# # ax1.errorbar(np.log10(galmass_my[ttype<0]),np.log10(nscmass_my[ttype<0]),fmt='o',xerr=[err,err],yerr=[err,err],markersize=10,color='red',label='This work')
# # ax1.errorbar(np.log10(galmass_my[ttype>0]),np.log10(nscmass_my[ttype>0]),fmt='o',xerr=[err,err],yerr=[err,err],markersize=10,color='blue',label='This work')

# ax1.errorbar(np.log10(galmass_my[qual==0][ttype[qual==0]<0]),np.log10(nscmass_my[qual==0][ttype[qual==0]<0]),yerr=[err_nsc1,err_nsc2],fmt='o',markersize=8,color='red')
# ax1.errorbar(np.log10(galmass_my[qual==0][ttype[qual==0]>0]),np.log10(nscmass_my[qual==0][ttype[qual==0]>0]),yerr=[err_nsc3,err_nsc4],fmt='o',markersize=8,color='blue')
# #ax1.plot(np.log10(galmass_my[qual==1]),np.log10(nscmass_my[qual==1]),'o',markersize=10,mfc='none',color='black')
# ax1.plot(np.log10(galmass_my[qual==1]),np.log10(nscmass_my[qual==1]),'o',mfc='none',markersize=10,color='black')



# #------- plot galmass - Reff ----------#
# ax2.plot(np.log10(galmass_g),(reff_g),'o',markersize=8,alpha=0.3,color='blue',label='This work')
# ax2.plot(acsvcs['M_star'],acsvcs['r_hz']*80,'o',color='red',markersize=8,alpha=0.5)
# ax2.plot(np.log10(galmass_d[ttype_d<0]),reff_d[ttype_d<0],'*',markersize=16,color='red',alpha=0.5)
# ax2.plot(np.log10(galmass_d[ttype_d>0]),reff_d[ttype_d>0],'*',markersize=16,color='blue',alpha=0.5)

# ax2.errorbar(np.log10(galmass_my[qual==0][ttype[qual==0]<0]),(reff_my[qual==0][ttype[qual==0]<0]),fmt='o',yerr=[err1,err1],markersize=8,color='red',label='This work')
# ax2.errorbar(np.log10(galmass_my[qual==0][ttype[qual==0]>0]),(reff_my[qual==0][ttype[qual==0]>0]),fmt='o',yerr=[err2,err2],markersize=8,color='blue',label='This work')

# ax2.plot(np.log10(galmass_my[qual==1]),(reff_my[qual==1]),'o',markersize=8,mfc='none',color='black')

# # ax2.errorbar(np.log10(galmass_my[qual==0]),(reff_my[qual==0]),fmt='o',yerr=[err1,err1],markersize=10,color='black',label='This work')
# # ax2.plot(np.log10(galmass_my[qual==1]),(reff_my[qual==1]),'o',markersize=10,mfc='none',color='black')


# # --------------- plot nscmass - Reff ------------------#
# ax.errorbar(np.log10(nscmass_my[qual==0][ttype[qual==0]<0]),(reff_my[qual==0][ttype[qual==0]<0]),fmt='o',yerr=[err1,err1],color='red',markersize=8)
# ax.errorbar(np.log10(nscmass_my[qual==0][ttype[qual==0]>0]),(reff_my[qual==0][ttype[qual==0]>0]),fmt='o',yerr=[err2,err2],color='blue',markersize=8)
# ax.plot(np.log10(nscmass_my[qual==1]),(reff_my[qual==1]),'o',markersize=8,mfc='none',color='black')

# ax.plot(np.log10(nscmass_d[ttype_d<0]),reff_d[ttype_d<0],'*',markersize=16,color='red',alpha=0.5)
# ax.plot(np.log10(nscmass_d[ttype_d>0]),reff_d[ttype_d>0],'*',markersize=16,color='blue',alpha=0.5)

# ax.plot(np.log10(nscmass_g),(reff_g),'o',color='blue',alpha=0.3,markersize=8)

# ax.plot(acsvcs['M_nuc'],acsvcs['r_hz']*80,'o',color='red',markersize=8,alpha=0.5)
# galmass_my1 = galmass_my[qual==0]
# nscmass_my1 = nscmass_my[qual==0]
# reff_my1 = reff_my[qual==0]



# # ---------- fit galmass - nscmass relation ----------#
# x = np.log10(galmass_my1[~np.isnan(nscmass_my1)])   #galmass 
# y = np.log10(nscmass_my1[~np.isnan(nscmass_my1)])	#nscmass
# z1 = np.log10(reff_my1[~np.isnan(nscmass_my1)])		#reff
# z = np.log10(reff_my[~np.isnan(reff_my)])           #reff
# x1 = np.log10(galmass_my[~np.isnan(reff_my)])		
# lm = linmix.LinMix(x, y, xsig=np.repeat(0.,len(x)), ysig=np.repeat(0.,len(y)),K=2)   # galmass - nscmass 
# lm.run_mcmc(silent=True)
# lm1 = linmix.LinMix(x1, z, xsig=np.repeat(0.113,len(x1)), ysig=np.repeat(0.07,len(x1)),K=2)  # galmass - reff
# lm1.run_mcmc(silent=True)
# lm2 = linmix.LinMix(y, z1, xsig=np.repeat(0.07,len(y)), ysig=np.repeat(0.04,len(z1)),K=2)   #nscmass - reff
# lm2.run_mcmc(silent=True)

# print('gal-nsc')
# print("{}, {}".format(lm.chain['alpha'].mean(), lm.chain['alpha'].std()))
# print("{}, {}".format(lm.chain['beta'].mean(), lm.chain['beta'].std()))
# print("{}, {}".format(lm.chain['sigsqr'].mean(), lm.chain['sigsqr'].std()))

# print('re-gal')
# print("{}, {}".format(lm1.chain['alpha'].mean(), lm1.chain['alpha'].std()))
# print("{}, {}".format(lm1.chain['beta'].mean(), lm1.chain['beta'].std()))
# print("{}, {}".format(lm1.chain['sigsqr'].mean(), lm1.chain['sigsqr'].std()))

# print('re-nsc')
# print("{}, {}".format(lm2.chain['alpha'].mean(), lm2.chain['alpha'].std()))
# print("{}, {}".format(lm2.chain['beta'].mean(), lm2.chain['beta'].std()))
# print("{}, {}".format(lm2.chain['sigsqr'].mean(), lm2.chain['sigsqr'].std()))

# # --------- galmass - nscmass ----------#
# print('gal-nsc')
# sx = np.repeat(0.113,len(x))
# sy = np.repeat(0.068,len(x))
# linear_model = Model(linear_func)
# data = RealData(x,y,sx=sx,sy=sy)
# odr = ODR(data,linear_model,beta0=[0.,1.])
# out = odr.run()
# out.pprint()
# be = out.beta
# popt = curve_fit(fit_line,np.log10(galmass_my[~np.isnan(nscmass_my)]),np.log10(nscmass_my[~np.isnan(nscmass_my)]))[0]

# # -------------galmass - reff -------------#
# print('gal-re')
# sx = np.repeat(0.113,len(x1))
# sy = np.repeat(-1,len(x1))
# linear_model = Model(linear_func)
# data = RealData(x1,z,sx=sx,sy=sy)
# odr = ODR(data,linear_model,beta0=[0.,1.])
# out = odr.run()
# out.pprint()
# be1 = out.beta
# popt1 = curve_fit(fit_line,np.log10(galmass_my),np.log10(reff_my))[0]

# # ---------------- nscmass - reff ---------#
# print('nsc-re')
# sx = np.repeat(0.068,len(x1))
# sy = np.repeat(-1,len(x1))
# linear_model = Model(linear_func)
# data = RealData(y,z1,sx=sx,sy=sy)
# odr = ODR(data,linear_model,beta0=[0.,1.])
# out = odr.run()
# out.pprint()
# be2 = out.beta

# popt2 = curve_fit(fit_line,np.log10(nscmass_my[~np.isnan(nscmass_my)]),np.log10(reff_my[~np.isnan(nscmass_my)]))[0]
# x = np.linspace(6.5,12)
# x1 = np.linspace(4,9)

# logre_l = 0.356*x - 0.356*np.log10(5.61e9)-0.012+np.log10(3.44) #late
# logre_e = 0.326*x - 0.326*np.log10(2.09e9)-0.011+np.log10(6.11) #early
# loggal_l = 1.001*x - 1.001*np.log10(3.94e9)+0.016+np.log10(2.78e6)
# loggal_e = 1.363*x - 1.363*np.log10(1.75e9)+0.01+np.log10(2.24e6)
# lognsc_l = 0.321*x1 - 0.321*np.log10(3.6e6)-0.011+np.log10(3.31)
# lognsc_e = 0.347*x1 - 0.347*np.log10(1.95e6)-0.024 +np.log10(6.27)

# # -------- ax1 - galmass -nscmass --------#
# ax1.set_xlabel(r'log (Galaxy stellar Mass)$[M_\odot]$')
# ax1.set_ylabel('log (NSC Mass)$[M_\odot]$')
# ax1.plot(x,loggal_l,color='black',linestyle=':',label='Late-Type')
# #ax1.plot(x,loggal_e,color='black',linestyle='dashed',label='Early-Type')
# #ax1.plot(x,fit_line(x,lm.chain['beta'].mean(),lm.chain['alpha'].mean()),color='black',label='This work')
# ax1.plot(x,fit_line(x,be[0],be[1]),color='black',label='This work')
# #ax1.plot(x,fit_line(x,popt[0],popt[1]),color='gray',label='This work')

# #ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

# ax1.set_xlim(7,11.5)
# ax1.legend(fontsize=15,loc=2,frameon=False)

# # ------- ax2 - galmass - Reff --------------#
# ax2.set_xlabel(r'log (Galaxy stellar Mass)$[M_\odot]$')
# ax2.set_ylabel('R$_{eff}$ (parsecs)')
# ax2.plot(x,10**logre_l,color='black',linestyle=':',label='Late-Type')
# ax2.plot(x,10**(be1[1]+be1[0]*x),color='black',label='This work')
# #ax2.plot(x,10**(popt1[1]+popt1[0]*x),color='gray',label='This work')

# ax2.set_yscale('log')
# ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax2.set_xlim(7,11.5)

# #----------ax - nscmass - Reff --------------------#
# ax.set_xlabel(r'log (NSC Mass)$[M_\odot]$')
# ax.set_ylabel('R$_{eff}$ (parsecs)')
# ax.plot(x1,10**lognsc_l,color='black',linestyle=':',label='Late-Type')
# ax.plot(x1,10**(be2[1]+be2[0]*x1),color='black',label='This work')
# #ax.plot(x1,10**(popt2[1]+popt2[0]*x1),color='gray',label='This work')

# ax.set_ylim(0.15,99)
# ax.set_yscale('log')
# ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

# plt.show()

#!!!!!!!!!----------------------------------------  imfit-images -------------------------------------------#
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':25,'weight':'normal'})
rc('axes',**{'linewidth':3})
plt.rcParams.update({'mathtext.default':'regular'})

model = fits.open('../Data/NGC3344/model_image_ser_sersic.fits')
model = model[0].data
img = fits.open('../Data/NGC3344/WFC3_UVIS/F814W/ngc3344_f814w_drc.fits')
img = img[1].data
img1 = gaussian_filter(img,0.7)
res = fits.open('../Data/NGC3344/res_image_ser_sersic.fits')
res = res[0].data
fig,(ax1,ax2,ax3) = plt.subplots(figsize=(15,5),nrows=1,ncols=3,sharey=True)
levels = np.array([0.269301, 0.271204, 0.273893, 0.27769 , 0.283054, 0.290631,
       0.301334, 0.316452, 0.337807, 0.367971, 0.410579, 0.470765,
       0.55578 , 0.675867, 0.84549 , 1.0851  , 1.42355 , 1.90162 ,
       2.57692 , 3.5308  , 4.87819 ])
levels = np.array([0,0.176317,0.441024,0.83843,1.43506,2.33079,3.67556,5.69447,8.72549,
	13.276,20.1077,30.3642,45.7624,68.88,103.587,155.692,233.918,351.36])

levels1 = np.array([1.88273e-03, 4.27183e-03, 7.30352e-03, 1.11506e-02, 1.60324e-02,
       2.22272e-02, 3.00882e-02, 4.00635e-02, 5.27217e-02, 6.87846e-02,
       8.91677e-02, 1.15033e-01, 1.47855e-01, 1.89505e-01, 2.42357e-01,
       3.09425e-01, 3.94531e-01, 5.02527e-01, 6.39570e-01, 8.13472e-01,
       1.03415e+00, 1.31417e+00, 1.66952e+00, 2.12044e+00, 2.69263e+00,
       3.41873e+00, 4.34012e+00, 5.50932e+00, 6.99300e+00])
levels1 = np.array([0,0.176317,0.441024,0.83843,1.43506,2.33079,3.67556,5.69447,8.72549,
	13.276,20.1077,30.3642,45.7624,68.88,103.587,155.692,233.918,351.36])
# model = -2.5*np.log10(model)+24.684+5*np.log10(0.04)
# img = -2.5*np.log10(img)+24.684+5*np.log10(0.04)
# img1 = -2.5*np.log10(img1)+24.684+5*np.log10(0.04)
ax2.contour(model,colors='limegreen',levels=levels,origin='lower',linewidths=1,extent=[-10,10,-10,10])
im = ax2.imshow(model,cmap='Greys_r',vmin=-0.9,vmax=20,origin='lower',extent=[-10,10,-10,10])
#ax2.colorbar()
fig.colorbar(im,ax=ax2,ticks=[0,10,20],label='Counts')
#cbar.ax2.set_yticklabels(['0','10','20'])
#ax3.contour(res[50:150,50:150],colors='limegreen',levels=levels,origin='lower',linewidths=0.8,extent=[-5,5,-5,5])
im1 = ax3.imshow(res/model,cmap='Greys_r',vmin=-0.2,vmax=0.2,origin='lower',extent=[-10,10,-10,10])
fig.colorbar(im1,ax=ax3,ticks=[-0.2,0,0.2],label ='%')
ax1.contour(img1[2474:2674,2039:2239],colors='limegreen',levels=levels1,origin='lower',linewidths=1,extent=[-10,10,-10,10])
im2 = ax1.imshow(img[2474:2674,2039:2239],cmap='Greys_r',vmin=-.9,vmax=20,origin='lower',extent=[-10,10,-10,10])
fig.colorbar(im2,ax=ax1,ticks=[0,10,20],label='Counts')
ax1.set_title('Data')
ax2.set_title('Model ')
ax3.set_title('Residual')
ax1.set_xlabel('Radius (")')
ax2.set_xlabel('Radius (")')
ax3.set_xlabel('Radius (")')
ax1.set_ylabel('Radius (")')
plt.show()

# !!!!!!!!!----------------------------------------- sersic index  -----------------------------------#

# a = np.loadtxt('../datafiles/table_modifiedsurf.csv',dtype=str,delimiter=',',skiprows=2)
# c = np.loadtxt('../datafiles/list.csv',dtype=str,skiprows = 1,delimiter=',')
# d = Table.read('../datafiles/table_sersicparams.csv')
# n_my = a[:,14].astype(float)
# galmass_my = a[:,8].astype(float)
# nscmass_my = a[:,12].astype(float)
# qual = a[:,-2].astype(int)

# galmass_d = c[:,12].astype(float)
# galmass_d = galmass_d[29:]
# n_d = c[:,10].astype(float)
# n_d = n_d[29:]
# nscmass_d = c[:,-5].astype(float)
# nscmass_d = nscmass_d[29:]

# fig,ax = plt.subplots()
# ax.plot((nscmass_d),(n_d),'*',markersize=13,color='gray',label='Carson 2015 + Ngyuen 2018')
# ye_lo = d['sersic_err_lo'][qual==0]
# ye_hi = d['sersic_err_hi'][qual==0]

# ax.errorbar((nscmass_my[qual==0]),(n_my[qual==0]),fmt='o',yerr=[ye_lo,ye_hi],markersize=8,capsize=2,color='black')
# ax.plot((nscmass_my[qual==0]),(n_my[qual==0]),'o',markersize=8,color='black',label='This work')

# ax.plot((nscmass_my[qual==1]),(n_my[qual==1]),'o',markersize=10,mfc='none',color='k')
# ax.axhline(6,linestyle='dashed',color='gray')
# xdata = np.linspace(0.1,15,100)
# #ax.plot(10**xdata,10**((xdata/2.69)-2.4262))
# #ax.plot(xdata,10**(2.69*np.log10(xdata/3)+7.81))
# #ax.plot(xdata,10**(3.7*np.log10(xdata/3)+7.98)-3.1*(np.log10(xdata/3))**2)

# ax.set_xscale('log')
# #galre_g = curve_fit(fit_line,np.log10(galmass_g),np.log10(reff_g))[0]
# ax.set_yscale('log')
# ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax.legend(fontsize=15,loc=3,frameon=False)
# 	#,bbox_to_anchor=(0.0, 0.2, 0.5, 0.5))
# ax.set_xlabel('NSC Mass$(M_\odot)$')
# ax.set_ylabel('Se$\'$rsic Index')

# plt.show()

# ----------------------- PA- ell plot -----------------------------#

# a = Table.read('../datafiles/re1.csv')
# ell_g = a['ell_gal']
# ell_n = a['ell']
# pa_diff = a['PA_diff']
# pa_diff = pa_diff[ell_n>0.2]
# pa_g = a['PA_gal']
# pa_g = pa_g[ell_g>0.15]
# pa_n = a['PA_nuc']
# x = np.arange(0,1,0.05)
# fig,(ax1,ax2) = plt.subplots(ncols=2, nrows=1)
# ax1.plot(ell_n,ell_g,'o',color='black',markersize=8)
# ax1.plot(x,x,color='gray')
# ax2.hist(pa_diff[~np.isnan(pa_diff)])
# ax1.set_xlim(0,0.6)
# ax1.set_ylim(0,0.9)
# ax1.set_xlabel('NSC ellipticity')
# ax2.set_ylabel('#')
# ax1.set_ylabel('Galaxy ellipticity')
# ax2.set_xlabel('$\delta$ PA (NSC - Galaxy)')
# plt.show()



