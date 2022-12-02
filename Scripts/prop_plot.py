import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
rc('axes',**{'linewidth':3})
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
a = np.loadtxt('../datafiles/prop_bh.dat',dtype=str)   #amy reines
b = np.loadtxt('../datafiles/prop_bh1.dat',dtype=str)  #jenny greene
c = np.loadtxt('../datafiles/prop_bh2.dat',dtype=str)  #saglia
e = np.loadtxt('../datafiles/prop_bh3.dat',dtype=str)  #amy reines dynamical bh

d = np.loadtxt('../datafiles/saglia2016_mbh_mbul_sigma_scaling.txt',dtype=str,skiprows=2)
x = np.arange(6.5,12.5,0.05)
logBH_line = 0.962*x-2.099 #saglia relation MBH-Mbu
logBH_line = 1.4*x-6.44 #kormendy ho
logBH_line2 = 1.05*x- 4.1#reines relation Mbh-Mstar agn
# dieu_bh = 10**np.array([5.7,5.94,5,6.1,5.176])
# dieu_mstar = 10**np.array([9.3,9.69 ,8.95,8.89,8.92])
dieu_bh = 10**np.array([6.09,5.6,5.94,5.5,4.84,3.85])
dieu_mstar = 10**np.array([8.89,9,9.77,9.38,8.97,8.98])

plt.plot(dieu_mstar,dieu_bh,'*',color='b',markersize=12,label='BH masses derived by our team')
plt.plot(10**e[:,-3].astype(float),10**e[:,-2].astype(float),'o',color='k',alpha=0.7,label='Dynamical BH mass estimates')

plt.plot(10**x,10**logBH_line,color='k',alpha=0.7, label='M$_\mathrm{\star}$ - M$_{BH}$ relation')
plt.plot(10**x,10**logBH_line2,color='k',linestyle='dashed',label='Relation derived for low-mass galaxies')


#plt.plot(10**b[:,8][0:11].astype(float), 10**b[:,4][0:11].astype(float),'o',markerfacecolor='None',color='k',label='Greene 2016')
#plt.plot(10**c[:,18].astype(float), 10**c[:,12].astype(float),'o',color='k',label='Saglia 2016')
#plt.plot(10**d[:,8].astype(float), 10**d[:,4].astype(float),'o',color='k',label='Saglia 2016')
#plt.plot(4*10**8,3*10**4,'*',color='r',markersize=20)

plt.xlabel(r'(M$_\mathrm{\star}$/M$_\odot$)')
plt.ylabel(r'(M$_\mathrm{BH}$ /M$_\odot$)')

plt.annotate('ESO 274-01',xy=(8.6,4.28), color='r')
plt.legend(fontsize=12)
plt.xlim(10**7,10**12.5)
plt.yscale('log')
plt.xscale('log')
plt.show()