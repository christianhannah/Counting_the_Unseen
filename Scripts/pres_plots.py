import numpy as np
from matplotlib import pyplot as plt


a = np.loadtxt('./Data_Sets/Pechetti20_datafiles/scaling_relation_table.dat',dtype=str,skiprows=1)

b = np.loadtxt('./Data_Sets/Pechetti20_datafiles/Table1.csv',delimiter=',',dtype=str,skiprows=1)

c = np.loadtxt('./Data_Sets/Pechetti20_datafiles/list.csv',dtype=str,skiprows = 1,delimiter=',')

d = np.loadtxt('./Data_Sets/Pechetti20_datafiles/disp.csv',dtype=str,delimiter=',',skiprows=1)

galmass = a[:,6].astype(float)
ncmass = a[:,7].astype(float)
reff = a[:,-2].astype(float)
ttype = a[:,1].astype(float)
dist = a[:,2].astype(float)
reffpc = reff*np.pi*dist/0.648
galmass_spec = a[:,6].astype(float)

galaxy_my = b[:,0]
ncmass_my = b[:,-4].astype(float)*b[:,-5].astype(float)
reff_my = b[:,-1].astype(float)
galmass_my = b[:,-6].astype(float)                     
n_my = b[:,-2].astype(float)
ttype_my=b[:,2].astype(int)

galmass_d = c[:,12].astype(float)
galmass_d = galmass_d[29:42]
n_d = c[:,9].astype(float)
n_d = n_d[29:42]
reff_d = c[:,11].astype(float)
reff_d = reff_d[29:42]
ncmass_d = c[:,-5].astype(float)
ncmass_d = ncmass_d[29:42]

lum = d[:,-4].astype(float)
zpt = d[:,-7].astype(float)
scale = d[:,-8].astype(float)
msun = d[:,-6].astype(float)
mu = -2.5*np.log10(lum)+zpt+5*np.log10(scale)
sb = 10**(-0.4*(mu-msun-21.5723))

c1,c2,alpha,beta = 3.44,5.61*10**9,0.356, -0.012#late
c1e,c2e,alphae,betae=6.11,2.09*10**9,0.326,-0.011  #early

x = np.linspace(10**7,10**13,10000)
re = c1*10**(alpha*np.log10(x/c2)+beta)

x1 = np.linspace(10**7,10**13,10000)
re1 = c1e*10**(alphae*np.log10(x1/c2e)+betae)

ncmass_type = a[:,-5]
ind = np.where(ncmass_type=='Walcher05')[0]

#------------ Reff - galmass --------------#
plt.figure()
plt.plot(galmass[ttype<0],reffpc[ttype<0],'o',alpha=0.3,color='r',markersize=8)
plt.plot(galmass[ttype>0],reffpc[ttype>0],'o',alpha=0.3,color='b',markersize=8)
plt.plot((galmass_my)[ttype_my<0],reff_my[ttype_my<0],'*',color='r',markersize=25)
plt.plot((galmass_my)[ttype_my>0],reff_my[ttype_my>0],'*',color='b',markersize=25)
#plt.plot(x,re,linestyle='dashed',color='b')
#plt.plot(x1,re1,linestyle='dashed',color='b')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Galaxy Mass (M$_\odot$)')
plt.ylabel('')
plt.show()

#--------------- n - ncmass -----------------#

plt.figure()
plt.plot(ncmass_d,(n_d),'*',markersize=15,label = 'Nguyen18 + Carson16')
plt.plot(ncmass_my,(n_my),'*',markersize=15, label = 'This work')
plt.xscale('log')
plt.yscale('log')
plt.show()

#----------------  sb -central pix -----------------#
plt.figure()
plt.plot(galmass_my,sb,'o')
plt.xscale('log')
plt.show()

#-------------------- ncmass - galmass ---------------------------#
plt.figure()
plt.plot(galmass,ncmass,'o',color='c',alpha=0.3,markersize=8,label='Photometric')
plt.plot(galmass_my,np.log10(ncmass_my),'*',color='orange',markersize=25,label='This Work')
plt.plot(galmass[ind],ncmass[ind],'*',color='green',markersize=25,label='Existing measurements')
#plt.plot(galmass_d,np.log10(ncmass_d),'*',color='green',markersize=25)
plt.xscale('log')
plt.xlabel('Galaxy Mass (M$_\odot$)')
plt.ylabel('log(NSC Mass) (M$_\odot$)')



