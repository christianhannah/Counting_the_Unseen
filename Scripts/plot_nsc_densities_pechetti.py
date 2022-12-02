import numpy as np
from matplotlib import pyplot as plt
import os
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from matplotlib import rc
import pdb
# =============================================================================
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':26,'weight':'normal'})
# rc('axes',**{'linewidth':4})
# plt.rcParams.update({'mathtext.default':'regular'})
# =============================================================================
#plt.close()

params = np.loadtxt('./Data_Sets/Pechetti20_datafiles/Table1.csv',delimiter=',',dtype=str,skiprows=1)  # galaxy,disp,disperr,vel,seeing,distance,scale,zpt,msun,extinctioin

params1 = np.loadtxt('./Data_Sets/Pechetti20_datafiles/disp.csv',delimiter=',',dtype=str,skiprows=1)

params2 = np.loadtxt('./Data_Sets/Pechetti20_datafiles/acs_wfc_data.csv',skiprows=1,dtype=str,delimiter=',')

outer_re = np.loadtxt('./Data_Sets/Pechetti20_datafiles/eqrad.txt',skiprows=1,dtype=str)
galmass = params[:,7].astype(float)      
print(galmass)
f = lambda x: sorted(range(len(x)), key=lambda i: x[i])
idx = f(galmass)
outer = outer_re[:,1].astype(float)[idx]
galmass2 = galmass[idx]
gal = params[:,0][idx]
ml  = params[:,-6].astype(float)[idx]
D = params1[:,5].astype(float)[idx]
lum = params1[:,-6].astype(float)[idx]
zpt = params1[:,7].astype(float)[idx]
scale = params1[:,6].astype(float)[idx]
msun = params1[:,8].astype(float)[idx]
mu = -2.5*np.log10(lum)+zpt+5*np.log10(scale)
sb = 10**(-0.4*(mu-msun-21.5723))
nsc_ml_col = params[:,-1].astype(float)[idx]
colors = cm.rainbow(np.linspace(0,1,len(params)))
dir1 = './Data_Sets/Pechetti20_datafiles/Data/'
plt.figure()
sc,n,r2,med_density,galmass1,gal1 = [],[],[],[],[],[]
fig,(ax,ax1) = plt.subplots(ncols=2,nrows=1,figsize=(7,5))
#fig,ax = plt.subplots()
a = np.loadtxt('./Data_Sets/Pechetti20_datafiles/fitparams_king.dat',dtype=str)

j0 = a[:,17].astype(float)
rc = a[:,5].astype(float)

j0 =10**j0
rc = 10**rc
r1 = np.linspace(1,30,1000)
# G = 6.67*10**-11
# c = psi/sig**2
# c = 12
# y = -8*np.pi*G*rho*r*((np.exp(c)*erf(np.sqrt(c)))-((np.sqrt(4*c/np.pi**2))*(1+(2*c/3))))

#%%

# rho = ((np.exp(c)*erf(np.sqrt(c)))-((np.sqrt(4*c/np.pi**2))*(1+(2*c/3))))
for i in range(0,153):

	j = j0[i*3]/((1+(r1/rc[i*3])**2)**1.5)

	ax1.plot(np.log10(r1),np.log10(j),color='grey',alpha=1)
#plt.xscale('log')

#plt.show()
pdb.set_trace()

#%%

for i in range(0,len(params)):
	
    pc = D[i]*np.pi/0.648
    radius = scale[i]*pc
    outer_pc = outer[i]*pc
    path = dir1+gal[i]+'/mge_profile_norm0.txt'
	
    try:
        if (os.stat(path).st_size==0):
            continue
    except FileNotFoundError:
        print (gal[i])
        print('not found')
        continue
	# print (i)
	# if (i == 3): #4
	# 	continue
	# elif (i==9):
	# 	continue
	# elif (i==16):
	# 	continue
	# elif (i==19): #25
	# 	continue
	# elif (i==27):
	# 	continue

    data = np.loadtxt(path)
    mass = data[:,0]*nsc_ml_col[i]
    sigma = data[:,1]*pc
    q = data[:,2]
    r = np.linspace(radius,30,1000)
    r1 = np.linspace(radius,30,1000)
    rho,rho1 = np.zeros(len(r)),np.zeros(len(r1))
    print(gal[i])
    inc = 90*np.pi/180
    q1 = (q**2 - (np.cos(inc)**2))/(np.sin(inc)**2)
	
    for k in range(0,len(r)):
        for j in range(0,len(data)):
            rho[k]+= mass[j]*np.exp(-(r[k]**2)/(2*sigma[j]**2))/(q[j]*(np.sqrt(2*np.pi)*sigma[j])**3)
            rho1[k]+= mass[j]*np.exp(-(r1[k]**2)/(2*sigma[j]**2))/(q1[j]*(np.sqrt(2*np.pi)*sigma[j])**3)
	
    sc.append(np.log10(rho[0]))
    print(np.log10(rho[0]))
    print(nsc_ml_col[i])
    n.append(i)
	# r1 = []
	# x = np.linspace(1,30,100)
	# y = np.linspace(1,30,100)
	# sig = np.zeros((len(x),len(y)))
	# for k in range(0,len(x)):
	# 	for l in range(0,len(y)):
	# 		for j in range(0,len(data)):
	# 			sig[k][l]+= mass[j]*np.exp(-(x[k]**2+(y[l]/q[j]**2))/(2*sigma[j]**2))/(q[j]**2*2*np.pi*sigma[j]**2)
	# 			r1.append(np.sqrt(x[k]**2+y[l]**2))
	# r2.append(np.log10(sig).max())
	#print(np.log10(rho[0]))
	# if (np.isnan(outer_pc)):
	# 	ax.plot(np.log10(r),np.log10(rho),color=colors[i],linewidth=4)
	# else:
	#ax.plot(np.log10(r),np.log10(rho),color=colors[i],linewidth=4)
    ax.plot(np.log10(r[r<outer_pc]),np.log10(rho[r<outer_pc]),color=colors[i],linewidth=4)
    ax1.plot(np.log10(r1[r<outer_pc]),np.log10(rho1[r<outer_pc]),color=colors[i],alpha=0.8,linewidth=4)
    if (np.isnan(rho[0])):
        continue
    else:
        f = interp1d(np.log10(r),np.log10(rho))
        med_density.append(f(0.7))
        #med_density.append(np.log10(rho[0]))
        galmass1.append(galmass2[i])
        gal1.append(gal[i])
	#plt.plot(np.log10(r1),np.log10(sig),color=colors[i],linewidth=2)

	# plt.colorbar(sc)

# plt.figure()
# plt.plot(galmass,(sb*ml/(np.pi*scale**2*pc**2)),'*',markersize=20,color='orange')
# plt.xscale('log')
# plt.ylabel(r'log (Central luminosity) (L$_\odot$/pc$^2$)')
# plt.xlabel('Galaxy mass (M$_\odot$)')
# plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
# color = sc.get_facecolor()[idx]
# color1 = sc.get_facecolor()[idx1]
# plt.plot(distance[idx][0],np.log10(mass[idx][0]),'*',color='black',markersize=18,markeredgewidth=2,mfc='none',label='Masses from Cook (2014)')
# plt.plot(distance[idx1][0],np.log10(mass[idx1][0]),'o',color='black',markersize=13,markeredgewidth=2,mfc='none',label='Masses from color-M/L relations')

# for j in range(0,len(color)):

# 	plt.plot(distance[idx][j],np.log10(mass[idx][j]),marker='*',markeredgecolor='black',markeredgewidth=1,markersize=23,color =color[j])
# for j in range(0,len(color1)):

# 	plt.plot(distance[idx1][j],np.log10(mass[idx1][j]),marker='o',markeredgecolor='black',markeredgewidth=1,markersize=16,color =color1[j])

# ax.errorbar(distance,np.log10(mass),yerr=[err1,err],fmt=None,marker=None,ecolor='black',elinewidth = 2,capsize=2,capthick=2,zorder=0)
# ax.set_ylabel(r'log (Galaxy stellar Mass)$[M_\odot]$')
# ax.set_xlabel('Distance (Mpc)')
# cbar = plt.colorbar(sc)
# cbar.set_label('T-type',fontsize=30)
# plt.text(12,7.6,'Early-Type',fontsize=30)
# plt.text(12,11.13,'Late-Type',fontsize=30)
# plt.legend(fontsize=20)
#%%
#---- density plot -----------#
x = np.linspace(-0.1,1.5,29)
y = np.linspace(-4,6,29)
#sc = ax.scatter(x,y,s=15,c=np.log10(np.sort(galmass)),cmap='rainbow')
cbar = plt.colorbar(sc)
cbar.set_label('log (Galaxy Mass) (M$_\odot$)',fontsize=20,rotation=270,labelpad=25)
ax1.set_xlabel('log (Radius) [pc]')
ax.set_ylabel('log (Density [M$_\odot$/pc$^3$]')
ax.set_xlabel('log (Radius) [pc]')
#ax1.set_ylim(-3,6.2)
#%%
#---------------median density plot -----------#
plt.figure()
plt.plot(galmass1,med_density,'o',markersize=10,color='black')
plt.xlabel('Galaxy mass (M$_\odot$)')
plt.ylabel('log (Density [M$_\odot$/pc$^3$]')
plt.xscale('log')

plt.show()
