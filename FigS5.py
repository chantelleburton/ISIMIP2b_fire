
#export PATH=/project/avd/sss/temp_envs/miniconda3/bin:$PATH
#export LD_LIBRARY_PATH=/usr/lib64/atlas:$LD_LIBRARY_PATH
#. activate default_legacy-2018_05_03

import os
import glob
import sys
sys.path.append('/net/home/h03/kwilliam/other_fcm/jules_py/trunk/jules/')
#~kwilliam/other_fcm/jules_py/trunk/jules/jules.py
import jules
import iris

import matplotlib as mpl
import pylab
import matplotlib.pyplot as plt
import numpy as np
from iris.plot import pcolormesh
import iris.quickplot as qplt
import iris.analysis.cartography
import iris.plot as iplt
import h5py
from iris.coord_systems import GeogCS
import cartopy.crs as ccrs
import cf_units
from collections import OrderedDict 
import numpy.ma as ma
import iris.coord_categorisation
from mpl_toolkits.axes_grid1 import make_axes_locatable

'''
var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
cube = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/GFDL/2000-2010/*2000*.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN) 
cube = cube[0] 


ax = plt.subplot(10,4,1)
mesh1=iplt.pcolormesh(cube)
plt.text(-0.3, 0.4, 'test1', transform=ax.transAxes)  

plt.subplot(10,4,2)
mesh2=iplt.pcolormesh(cube)

plt.subplot(10,4,3)
plt.subplot(10,4,4)
plt.subplot(10,4,5)
plt.subplot(10,4,6)
plt.subplot(10,4,7)
plt.subplot(10,4,8)
plt.subplot(10,4,9)
plt.subplot(10,4,10)
plt.subplot(10,4,11)
plt.subplot(10,4,12)
plt.subplot(10,4,13)
plt.subplot(10,4,14)
plt.subplot(10,4,15)
plt.subplot(10,4,16)
plt.subplot(10,4,17)
plt.subplot(10,4,18)
plt.subplot(10,4,19)
plt.subplot(10,4,20)
plt.subplot(10,4,21)
plt.subplot(10,4,22)
plt.subplot(10,4,23)
plt.subplot(10,4,24)
plt.subplot(10,4,25)
plt.subplot(10,4,26)
plt.subplot(10,4,27)
plt.subplot(10,4,28)
plt.subplot(10,4,29)
plt.subplot(10,4,30)
plt.subplot(10,4,31)
plt.subplot(10,4,32)
plt.subplot(10,4,33)
plt.subplot(10,4,34)
plt.subplot(10,4,35)
plt.subplot(10,4,36)
plt.subplot(10,4,37)
plt.subplot(10,4,38)
plt.subplot(10,4,39)
plt.subplot(10,4,40)

#Shared colorbar
#plt.colorbar(orientation='horizontal', pad=0.08, )
                      #Position: Left, Bottom, width, height
#colorbar_axes = plt.gcf().add_axes([0.15, 0.14, 0.45, 0.03])
#colorbar = plt.colorbar(mesh1, location='right', orientation='vertical') #, pad=0.08, )

      #Position: Right,Bottom, width, height
cax = plt.axes([0.85, 0.5, 0.05, 0.35])
plt.colorbar(cax=cax)

      #Position: Right,Bottom, width, height
cax = plt.axes([0.85, 0.05, 0.05, 0.25])
plt.colorbar(cax=cax)

#cax = plt.axes([0.85, 0.1, 0.075, 0.4])
#plt.colorbar(cax=cax)

#      L to R, Top to Bottom    
#plt.text(-0.1, 0.4, '1.5', transform=ax.transAxes, fontsize=14)  


plt.text(-0.3, 0.4, 'test1', transform=ax.transAxes)  

plt.text(-0.3, 3.2, 'test2')
plt.text(-0.3, 2.9, 'test3')
plt.text(-0.3, 2.6, 'test4')
plt.text(-0.3, 2.3, 'test5')
plt.text(-0.3, 2.0, 'test6')
plt.text(-0.3, 1.9, 'test7')
plt.text(-0.3, 1.6, 'test8')
plt.text(-1.7, 1.35, 'test9')
plt.text(-1.6, 1.0, 'test10')
plt.text(-1.5, 0.6, 'test11')
plt.text(-2.4, 0.25, 'test12')
plt.text(-1.4, 0.1, 'test13')



plt.show()

raise Exception('stop')
'''




### Get  CCI data
var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
CCI = iris.load_cube('/project/LandCoverCCI/V1/ancils/n96eFracCover_v1.6.1_9PFT_2010/LC_CCI/vegfrac_lc-cci_9PFTs_2010_v1.6.1.nc')
cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/GFDL/2010-2019/PFTYears.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN)
cs_new = iris.coord_systems.GeogCS(6371229.)
CCI.coord('latitude').coord_system = cs_new
CCI.coord('longitude').coord_system = cs_new
cube.coord('latitude').coord_system = cs_new
cube.coord('longitude').coord_system = cs_new
CCI = CCI.regrid(cube, iris.analysis.Linear()) 

#raise Exception('stop')
CCI_BETr = CCI[0,:,:]
CCI_BETe = CCI[1,:,:]
CCI_BLD = CCI[2,:,:]
CCI_NLE = CCI[3,:,:]
CCI_NLD = CCI[4,:,:]
CCI_C3 = CCI[5,:,:]
CCI_C4 = CCI[6,:,:]
CCI_SBE = CCI[7,:,:]
CCI_SBD = CCI[8,:,:]
CCI_BS = CCI[11,:,:]

#Present Day
models = ('GFDL', 'IPSL', 'MIROC','HadGEM2')
PD = {}
for model in models:
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
    cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/'+model+'/2010-2019/PFTYears.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN) 
    PD['BETr'+model] = cube[1,:,:]
    PD['BETe'+model] = cube[2,:,:]
    PD['BLD'+model] = cube[0,:,:]
    PD['NLE'+model] = cube[4,:,:]
    PD['NLD'+model] = cube[3,:,:]
    PD['C3'+model] = cube[5,:,:]
    PD['C4'+model] = cube[8,:,:]
    PD['SBE'+model] = cube[12,:,:]
    PD['SBD'+model] = cube[11,:,:]
    PD['BS'+model] = cube[15,:,:]

PD_BETr = (PD['BETrGFDL']+PD['BETrIPSL']+PD['BETrMIROC']+PD['BETrHadGEM2'])/4
PD_BETe = (PD['BETeGFDL']+PD['BETeIPSL']+PD['BETeMIROC']+PD['BETeHadGEM2'])/4
PD_BLD = (PD['BLDGFDL']+PD['BLDIPSL']+PD['BLDMIROC']+PD['BLDHadGEM2'])/4
PD_NLE = (PD['NLEGFDL']+PD['NLEIPSL']+PD['NLEMIROC']+PD['NLEHadGEM2'])/4
PD_NLD = (PD['NLDGFDL']+PD['NLDIPSL']+PD['NLDMIROC']+PD['NLDHadGEM2'])/4
PD_C3 = (PD['C3GFDL']+PD['C3IPSL']+PD['C3MIROC']+PD['C3HadGEM2'])/4
PD_C4 = (PD['C4GFDL']+PD['C4IPSL']+PD['C4MIROC']+PD['C4HadGEM2'])/4
PD_SBE = (PD['SBEGFDL']+PD['SBEIPSL']+PD['SBEMIROC']+PD['SBEHadGEM2'])/4
PD_SBD = (PD['SBDGFDL']+PD['SBDIPSL']+PD['SBDMIROC']+PD['SBDHadGEM2'])/4
PD_BS = (PD['BSGFDL']+PD['BSIPSL']+PD['BSMIROC']+PD['BSHadGEM2'])/4

#1.5 deg
Fut15 = {}
for model in models:
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
    cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/'+model+'/RCP60_1p5/Years_PFTlayers.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN) 
    Fut15['BETr'+model] = cube[1,:,:]
    Fut15['BETe'+model] = cube[2,:,:]
    Fut15['BLD'+model] = cube[0,:,:]
    Fut15['NLE'+model] = cube[4,:,:]
    Fut15['NLD'+model] = cube[3,:,:]
    Fut15['C3'+model] = cube[5,:,:]
    Fut15['C4'+model] = cube[8,:,:]
    Fut15['SBE'+model] = cube[12,:,:]
    Fut15['SBD'+model] = cube[11,:,:]
    Fut15['BS'+model] = cube[15,:,:]

Fut15_BETr = (Fut15['BETrGFDL']+Fut15['BETrIPSL']+Fut15['BETrMIROC']+Fut15['BETrHadGEM2'])/4
Fut15_BETe = (Fut15['BETeGFDL']+Fut15['BETeIPSL']+Fut15['BETeMIROC']+Fut15['BETeHadGEM2'])/4
Fut15_BLD = (Fut15['BLDGFDL']+Fut15['BLDIPSL']+Fut15['BLDMIROC']+Fut15['BLDHadGEM2'])/4
Fut15_NLE = (Fut15['NLEGFDL']+Fut15['NLEIPSL']+Fut15['NLEMIROC']+Fut15['NLEHadGEM2'])/4
Fut15_NLD = (Fut15['NLDGFDL']+Fut15['NLDIPSL']+Fut15['NLDMIROC']+Fut15['NLDHadGEM2'])/4
Fut15_C3 = (Fut15['C3GFDL']+Fut15['C3IPSL']+Fut15['C3MIROC']+Fut15['C3HadGEM2'])/4
Fut15_C4 = (Fut15['C4GFDL']+Fut15['C4IPSL']+Fut15['C4MIROC']+Fut15['C4HadGEM2'])/4
Fut15_SBE = (Fut15['SBEGFDL']+Fut15['SBEIPSL']+Fut15['SBEMIROC']+Fut15['SBEHadGEM2'])/4
Fut15_SBD = (Fut15['SBDGFDL']+Fut15['SBDIPSL']+Fut15['SBDMIROC']+Fut15['SBDHadGEM2'])/4
Fut15_BS = (Fut15['BSGFDL']+Fut15['BSIPSL']+Fut15['BSMIROC']+Fut15['BSHadGEM2'])/4

#2.0 deg
models = ('GFDL', 'IPSL', 'MIROC','HadGEM2')
Fut20 = {}
for model in models:
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
    cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/'+model+'/RCP60_2p0/Years_PFTlayers.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN) 
    Fut20['BETr'+model] = cube[1,:,:]
    Fut20['BETe'+model] = cube[2,:,:]
    Fut20['BLD'+model] = cube[0,:,:]
    Fut20['NLE'+model] = cube[4,:,:]
    Fut20['NLD'+model] = cube[3,:,:]
    Fut20['C3'+model] = cube[5,:,:]
    Fut20['C4'+model] = cube[8,:,:]
    Fut20['SBE'+model] = cube[12,:,:]
    Fut20['SBD'+model] = cube[11,:,:]
    Fut20['BS'+model] = cube[15,:,:]

Fut20_BETr = (Fut20['BETrGFDL']+Fut20['BETrIPSL']+Fut20['BETrMIROC']+Fut20['BETrHadGEM2'])/4
Fut20_BETe = (Fut20['BETeGFDL']+Fut20['BETeIPSL']+Fut20['BETeMIROC']+Fut20['BETeHadGEM2'])/4
Fut20_BLD = (Fut20['BLDGFDL']+Fut20['BLDIPSL']+Fut20['BLDMIROC']+Fut20['BLDHadGEM2'])/4
Fut20_NLE = (Fut20['NLEGFDL']+Fut20['NLEIPSL']+Fut20['NLEMIROC']+Fut20['NLEHadGEM2'])/4
Fut20_NLD = (Fut20['NLDGFDL']+Fut20['NLDIPSL']+Fut20['NLDMIROC']+Fut20['NLDHadGEM2'])/4
Fut20_C3 = (Fut20['C3GFDL']+Fut20['C3IPSL']+Fut20['C3MIROC']+Fut20['C3HadGEM2'])/4
Fut20_C4 = (Fut20['C4GFDL']+Fut20['C4IPSL']+Fut20['C4MIROC']+Fut20['C4HadGEM2'])/4
Fut20_SBE = (Fut20['SBEGFDL']+Fut20['SBEIPSL']+Fut20['SBEMIROC']+Fut20['SBEHadGEM2'])/4
Fut20_SBD = (Fut20['SBDGFDL']+Fut20['SBDIPSL']+Fut20['SBDMIROC']+Fut20['SBDHadGEM2'])/4
Fut20_BS = (Fut20['BSGFDL']+Fut20['BSIPSL']+Fut20['BSMIROC']+Fut20['BSHadGEM2'])/4

vmin=0
vmax=0.5


'''
a = np.array([[0,1]])
b = np.array([[0,2]])

fig, axs = plt.subplots(10, 4, constrained_layout=True)

ax = axs[0,0]
mesh1=ax.pcolormesh(a)
ax = axs[0,1]
mesh2=ax.pcolormesh(b)

plt.colorbar(mesh1, ax=axs[:6, :], location='right')
plt.colorbar(mesh2, ax=axs[6:, :], location='right')
plt.show()
'''
import cartopy.feature as cfeature
countries = cfeature.NaturalEarthFeature(category='cultural',
                                         name='admin_0_countries',
                                         scale='50m',
                                         facecolor='none')



ax = plt.subplot(10,4,1)
mesh1=iplt.pcolormesh(CCI_BETr, vmin=vmin, vmax=vmax)
plt.text(-0.3, 0.4, 'BETr', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 
plt.title('a) PD CCI')

plt.subplot(10,4,2)
iplt.pcolormesh(PD_BETr, vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 
plt.title('b) PD MMM')

plt.subplot(10,4,3)
iplt.pcolormesh(Fut15_BETr-PD_BETr, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 
plt.title('c) 1.5C - PD')

plt.subplot(10,4,4)
iplt.pcolormesh(Fut20_BETr-PD_BETr, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 
plt.title('d) 2.0C - PD')


ax=plt.subplot(10,4,5)
iplt.pcolormesh(CCI_NLE, vmin=vmin, vmax=vmax)
plt.text(-0.3, 0.4, 'NLE', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,6)
iplt.pcolormesh(PD_NLE, vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,7)
iplt.pcolormesh(Fut15_NLE-PD_NLE, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,8)
iplt.pcolormesh(Fut20_NLE-PD_NLE, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 


ax=plt.subplot(10,4,9)
iplt.pcolormesh(CCI_NLD, vmin=vmin, vmax=vmax)
plt.text(-0.3, 0.4, 'NLD', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,10)
iplt.pcolormesh(PD_NLD, vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,11)
iplt.pcolormesh(Fut15_NLD-PD_NLD, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,12)
iplt.pcolormesh(Fut20_NLD-PD_NLD, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 



ax=plt.subplot(10,4,13)
iplt.pcolormesh(CCI_C3, vmin=vmin, vmax=vmax)
plt.text(-0.3, 0.4, 'C3G', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,14)
iplt.pcolormesh(PD_C3, vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,15)
iplt.pcolormesh(Fut15_C3-PD_C3, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,16)
iplt.pcolormesh(Fut20_C3-PD_C3, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 



ax=plt.subplot(10,4,17)
iplt.pcolormesh(CCI_C4, vmin=vmin, vmax=vmax)
plt.text(-0.3, 0.4, 'C4G', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,18)
iplt.pcolormesh(PD_C4, vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,19)
iplt.pcolormesh(Fut15_C4-PD_C4, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,20)
iplt.pcolormesh(Fut20_C4-PD_C4, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 



ax=plt.subplot(10,4,21)
iplt.pcolormesh(CCI_BS, vmin=vmin, vmax=vmax)
plt.text(-0.3, 0.4, 'BS', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,22)
iplt.pcolormesh(PD_BS, vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,23)
iplt.pcolormesh(Fut15_BS-PD_BS, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,24)
mesh3=iplt.pcolormesh(Fut20_BS-PD_BS, cmap='bwr', vmin=-0.1, vmax=0.1)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 



vmin=0
vmax2=0.1

ax=plt.subplot(10,4,25)
iplt.pcolormesh(CCI_BETe, vmin=vmin, vmax=vmax2)
plt.text(-0.3, 0.4, 'BETe', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,26)
iplt.pcolormesh(PD_BETe, vmin=vmin, vmax=vmax2)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,27)
iplt.pcolormesh(Fut15_BETe-PD_BETe, cmap='bwr', vmin=-0.01, vmax=0.01)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,28)
iplt.pcolormesh(Fut20_BETe-PD_BETe, cmap='bwr', vmin=-0.01, vmax=0.01)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 


ax=plt.subplot(10,4,29)
iplt.pcolormesh(CCI_BLD, vmin=vmin, vmax=vmax2)
plt.text(-0.3, 0.4, 'BLD', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,30)
iplt.pcolormesh(PD_BLD, vmin=vmin, vmax=vmax2)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,31)
iplt.pcolormesh(Fut15_BLD-PD_BLD, cmap='bwr', vmin=-0.01, vmax=0.01)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,32)
iplt.pcolormesh(Fut20_BLD-PD_BLD, cmap='bwr', vmin=-0.01, vmax=0.01)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 


ax=plt.subplot(10,4,33)
iplt.pcolormesh(CCI_SBE, vmin=vmin, vmax=vmax2)
plt.text(-0.3, 0.4, 'SbE', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,34)
iplt.pcolormesh(PD_SBE, vmin=vmin, vmax=vmax2)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,35)
iplt.pcolormesh(Fut15_SBE-PD_SBE, cmap='bwr', vmin=-0.01, vmax=0.01)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,36)
iplt.pcolormesh(Fut20_SBE-PD_SBE, cmap='bwr', vmin=-0.01, vmax=0.01)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 


ax=plt.subplot(10,4,37)
iplt.pcolormesh(CCI_SBD, vmin=vmin, vmax=vmax2)
plt.text(-0.3, 0.4, 'SbD', transform=ax.transAxes)  
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,38)
mesh2=iplt.pcolormesh(PD_SBD, vmin=vmin, vmax=vmax2)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,39)
iplt.pcolormesh(Fut15_SBD-PD_SBD, cmap='bwr', vmin=-0.01, vmax=0.01)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

plt.subplot(10,4,40)
mesh4=iplt.pcolormesh(Fut20_SBD-PD_SBD, cmap='bwr', vmin=-0.01, vmax=0.01)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 

#ax.add_feature(cfeature.BORDERS, linewidth=0.5)

'''
      #Position: Right,Bottom, width, height
cax = plt.axes([0.85, 0.45, 0.02, 0.4])
plt.colorbar(cax=cax)

      #Position: Right,Bottom, width, height
cax = plt.axes([0.85, 0.145, 0.02, 0.25])
plt.colorbar(cax=cax)

'''
      #Position: Right,Bottom, width, height
cax = plt.axes([0.005, 0.45, 0.02, 0.4])
plt.colorbar(mesh1,cax=cax, ticks=[0,0.1,0.2,0.3,0.4,0.5])

      #Position: Right,Bottom, width, height
cax = plt.axes([0.005, 0.145, 0.02, 0.25])
plt.colorbar(mesh2,cax=cax, ticks=[0,0.1])

      #Position: Right,Bottom, width, height
cax = plt.axes([0.85, 0.45, 0.02, 0.4])
plt.colorbar(mesh3,cax=cax, ticks=[-0.1,-0.05,0,0.05,0.1])

      #Position: Right,Bottom, width, height
cax = plt.axes([0.85, 0.145, 0.02, 0.25])
plt.colorbar(mesh3,cax=cax, ticks=[-0.1,-0.05,0,0.05,0.1])

plt.show()



'''
plt.subplot(10,4,3)
plt.subplot(10,4,4)
plt.subplot(10,4,5)
plt.subplot(10,4,6)
plt.subplot(10,4,7)
plt.subplot(10,4,8)
plt.subplot(10,4,9)
plt.subplot(10,4,10)
plt.subplot(10,4,11)
plt.subplot(10,4,12)
plt.subplot(10,4,13)
plt.subplot(10,4,14)
plt.subplot(10,4,15)
plt.subplot(10,4,16)
plt.subplot(10,4,17)
plt.subplot(10,4,18)
plt.subplot(10,4,19)
plt.subplot(10,4,20)
plt.subplot(10,4,21)
plt.subplot(10,4,22)
plt.subplot(10,4,23)
plt.subplot(10,4,24)
plt.subplot(10,4,25)
plt.subplot(10,4,26)
plt.subplot(10,4,27)
plt.subplot(10,4,28)
plt.subplot(10,4,29)
plt.subplot(10,4,30)
plt.subplot(10,4,31)
plt.subplot(10,4,32)
plt.subplot(10,4,33)
plt.subplot(10,4,34)
plt.subplot(10,4,35)
plt.subplot(10,4,36)
plt.subplot(10,4,37)
plt.subplot(10,4,38)
plt.subplot(10,4,39)
plt.subplot(10,4,40)



                         #Position: Left, Bottom, width, height
colorbar_axes = plt.gcf().add_axes([0.3, 0.45, 0.45, 0.01])
colorbar = plt.colorbar(mesh1, colorbar_axes, orientation='vertical') 

                         #Position: Left, Bottom, width, height
colorbar_axes = plt.gcf().add_axes([0.3, 0.06, 0.45, 0.01])
colorbar = plt.colorbar(mesh2, colorbar_axes, orientation='vertical') 

'''
