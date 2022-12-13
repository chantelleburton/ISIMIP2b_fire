
import os
import glob
import sys
sys.path.append('/net/home/h03/kwilliam/other_fcm/jules_py/trunk/jules/')
#~kwilliam/other_fcm/jules_py/trunk/jules/jules.py
import jules
import iris
import cf_units
import numpy as np
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from scipy import stats
import iris.analysis.stats
import scipy.stats
import iris.plot as iplt
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib


###Six Maps
countries = cfeature.NaturalEarthFeature(category='cultural',
                                         name='admin_0_countries',
                                         scale='50m',
                                         facecolor='none')


folder = '/scratch/cburton/scratch/ISIMIP_PAPER/HadGEM2'

###### Plot 6 maps ####

TS2 = jules.load_cube(folder+'/RCP60_2090s.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'tstar_gb'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
TS1 = jules.load_cube(folder+'/2010-2019/Years.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'tstar_gb'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
TSchange = TS2-TS1


PR2= jules.load_cube(folder+'/RCP60_2090s.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'precip'),missingdata=np.ma.masked)*86400*360
PR2 = PR2.collapsed(['time'], iris.analysis.MEAN)
PR1 = jules.load_cube(folder+'/2010-2019/Years.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'precip'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)*86400*360
PRchange = PR2-PR1

cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/HadGEM2/hadgem2-es_c20c.gen_ann_gb.1921.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'tstar_gb'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)


QW2= iris.load('/scratch/cburton/scratch/ISIMIP_PAPER/DRIVE/hurs_day_HadGEM2-ES_rcp60_r1i1p1_EWEMBI_20910101-20991231.nc')
QW2=QW2[1].collapsed(['time'], iris.analysis.MEAN)
QW1 = iris.load('/scratch/cburton/scratch/ISIMIP_PAPER/DRIVE/hurs_day_HadGEM2-ES_rcp60_r1i1p1_EWEMBI_20110101-20201231.nc')
QW1=QW1[1].collapsed(['time'], iris.analysis.MEAN)

QWchange = QW2-QW1
QWchange = QWchange.regrid(cube, iris.analysis.Linear())
RHchange = cube.copy()
RHchange.data = QWchange.data


#fsmc_gb = available soil moisture / smc_tot = total column SM / sthu = unfrozen soil moisture frac 
SM2= jules.load_cube(folder+'/RCP60_2090s_PFTlayer.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'sthu'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
SM2 = SM2[0,:,:]
SM1 = jules.load_cube(folder+'/2010-2019/PFTYears.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'sthu'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
SM1 = SM1[0,:,:]
SMchange = SM2-SM1


CV2= jules.load_cube(folder+'/RCP60_2090s.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'cv'),missingdata=np.ma.masked)
CV2 = CV2.collapsed(['time'], iris.analysis.MEAN)
CV1 = jules.load_cube(folder+'/2010-2019/Years.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'cv'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
CVchange = CV2-CV1


Ag2= jules.load_cube(folder+'/RCP60_2090s_PFTlayer.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'frac'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
Ag2 = Ag2[6,:,:]+Ag2[9,:,:]

Ag1 = jules.load_cube(folder+'/2010-2019/PFTYears.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'frac'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
Ag1 = Ag1[6,:,:]+Ag1[9,:,:]

Agchange = Ag2-Ag1
iris.save(Agchange, '/scratch/cburton/scratch/ISIMIP_PAPER/Agchange.nc')


#####~~~~~~~~PLOTS~~~~~~~#######

plt.subplot(2,3,1)
iplt.pcolormesh(TSchange, vmin=0, vmax=6, cmap='Reds')
plt.gca().coastlines()
plt.colorbar(orientation='horizontal', label='Celsius',ticks=[0,2,4,6])
plt.title('a) Temperature')

plt.subplot(2,3,2)
iplt.pcolormesh(PRchange, vmin=-200, vmax=200, cmap='coolwarm_r')
plt.gca().coastlines()
plt.colorbar(orientation='horizontal', label='mm/yr',ticks=[-200,-100,0,100,200])
plt.title('b) Precipitation')

plt.subplot(2,3,3)
iplt.pcolormesh(RHchange, vmin=-10, vmax=10, cmap='coolwarm_r')
plt.gca().coastlines()
plt.colorbar(orientation='horizontal', label='kg/kg', ticks=[-10,-5,0,5,10])
plt.title('c) Relative Humidity')

plt.subplot(2,3,4)
iplt.pcolormesh(SMchange, vmin=-0.2, vmax=0.2, cmap='coolwarm_r')
plt.gca().coastlines()
plt.colorbar(orientation='horizontal', label='Frac of saturation', ticks=[-0.2,-0.1,0,0.1,0.2])
plt.title('d) Top layer soil moisture')

plt.subplot(2,3,5)
iplt.pcolormesh(CVchange, vmin=-10, vmax=10, cmap='BrBG')
plt.gca().coastlines()
plt.colorbar(orientation='horizontal',label='kgC/m$^2$', ticks=[-10,-5,0,5,10])
plt.title('e) Veg Carbon')

plt.subplot(2,3,6)
iplt.pcolormesh(Agchange, vmin=-0.2, vmax=0.2, cmap='bwr')
plt.gca().coastlines()
plt.colorbar(orientation='horizontal', label='Gridbox fraction',ticks=[-0.2,-0.1,0,0.1,0.2])
plt.title('f) Agriculture')


plt.show()
raise Exception('Stop')


'''
plt.subplot(2,4,3)
iplt.pcolormesh(ETchange, vmin=-0.8, vmax=0.8,cmap='coolwarm')
plt.gca().coastlines()
plt.gca().add_feature(countries, lw=0.3) 
plt.colorbar(orientation='horizontal', label='kg/m$^2$/day',ticks=[-0.8,-0.4,0,0.4,0.8])
plt.title('ET')

'''


