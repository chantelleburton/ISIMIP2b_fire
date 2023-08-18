
import os
import glob
import sys
sys.path.append('/net/home/h03/kwilliam/other_fcm/jules_py/trunk/jules/')
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
### Global totals
### Model Data
var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'cv')
cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/GFDL/2010-2019/GenAnn.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN) 
cube.data[np.where(np.isnan(cube.data))] = 0.0
coords = ('longitude', 'latitude')
for coord in coords:
    if not cube.coord(coord).has_bounds():
        cube.coord(coord).guess_bounds()
grid_weights = iris.analysis.cartography.area_weights(cube)
cube = cube.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights = grid_weights) /1E12
print(cube.data)

# Get cveg data
AGB = iris.load_cube('/data/users/dkelley/jules_benchmarking/outputs/CCI0.5/ESACCI-BIOMASS-L4-AGB-MERGED-0.5Degree-2010-fv3.0-ensembles-3.nc', 'agb')/10 # convert Mg/ha to kg/m2
print(AGB)
for coord in coords:
    if not AGB.coord(coord).has_bounds():
        AGB.coord(coord).guess_bounds()
grid_weights = iris.analysis.cartography.area_weights(AGB)
AGB = AGB.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights = grid_weights) /1E12
print(AGB.data)

exit()
'''
'''
### Obs Data

# Get  CCI data
var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
CCI = iris.load_cube('~/PAPERS/3.FutureFiresGlobal/Data/vegfrac_lc-cci_9PFTs_2010_v1.6.1.nc')
cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/GFDL/2010-2019/PFTYears.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN)
cs_new = iris.coord_systems.GeogCS(6371229.)
CCI.coord('latitude').coord_system = cs_new
CCI.coord('longitude').coord_system = cs_new
cube.coord('latitude').coord_system = cs_new
cube.coord('longitude').coord_system = cs_new
CCI = CCI.regrid(cube, iris.analysis.Linear()) 

CCITrees = CCI[0,:,:]+CCI[1,:,:]+CCI[2,:,:]+CCI[3,:,:]+CCI[4,:,:]
'''


### Model Data
models = ('GFDL', 'IPSL', 'MIROC','HadGEM2')

#Pre-Industrial
PI = {}
for model in models:
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
    cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/'+model+'/1860-1900/Years.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN) 
    PI['DTree'+model] = cube[0,:,:]+cube[3,:,:]
    PI['ETree'+model] = cube[1,:,:]+cube[2,:,:]+cube[4,:,:]
    PI['DSb'+model] = cube[11,:,:]
    PI['ESb'+model] = cube[12,:,:]
    PI['C3'+model] = cube[5,:,:]+cube[6,:,:]+cube[7,:,:]
    PI['C4'+model] = cube[8,:,:]+cube[9,:,:]+cube[10,:,:]


PI_DTree = (PI['DTreeGFDL']+PI['DTreeIPSL']+PI['DTreeMIROC']+PI['DTreeHadGEM2'])/4
PI_ETree = (PI['ETreeGFDL']+PI['ETreeIPSL']+PI['ETreeMIROC']+PI['ETreeHadGEM2'])/4
PI_DSb = (PI['DSbGFDL']+PI['DSbIPSL']+PI['DSbMIROC']+PI['DSbHadGEM2'])/4
PI_ESb = (PI['ESbGFDL']+PI['ESbIPSL']+PI['ESbMIROC']+PI['ESbHadGEM2'])/4
PI_C3 = (PI['C3GFDL']+PI['C3IPSL']+PI['C3MIROC']+PI['C3HadGEM2'])/4
PI_C4 = (PI['C4GFDL']+PI['C4IPSL']+PI['C4MIROC']+PI['C4HadGEM2'])/4


#1.5 deg
Fut15 = {}
for model in models:
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
    cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/'+model+'/RCP60_1p5/PFTlayers.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN) 

    Fut15['DTree'+model] = cube[0,:,:]+cube[3,:,:]
    Fut15['ETree'+model] = cube[1,:,:]+cube[2,:,:]+cube[4,:,:]
    Fut15['DSb'+model] = cube[11,:,:]
    Fut15['ESb'+model] = cube[12,:,:]
    Fut15['C3'+model] = cube[5,:,:]+cube[6,:,:]+cube[7,:,:]
    Fut15['C4'+model] = cube[8,:,:]+cube[9,:,:]+cube[10,:,:]


Fut15_DTree = (Fut15['DTreeGFDL']+Fut15['DTreeIPSL']+Fut15['DTreeMIROC']+Fut15['DTreeHadGEM2'])/4
Fut15_ETree = (Fut15['ETreeGFDL']+Fut15['ETreeIPSL']+Fut15['ETreeMIROC']+Fut15['ETreeHadGEM2'])/4
Fut15_DSb = (Fut15['DSbGFDL']+Fut15['DSbIPSL']+Fut15['DSbMIROC']+Fut15['DSbHadGEM2'])/4
Fut15_ESb = (Fut15['ESbGFDL']+Fut15['ESbIPSL']+Fut15['ESbMIROC']+Fut15['ESbHadGEM2'])/4
Fut15_C3 = (Fut15['C3GFDL']+Fut15['C3IPSL']+Fut15['C3MIROC']+Fut15['C3HadGEM2'])/4
Fut15_C4 = (Fut15['C4GFDL']+Fut15['C4IPSL']+Fut15['C4MIROC']+Fut15['C4HadGEM2'])/4


import cartopy.feature as cfeature
countries = cfeature.NaturalEarthFeature(category='cultural',
                                         name='admin_0_countries',
                                         scale='50m',
                                         facecolor='none')

vmin=-0.4
vmax=0.4
#, vmin=vmin, vmax=vmax
ax = plt.subplot(2,2,1)
iplt.pcolormesh(Fut15_DTree-PI_DTree, cmap='bwr', vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2) 
plt.title('a) Deciduous tree')

plt.subplot(2,2,2)
iplt.pcolormesh(Fut15_ETree-PI_ETree, cmap='bwr', vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2)
plt.title('b) Evergreen tree')

ax=plt.subplot(2,2,3)
iplt.pcolormesh(Fut15_C3-PI_C3, cmap='bwr', vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2)
plt.title('c) C3')

plt.subplot(2,2,4)
mesh = iplt.pcolormesh(Fut15_C4-PI_C4, cmap='bwr', vmin=vmin, vmax=vmax)
plt.gca().coastlines(lw=0.2)
plt.gca().add_feature(countries, lw=0.2)
plt.title('d) C4')

      #Position: Left,Bottom, width, height
cax = plt.axes([0.40, 0.15, 0.25, 0.05])
plt.colorbar(mesh, orientation='horizontal',cax=cax, label='frac')
plt.show()




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


'''
###Test layout

var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'frac')
CCI = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/GFDL/2010-2019/PFTYears.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN)
CCI = CCI[0,:,:]+CCI[1,:,:]+CCI[2,:,:]+CCI[3,:,:]+CCI[4,:,:]

for n in np.arange(8):
    ax = plt.subplot(4,2,n+1)
    mesh1=iplt.pcolormesh(CCI,cmap='Greens', vmin=0, vmax=2)
    plt.gca().coastlines(lw=0.2)
    plt.title('a) PD CCI Tree Frac')
    plt.colorbar()

cax = plt.axes([0.4, 0.7, 0.3, 0.02])
plt.colorbar(mesh1, orientation='horizontal', cax=cax, ticks=[0,0.2, 0.4, 0.6, 0.8, 1.0])

      #Position: Left,Bottom, width, height
cax = plt.axes([0.4, 0.5, 0.3, 0.02])
plt.colorbar(mesh1, orientation='horizontal',cax=cax, ticks=[-0.1,-0.05,0,0.05,0.1])

      #Position: Left,Bottom, width, height
cax = plt.axes([0.4, 0.3, 0.3, 0.02])
plt.colorbar(mesh1, orientation='horizontal',cax=cax, ticks=[0,10,20,30,40])

      #Position: Left,Bottom, width, height
cax = plt.axes([0.4, 0.1, 0.3, 0.02])
plt.colorbar(mesh1, orientation='horizontal',cax=cax, ticks=[-5, -2.5, 0, 2.5, 5])

plt.show()


exit()
'''


