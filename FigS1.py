
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
import cftime


def Regrid(cube):
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb')
    file1 = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/GFDL/gfdl-esm2m_c20c.gen_ann_gb.1861.nc', var_constraint,  missingdata=np.ma.masked)*86400*360
    cs_new = iris.coord_systems.GeogCS(6371229.)
    file1.coord('latitude').coord_system = cs_new
    file1.coord('longitude').coord_system = cs_new
    cube.coord('latitude').coord_system = cs_new
    cube.coord('longitude').coord_system = cs_new
    cube = cube.regrid(file1, iris.analysis.Linear())
    return cube


### Get GFED data
GFED = iris.load_cube('/data/cr1/cburton/GFED/GFED4s_AnnualTotalBA_2001-2016.nc')
GFED = Regrid(GFED)
GFED = GFED.collapsed(['time'], iris.analysis.MEAN)

### Get Fire CCI data
CCI = iris.load_cube('/data/cr1/cburton/FireCCI/2001-2016_Annual_BAFrac_CB.nc')
CCI = Regrid(CCI)
CCI = CCI.collapsed(['time'], iris.analysis.MEAN)




####### Get JULES data ######
folder = '/scratch/cburton/scratch/ISIMIP_PAPER'
SectoAnn = 86400*360

var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb')

cube = jules.load_cube(folder+'/HadGEM2/2001-2016/Years.nc', var_constraint, missingdata=np.ma.masked)
HG2_PD = (cube.collapsed(['time'], iris.analysis.MEAN)*SectoAnn)

cube = jules.load_cube(folder+'/MIROC/2001-2016/Years.nc', var_constraint, missingdata=np.ma.masked)
MIROC_PD = (cube.collapsed(['time'], iris.analysis.MEAN)*SectoAnn)

cube = jules.load_cube(folder+'/IPSL/2001-2016/Years.nc', var_constraint, missingdata=np.ma.masked)
IPSL_PD = (cube.collapsed(['time'], iris.analysis.MEAN)*SectoAnn)

cube = jules.load_cube(folder+'/GFDL/2001-2016/Years.nc', var_constraint, missingdata=np.ma.masked)
GFDL_PD = (cube.collapsed(['time'], iris.analysis.MEAN)*SectoAnn)
PD = (HG2_PD+MIROC_PD+IPSL_PD+GFDL_PD)/4


cube = jules.load_cube(folder+'/HadGEM2/RCP26_1p5/Years.nc', var_constraint, missingdata=np.ma.masked)
HG2_RCP26 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

cube = jules.load_cube(folder+'/MIROC/RCP26_1p5/Years.nc', var_constraint, missingdata=np.ma.masked)
MIROC_RCP26 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

cube = jules.load_cube(folder+'/IPSL/RCP26_1p5/Years.nc', var_constraint, missingdata=np.ma.masked)
IPSL_RCP26 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

RCP26_15 = ((HG2_RCP26+MIROC_RCP26+IPSL_RCP26)/3)


cube = jules.load_cube(folder+'/HadGEM2/RCP60_1p5/Years.nc', var_constraint, missingdata=np.ma.masked)
HG2_RCP60 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

cube = jules.load_cube(folder+'/MIROC/RCP60_1p5/Years.nc', var_constraint, missingdata=np.ma.masked)
MIROC_RCP60 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

cube = jules.load_cube(folder+'/IPSL/RCP60_1p5/Years.nc', var_constraint, missingdata=np.ma.masked)
IPSL_RCP60 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

cube = jules.load_cube(folder+'/GFDL/RCP60_1p5/Years.nc', var_constraint, missingdata=np.ma.masked)
GFDL_RCP60 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

RCP60_15 = ((HG2_RCP60+MIROC_RCP60+IPSL_RCP60+GFDL_RCP60)/4) 
RCP60_15change = RCP60_15 - PD


cube = jules.load_cube(folder+'/HadGEM2/RCP60_2p0/Years.nc', var_constraint, missingdata=np.ma.masked)
HG2_RCP60_20 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

cube = jules.load_cube(folder+'/MIROC/RCP60_2p0/Years.nc', var_constraint, missingdata=np.ma.masked)
MIROC_RCP60_20 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

cube = jules.load_cube(folder+'/IPSL/RCP60_2p0/Years.nc', var_constraint, missingdata=np.ma.masked)
IPSL_RCP60_20 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

cube = jules.load_cube(folder+'/GFDL/RCP60_2p0/Years.nc', var_constraint, missingdata=np.ma.masked)
GFDL_RCP60_20 = (cube.collapsed(['time'], iris.analysis.MEAN))*SectoAnn 

RCP60_20 = ((HG2_RCP60_20+MIROC_RCP60_20+IPSL_RCP60_20+GFDL_RCP60_20)/4) 
RCP60_20change = RCP60_15 - PD


RCP60_2015change = RCP60_20 - RCP60_15


####### Make Plots ######
plt.subplot(3,3,1)
iplt.pcolormesh(GFED, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('a) Present Day: GFED4s')
plt.colorbar(orientation='horizontal')

plt.subplot(3,3,2)
iplt.pcolormesh(CCI, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('b) Present Day: Fire CCI')
plt.colorbar(orientation='horizontal')

plt.subplot(3,3,3)
iplt.pcolormesh(PD, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('c) Present Day: MMM')
plt.colorbar(orientation='horizontal')


plt.subplot(3,3,4)
iplt.pcolormesh(RCP26_15, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('d) RCP2.6 1.5$^o$C')
plt.colorbar(orientation='horizontal')

plt.subplot(3,3,5)
iplt.pcolormesh(RCP60_15, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('e) RCP6.0 1.5$^o$C')
plt.colorbar(orientation='horizontal')

plt.subplot(3,3,6)
iplt.pcolormesh(RCP60_15change, vmin=-0.04, vmax=0.04, cmap='bwr')
plt.gca().coastlines()
plt.title('f) RCP6.0 1.5$^o$C - PD')
plt.colorbar(orientation='horizontal', ticks=[-0.04,-0.02,0.0,0.02,0.04])


plt.subplot(3,3,7)
iplt.pcolormesh(RCP60_20, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('g) RCP6.0 2.0$^o$C')
plt.colorbar(label='frac', orientation='horizontal')

plt.subplot(3,3,8)
iplt.pcolormesh(RCP60_2015change, vmin=-0.02, vmax=0.02, cmap='bwr')
plt.gca().coastlines()
plt.title('h) RCP6.0 2.0$^o$C - 1.5$^o$C')
plt.colorbar(label='frac', orientation='horizontal')

plt.subplot(3,3,9)
iplt.pcolormesh(RCP60_20change, vmin=-0.04, vmax=0.04, cmap='bwr')
plt.gca().coastlines()
plt.title('i) RCP6.0 2.0$^o$C - PD')
plt.colorbar(label='frac', orientation='horizontal', ticks=[-0.04,-0.02,0.0,0.02,0.04])



plt.show()


