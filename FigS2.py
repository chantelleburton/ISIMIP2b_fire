
import os
import glob
import sys
sys.path.append('~kwilliam/other_fcm/jules_py/trunk/jules/')
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
from collections import OrderedDict 
import pandas as pd
from matplotlib.ticker import FormatStrFormatter


def Regrid(template, data):
    cs_new = iris.coord_systems.GeogCS(6371229.)
    template.coord('latitude').coord_system = cs_new
    template.coord('longitude').coord_system = cs_new
    data.coord('latitude').coord_system = cs_new
    data.coord('longitude').coord_system = cs_new
    data = data.regrid(template, iris.analysis.Linear()) 
    return data

def CalcRegion(mydata, GFEDregions, key):
    region = mydata.copy()

    region.data = np.ma.masked_where(GFEDregions.data != key, mydata.data)

    coords = ('longitude', 'latitude')
    for coord in coords:
        if not region.coord(coord).has_bounds():
            region.coord(coord).guess_bounds()
    grid_weights = iris.analysis.cartography.area_weights(region)

    TOTAL = region.collapsed(coords, iris.analysis.SUM, weights=grid_weights)/1E12 
    result = TOTAL.data
    print (result)
    return result


def main(): 

    Region_dict = OrderedDict([
        ('a) BONA',1),
        ('b) TENA',2),
        ('c) CEAM',3),
        ('d) NHSA',4),
        ('e) SHSA',5),
        ('f) EURO',6),
        ('g) MIDE',7),
        ('h) NHAF',8),
        ('i) SHAF',9),
        ('j) BOAS',10),
        ('k) CEAS',11),
        ('l) SEAS',12),
        ('m) EQAS',13),
        ('n) AUST',14)])

    mydata = jules.load_cube('/MyScratchFolder/IPSL/ipsl-cm5a-lr_c20c.gen_ann_gb.2001.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb')).collapsed(['time'], iris.analysis.MEAN)
  
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions')
    GFEDregions = iris.load_cube('$DATADIR/GFED/GFED4.1s/GFED4.1s_1997-2016.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN)
    GFEDregions = Regrid(template=mydata, data=GFEDregions)

    d = {}
    
    models = ['HadGEM2','MIROC','GFDL','IPSL']
    RCPS = ['2010-2019','RCP60_1p5','RCP60_2p0']
    for model in models:
        for RCP in RCPS:
            for key in Region_dict:

                folder = '/MyScratchFolder/'
                BA = jules.load_cube(folder+model+'/'+RCP+'/Years.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb')).collapsed(['time'], iris.analysis.MEAN)*86400*360
                BA.data[np.where(np.isnan(BA.data))] = 0.0

                var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'LSM')
                new_landmask = landmask = iris.load_cube('/MyPaperFolder/landfrac_isimip_0p5.nc', var_constraint)
                landmask = landmask.regrid(BA, iris.analysis.Linear())
                BA = BA*landmask

                #Constrain by region
                d[model+RCP+key] = CalcRegion(BA, GFEDregions, Region_dict[key])

    print('done models')


### Separate bars for each model 
    n = 1
    for key in Region_dict:
        index = np.arange(4)
        bar_width = 0.2
        fontsize=10

        PD = [ d["HadGEM22010-2019"+key],d["MIROC2010-2019"+key],d["GFDL2010-2019"+key],d["IPSL2010-2019"+key] ]
        Fut1pt5 = [ d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key] ]
        Fut2pt0 = [ d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key] ]

        plt.subplot(7,2,n) 
        ax = plt.gca()
        plt.bar(index,PD, bar_width, color='blue')
        plt.bar(index + bar_width, Fut1pt5,bar_width,  color='orange')
        plt.bar(index + bar_width*2, Fut2pt0, bar_width, color='grey')
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        plt.yticks(fontsize=fontsize)
        if n == 7:
            plt.ylabel('Burnt Area (Mkm$^2$)', fontsize=fontsize)
            plt.xticks(index + bar_width, ('','','',''))
        elif n == 13:
            plt.xticks(index + bar_width, ('HG2','MIROC','GFDL','IPSL'),fontsize=fontsize)
        elif n == 14:
            plt.xticks(index + bar_width, ('HG2','MIROC','GFDL','IPSL'),fontsize=fontsize)
        else:
            plt.xticks(index + bar_width, ('','','',''))
   
        plt.title(str(key), fontsize=fontsize)
        n = n+1

    plt.show()



if __name__ == '__main__':
    main()






