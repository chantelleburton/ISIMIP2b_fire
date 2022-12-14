
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
    #Save text out to a file
    f = open('/MyScratchFolder/OctUpdates/Fig1.dat','a')
    np.savetxt(f,('Region','MeanPD','Mean15','Mean20'),newline=' ',fmt='  %s')
    f.write('\n')
    f.close()

    Region_dict = OrderedDict([
        ('BONA',1),
        ('TENA',2),
        ('CEAM',3),
        ('NHSA',4),
        ('SHSA',5),
        ('EURO',6),
        ('MIDE',7),
        ('NHAF',8),
        ('SHAF',9),
        ('BOAS',10),
        ('CEAS',11),
        ('SEAS',12),
        ('EQAS',13),
        ('AUST',14)])

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


#### Error bars for mean, min, max, WITHOUT observations 

    for key in Region_dict:
    # setup the dataframe
        Exps = ['PD', '1.5', '2.0']
        Mean = [ 
            np.mean([d["HadGEM22010-2019"+key],d["MIROC2010-2019"+key],d["GFDL2010-2019"+key],d["IPSL2010-2019"+key]]),
            np.mean([d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key]]),
            np.mean([d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key]]) ]

        _min = [ 
            np.min([d["HadGEM22010-2019"+key],d["MIROC2010-2019"+key],d["GFDL2010-2019"+key],d["IPSL2010-2019"+key]]),
            np.min([d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key]]),
            np.min([d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key]]) ]

        _max = [ 
            np.max([d["HadGEM22010-2019"+key],d["MIROC2010-2019"+key],d["GFDL2010-2019"+key],d["IPSL2010-2019"+key]]),
            np.max([d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key]]),
            np.max([d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key]]) ]


        MeanPD = np.mean([d["HadGEM22010-2019"+key],d["MIROC2010-2019"+key],d["GFDL2010-2019"+key],d["IPSL2010-2019"+key]])
        Mean15 = np.mean([d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key]])
        Mean20 = np.mean([d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key]])
        print (key, MeanPD, Mean15, Mean20)

        #Save text out to a file
        f = open('/MyScratchFolder/OctUpdates/Fig1.dat','a')
        np.savetxt(f,(key,MeanPD,Mean15,Mean20),newline=' ',fmt='  %s')
        f.write('\n')
        f.close()

        
        fontsize=31
        df = pd.DataFrame({'Exps':Exps,'Mean':Mean, 'Min': _min, 'Max': _max})
        df = (df.assign(yerr_min = df.Mean-df.Min)
            .assign(yerr_max=df.Max-df.Mean))
        plt.figure(figsize=(8, 6))
        ax = plt.gca()
        plt.bar(x='Exps', height='Mean', width=0.5, yerr=df[['yerr_min', 'yerr_max']].T.values, capsize=10, data=df, color=['blue', 'orange', 'lightgrey'])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        plt.yticks(fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        #plt.title(str(key), fontsize=fontsize)
        plt.savefig('/MyScratchFolder/OctUpdates/'+str(key)+'_Range.png')
        



if __name__ == '__main__':
    main()




