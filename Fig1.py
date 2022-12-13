
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
    f = open('/scratch/cburton/scratch/ISIMIP_PAPER/OctUpdates/Fig1.dat','a')
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

    mydata = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/IPSL/ipsl-cm5a-lr_c20c.gen_ann_gb.2001.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb')).collapsed(['time'], iris.analysis.MEAN)
  
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions')
    GFEDregions = iris.load_cube('/data/cr1/cburton/GFED/GFED4.1s/GFED4.1s_1997-2016.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN)
    GFEDregions = Regrid(template=mydata, data=GFEDregions)

    d = {}
    '''
    GFED = iris.load_cube('/data/cr1/cburton/GFED/GFED4s_AnnualTotalBA_2001-2016.nc').collapsed(['time'], iris.analysis.MEAN)
    GFED = Regrid(template=mydata, data=GFED)

    for key in Region_dict:
        d["GFED"+key] = CalcRegion(GFED, GFEDregions, Region_dict[key])
    print('done GFED')

    CCI = iris.load_cube('/data/cr1/cburton/FireCCI/2001-2016_Annual_BAFrac_CB.nc').collapsed(['time'], iris.analysis.MEAN)
    CCI = Regrid(template=mydata, data=CCI)
    for key in Region_dict:
        d["CCI"+key] = CalcRegion(CCI, GFEDregions, Region_dict[key])
    print('done CCI')
    '''
    models = ['HadGEM2','MIROC','GFDL','IPSL']
    RCPS = ['2010-2019','RCP60_1p5','RCP60_2p0']
    for model in models:
        for RCP in RCPS:
            for key in Region_dict:

                folder = '/scratch/cburton/scratch/ISIMIP_PAPER/'
                BA = jules.load_cube(folder+model+'/'+RCP+'/Years.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb')).collapsed(['time'], iris.analysis.MEAN)*86400*360
                BA.data[np.where(np.isnan(BA.data))] = 0.0

                var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'LSM')
                new_landmask = landmask = iris.load_cube('/data/users/hadea/jules_ancils/isimip2/landfrac_isimip_0p5.nc', var_constraint)
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
        f = open('/scratch/cburton/scratch/ISIMIP_PAPER/OctUpdates/Fig1.dat','a')
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
        plt.savefig('/scratch/cburton/scratch/ISIMIP_PAPER/OctUpdates/'+str(key)+'_Range.png')
        



if __name__ == '__main__':
    main()







'''





### Separate bars for each model 
    for key in Region_dict:
        index = np.arange(4)
        bar_width = 0.2
        fontsize=27

        PD = [ d["HadGEM22001-2016"+key],d["MIROC2001-2016"+key],d["GFDL2001-2016"+key],d["IPSL2001-2016"+key] ]
        Fut1pt5 = [ d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key] ]
        Fut2pt0 = [ d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key] ]

        plt.figure(figsize=(8, 6))
        ax = plt.gca()
        plt.bar(index,PD, bar_width, color='blue')
        plt.bar(index + bar_width, Fut1pt5,bar_width,  color='orange')
        plt.bar(index + bar_width*2, Fut2pt0, bar_width, color='grey')
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        plt.yticks(fontsize=fontsize)
        plt.xticks(index + bar_width, ('HadGEM2','MIROC','GFDL','IPSL'),fontsize=fontsize)
        plt.title(str(key), fontsize=fontsize)
        plt.savefig('/scratch/cburton/scratch/ISIMIP_PAPER/SeptUpdates/'+str(key)+'_Bars.png')




### Error bars for mean, min, max, WITH observations

    for key in Region_dict:
    # setup the dataframe
        Exps = ['Obs', 'PD', '1.5', '2.0']
        Mean = [ np.mean([d["GFED"+key], d["CCI"+key]]),
            np.mean([d["HadGEM22001-2016"+key],d["MIROC2001-2016"+key],d["GFDL2001-2016"+key],d["IPSL2001-2016"+key]]),
            np.mean([d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key]]),
            np.mean([d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key]]) ]

        _min = [ np.min([d["GFED"+key], d["CCI"+key]]),
            np.min([d["HadGEM22001-2016"+key],d["MIROC2001-2016"+key],d["GFDL2001-2016"+key],d["IPSL2001-2016"+key]]),
            np.min([d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key]]),
            np.min([d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key]]) ]

        _max = [ np.max([d["GFED"+key], d["CCI"+key]]),
            np.max([d["HadGEM22001-2016"+key],d["MIROC2001-2016"+key],d["GFDL2001-2016"+key],d["IPSL2001-2016"+key]]),
            np.max([d["HadGEM2RCP60_1p5"+key],d["MIROCRCP60_1p5"+key],d["GFDLRCP60_1p5"+key],d["IPSLRCP60_1p5"+key]]),
            np.max([d["HadGEM2RCP60_2p0"+key],d["MIROCRCP60_2p0"+key],d["GFDLRCP60_2p0"+key],d["IPSLRCP60_2p0"+key]]) ]


        df = pd.DataFrame({'Exps':Exps,'Mean':Mean, 'Min': _min, 'Max': _max})
        df = (df.assign(yerr_min = df.Mean-df.Min)
            .assign(yerr_max=df.Max-df.Mean))
        plt.figure(figsize=(8, 6))
        plt.bar(x='Exps', height='Mean', yerr=df[['yerr_min', 'yerr_max']].T.values, capsize=10, data=df, color=['grey', 'blue', 'orange', 'lightgrey'])
        plt.title(str(key))
        #plt.show()
        plt.savefig('/scratch/cburton/scratch/ISIMIP_PAPER/SeptUpdates/'+str(key)+'.png')


'''


'''
#Exps = ('HadGEM2','MIROC','GFDL','IPSL')
index = np.arange(4)
PD = [0.5,0.4,0.5,0.6]
pot15 = [ 1.0,1.2,1.1,0.9 ]
pot20 = [ 1.5,1.6,1.7,1.8 ]
bar_width = 0.2
#df = pd.DataFrame({'Exps':Exps,'PD':PD, '1.5': pot15, '2.0': pot20})
plt.figure(figsize=(8, 6))
plt.bar(index,PD, bar_width, color='blue')
plt.bar(index + bar_width, pot15,bar_width,  color='orange')
plt.bar(index + bar_width*2, pot20, bar_width, color='grey')
plt.yticks(fontsize=18)
plt.xticks(index + bar_width, ('HadGEM2','MIROC','GFDL','IPSL'),fontsize=18)
plt.show()
exit()


'''
