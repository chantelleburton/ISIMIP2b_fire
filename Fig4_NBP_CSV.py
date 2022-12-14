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
from scipy import stats
import iris.analysis.stats
import scipy.stats
import iris.plot as iplt
import numpy as np
from scipy.stats import ks_2samp
import iris.coords as icoords
import datetime
import time
import numpy.ma as ma
from collections import OrderedDict 
import iris.coord_categorisation
import matplotlib
import matplotlib.pyplot as plt
#For using IRIS with SPICE:
#matplotlib.use('Agg')
import pandas as pd
import csv




'''
# First calcualtion for Temperature, saved out
#PI = iris.load_cube('/scratch/cburton/ISIMIP_PAPER/DRIVE/HADGEM2-ES/*historical*.nc')
    files = '/scratch/cburton/ISIMIP_PAPER/DRIVE/'+model+'*/tas_day_'+model+'*_historical_r1i1p1_EWEMBI_18*.nc'
    cubelist  = iris.load(files)
    iris.util.unify_time_units(cubelist)
    PI = cubelist.concatenate_cube()

    PI.convert_units('celsius')
    PI = PI.collapsed(('time'), iris.analysis.MEAN) 

    coords = ('longitude', 'latitude')
    for coord in coords:
        if not PI.coord(coord).has_bounds():
            PI.coord(coord).guess_bounds()
    weights = iris.analysis.cartography.area_weights(PI)
    PIOut = PI.collapsed(coords, iris.analysis.MEAN, weights = weights) 

    Future = iris.load_cube('/scratch/cburton/ISIMIP_PAPER/DRIVE/'+model+'*/Years.nc')
    Future.convert_units('celsius')
    #Change daily data to monthly, within several years
    iris.coord_categorisation.add_month(Future, 'time', name='month')
    iris.coord_categorisation.add_year(Future, 'time', name='year')
    FutureOut = Future.aggregated_by(['year'],iris.analysis.MEAN)

    Temp = FutureOut-PIOut
    iris.save(Temp, '/scratch/cburton/ISIMIP_PAPER/DRIVE/'+model+'_TempCAbovePI.nc')
    print ('saved',model)
'''



'''
### NCO commands to concatenate all NBP cube files for each model (add time dimension)
years = np.arange(1861, 2100)
for year in years:
    print("ncecat -O HadGEM2_FireOff"+str(year)+".nc HadGEM2_FireOff"+str(year)+"new.nc")
    print("ncpdq -O -a time,record HadGEM2_FireOff"+str(year)+"new.nc HadGEM2_FireOff"+str(year)+"new.nc")

#ncrcat IPSL*new.nc IPSL_ALL.nc
raise Exception ('Done')


'''


def CalcWorld(cube):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    grid_weights = iris.analysis.cartography.area_weights(cube)

    var_constraint2 = iris.Constraint(cube_func=lambda x: x.var_name == 'field36')
    landmask = iris.load_cube('/hpc/data/d05/hadaw/HadGEM3-GA6/Ancils/GA7_Ancil/CRU-NCEPv7.landfrac.nc', var_constraint2)
    landmask = landmask.regrid(cube, iris.analysis.Linear())
    cube = cube*landmask

    TOT = cube.collapsed(coords, iris.analysis.SUM, weights=grid_weights)/1E12 
    if len(TOT.data) == 239:
        TOT = TOT[:-1]
    print (len(TOT.data))
    return TOT
  


def CalcRegion(mydata, model, fire):
    
    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions')
    GFED1 = iris.load_cube('/home/h01/cburton/PAPERS/3.FutureFiresGlobal/Data/GFED.nc')

    region = mydata.copy()
    GFED = GFED1.regrid(mydata, iris.analysis.Linear()) 
    
    file_dict = OrderedDict([
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
    
    # regions for each RCP
    cube_dict = {}

    for key in file_dict:
        #print (key)
        region.data = np.ma.masked_where(GFED.data != file_dict[key],mydata.data)
        coords = ('longitude', 'latitude')
        for coord in coords:
            if not region.coord(coord).has_bounds():
                region.coord(coord).guess_bounds()
        grid_weights = iris.analysis.cartography.area_weights(region)

        TOT = region.collapsed(coords, iris.analysis.SUM, weights=grid_weights)/1E12
        if len(TOT.data) == 239:
            TOT = TOT[:-1]
        cube_dict[key] = TOT

        f = open('/home/h01/cburton/PAPERS/3.FutureFiresGlobal/Data/NBP.dat','a')
        np.savetxt(f,(model,key,fire,TOT.data),newline=' ',fmt='  %s')
        f.write('\n')
        f.close()
  
    return cube_dict


def CalcTemp(Temp):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not Temp.coord(coord).has_bounds():
            Temp.coord(coord).guess_bounds()
    weights = iris.analysis.cartography.area_weights(Temp)
    TEMP = Temp.collapsed(coords, iris.analysis.MEAN, weights = weights) 
    if len(TEMP.data) == 239:
        TEMP = TEMP[:-1]
    return TEMP


models = ('H', 'G', 'I', 'M')
NBPdict = {}
for model in models:

        #Calc Global Mean Temp from driving data
        Temp = iris.load_cube('~/PAPERS/3.FutureFiresGlobal/Data/DRIVE/'+model+'_TempCAbovePI.nc')
        TEMP = CalcTemp(Temp)
        NBPdict['_TEMP_'+model] = TEMP

        f = open('/home/h01/cburton/PAPERS/3.FutureFiresGlobal/Data/NBP.dat','a')
        np.savetxt(f,(model,'Glob','Temp',TEMP.data),newline=' ',fmt='  %s')
        f.write('\n')
        f.close()

        # NBP Fire ON
        print ('starting FIRE ON')
        NBP_ON = iris.load_cube('/home/h01/cburton/PAPERS/3.FutureFiresGlobal/Data/NBPyears/'+model+'*_ALL.nc')
        print(NBP_ON)

        NBP_ON=NBP_ON[:,0,:,:]
        NBP_ON.data[np.where(np.isnan(NBP_ON.data))] = 0.0
        ### Calc world
        GLOB = CalcWorld(NBP_ON)
        NBPdict['_GlobON_'+model] = GLOB

        f = open('/home/h01/cburton/PAPERS/3.FutureFiresGlobal/Data/NBP.dat','a')
        np.savetxt(f,(model,'Glob','NBP F-ON',GLOB.data),newline=' ',fmt='  %s')
        f.write('\n')
        f.close()

        ### Calc regions
        NBPdict['_ON_'+model]= CalcRegion(NBP_ON, model, fire='ON')


        # NBP Fire OFF
        print ('starting FIRE OFF')
        NBP_OFF = iris.load_cube('/home/h01/cburton/PAPERS/3.FutureFiresGlobal/Data/NBPyears_FireOff/'+model+'*_ALL.nc')
        NBP_OFF=NBP_OFF[:,0,:,:]
        NBP_OFF.data[np.where(np.isnan(NBP_OFF.data))] = 0.0
        ### Calc world
        GLOB_OFF = CalcWorld(NBP_OFF)
        NBPdict['_GlobOFF_'+model]= GLOB_OFF

        f = open('/home/h01/cburton/PAPERS/3.FutureFiresGlobal/Data/NBP.dat','a')
        np.savetxt(f,(model,'Glob','NBP F-OFF',GLOB_OFF.data),newline=' ',fmt='  %s')
        f.write('\n')
        f.close()

        ### Calc regions
        NBPdict['_OFF_'+model] = CalcRegion(NBP_OFF, model, fire='OFF')


raise Exception ('Done')




## T = 239 points
## NBP = 238 points

#Find multi-model means
Temp = (NBPdict['_TEMP_H'].data+NBPdict['_TEMP_G'].data+NBPdict['_TEMP_I'].data+NBPdict['_TEMP_M'].data)/4.0

Glob_FireON = (NBPdict['_GlobON_H'].data+NBPdict['_GlobON_G'].data+NBPdict['_GlobON_I'].data+NBPdict['_GlobON_M'].data)/4.0
Glob_FireOFF = (NBPdict['_GlobOFF_H'].data+NBPdict['_GlobOFF_G'].data+NBPdict['_GlobOFF_I'].data+NBPdict['_GlobOFF_M'].data)/4.0

regions = ('BONA', 'TENA','CEAM','NHSA','SHSA','EURO','MIDE','NHAF','SHAF', 'BOAS', 'CEAS','SEAS','EQAS', 'AUST')
for region in regions: 
    exec(f'{region}_FireON = (NBPdict["_ON_H"][region].data+NBPdict["_ON_G"][region].data+NBPdict["_ON_I"][region].data+NBPdict["_ON_M"][region].data)/4.0')
    exec(f'{region}_FireOFF = (NBPdict["_OFF_H"][region].data+NBPdict["_OFF_G"][region].data+NBPdict["_OFF_I"][region].data+NBPdict["_OFF_M"][region].data)/4.0')


##Alternative script for saving out data into text file
models = ('H','G', 'I', 'M')
for model in models:
    with open('/scratch/cburton/ISIMIP_PAPER/KSscripts/NBP/RESULTS_NBP.txt', 'a') as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(model.split())
        #writer.writerow("Year".split())
        #writer.writerow(NBPdict['_TEMP_'+model].coord('time'))
        writer.writerow("GlobTemp".split())
        writer.writerow((NBPdict['_TEMP_'+model].data))
        #writer.writerow("Year".split())
        #writer.writerow(NBPdict['_GlobON_'+model].coord('time'))
        writer.writerow("NBPGlobFON".split())
        writer.writerow((NBPdict['_GlobON_'+model].data))
        #writer.writerow("Year".split())
        #writer.writerow(NBPdict['_GlobOFF_'+model].coord('time'))
        writer.writerow("NBPGlobFOFF".split())
        writer.writerow((NBPdict['_GlobOFF_'+model].data))
        regions = ('BONA', 'TENA','CEAM','NHSA','SHSA','EURO','MIDE','NHAF','SHAF', 'BOAS', 'CEAS','SEAS','EQAS', 'AUST')
        for region in regions:
            writer.writerow(region.split())
            #writer.writerow(NBPdict["_ON_"+model][region].coord('time'))
            writer.writerow((NBPdict["_ON_"+model][region].data))
            #writer.writerow(NBPdict["_OFF_"+model][region].coord('time'))
            writer.writerow((NBPdict["_OFF_"+model][region].data))
f.close()
raise Exception ('Done')









