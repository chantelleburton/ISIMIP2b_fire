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
import matplotlib
from collections import OrderedDict 
import numpy as np
import iris
from scipy.stats import ks_2samp
import cftime
from datetime import datetime, timedelta
import time
from datetime import date
from time import strftime, strptime


def CalcNBP(model,year):
        NPP = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'npp_n_gb'))     
        NPP.units =  cf_units.Unit(1) 
        SoilResp =jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'resp_s_to_atmos_gb'))
        SoilResp.units =  cf_units.Unit(1) 
        WPf = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'WP_fast_out'))
        WPf.units =  cf_units.Unit(1) 
        WPm = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'WP_med_out'))
        WPm.units =  cf_units.Unit(1) 
        WPs = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'WP_slow_out') ) 
        WPs.units =  cf_units.Unit(1) 
        Hv = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'harvest_gb'))
        Hv.units =  cf_units.Unit(1) 
        
        FireV = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.gen_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'veg_c_fire_emission_gb'))
        FireV.units =  cf_units.Unit(1) 
        FireDPM = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.gen_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_carbon_dpm'))
        FireDPM.units =  cf_units.Unit(1) 
        FireRPM = jules.load_cube('/scratch/cburton/ISIMIP_PAPER/'+model+'/*.gen_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_carbon_rpm'))
        FireRPM.units =  cf_units.Unit(1) 
        
        #NBP = (NPP - SoilResp - WPf - WPm - WPs)/1000000000000

        NBP = (NPP - SoilResp - WPf - WPm - WPs - Hv - FireV - FireDPM - FireRPM)
        return NBP




models = ('GFDL', 'IPSL', 'MIROC','HadGEM2')
for model in models:
    years = np.arange(1861,2100)
    for year in years:
        NBP  = CalcNBP(model,year) 
        iris.save(NBP, '/scratch/cburton/ISIMIP_PAPER/KScubes/NBPyears/'+model+str(year)+'.nc')


print ("*******************calculation done **************************")

#https://web.yammer.com/main/org/metoffice.gov.uk/search/threads?search=p_value






raise Exception ('Done')



'''
## Notes
t_value = mean1.copy(data=None)
p_value = mean1.copy(data=None)
t_value.data, p_value.data = ks_2samp(x,y)

for i in range(0, nlat-1):
    for j in range(0, nlon-1): 
        print (i,j, p_value.data[i,j])

x = [3]
y = [216]

t_value, p_value = ks_2samp(x,y)
test = ks_2samp(x,y)
print (test)
print (t_value)
print (p_value)


raise Exception ('Done')
'''

