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
        NPP = jules.load_cube('/MyScratchFolder/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'npp_n_gb'))     
        NPP.units =  cf_units.Unit(1) 
        SoilResp =jules.load_cube('/MyScratchFolder/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'resp_s_to_atmos_gb'))
        SoilResp.units =  cf_units.Unit(1) 
        WPf = jules.load_cube('/MyScratchFolder/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'WP_fast_out'))
        WPf.units =  cf_units.Unit(1) 
        WPm = jules.load_cube('/MyScratchFolder/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'WP_med_out'))
        WPm.units =  cf_units.Unit(1) 
        WPs = jules.load_cube('/MyScratchFolder/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'WP_slow_out') ) 
        WPs.units =  cf_units.Unit(1) 
        Hv = jules.load_cube('/MyScratchFolder/'+model+'/*.c_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'harvest_gb'))
        Hv.units =  cf_units.Unit(1) 
        
        FireV = jules.load_cube('/MyScratchFolder/'+model+'/*.gen_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'veg_c_fire_emission_gb'))
        FireV.units =  cf_units.Unit(1) 
        FireDPM = jules.load_cube('/MyScratchFolder/'+model+'/*.gen_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_carbon_dpm'))
        FireDPM.units =  cf_units.Unit(1) 
        FireRPM = jules.load_cube('/MyScratchFolder/'+model+'/*.gen_ann_gb.'+str(year)+'.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_carbon_rpm'))
        FireRPM.units =  cf_units.Unit(1) 
        
        #NBP = (NPP - SoilResp - WPf - WPm - WPs)/1000000000000

        NBP = (NPP - SoilResp - WPf - WPm - WPs - Hv - FireV - FireDPM - FireRPM)
        return NBP




models = ('GFDL', 'IPSL', 'MIROC','HadGEM2')
for model in models:
    years = np.arange(1861,2100)
    for year in years:
        NBP  = CalcNBP(model,year) 
        iris.save(NBP, '/MyScratchFolder/NBPyears/'+model+str(year)+'.nc')


print ("*******************calculation done **************************")






raise Exception ('Done')



