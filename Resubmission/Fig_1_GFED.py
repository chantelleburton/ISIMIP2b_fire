import sys
sys.path.append('/net/home/h03/kwilliam/other_fcm/jules_py/trunk/jules/')
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
import numpy.ma as ma
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


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



var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions')
GFEDregions = iris.load_cube('/data/cr1/cburton/GFED/GFED4.1s/GFED4.1s_1997-2016.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN)

folder = '/hpc/data/d05/cburton/jules_output/u-cf137/H*/'
cube = jules.load_cube(folder+'/hadgem2-es_c20c.gen_ann_gb.1921.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'tstar_gb'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
#GFEDregions = GFEDregions.regrid(cube, iris.analysis.Linear())
GFEDregions=GFEDregions.extract(iris.Constraint(latitude=lambda cell: (-56.0) < cell < (90.0))) #Cut South America

#var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'LSM')
#new_landmask = landmask = iris.load_cube('/data/users/hadea/jules_ancils/isimip2/landfrac_isimip_0p5.nc', var_constraint)
#landmask = landmask.regrid(GFEDregions, iris.analysis.Linear())
#GFEDregions = GFEDregions*landmask

GFEDregions.data[GFEDregions.data == 1.0] = 98.3
GFEDregions.data[GFEDregions.data == 2.0] = 96.6
GFEDregions.data[GFEDregions.data == 3.0] = 81.0
GFEDregions.data[GFEDregions.data == 4.0] = 89.7
GFEDregions.data[GFEDregions.data == 5.0] = 63.5
GFEDregions.data[GFEDregions.data == 6.0] = 45.1
GFEDregions.data[GFEDregions.data == 7.0] = 96.9
GFEDregions.data[GFEDregions.data == 8.0] = 38.7
GFEDregions.data[GFEDregions.data == 9.0] = 97.2
GFEDregions.data[GFEDregions.data == 10.0] = 92.1
GFEDregions.data[GFEDregions.data == 11.0] = 39.0
GFEDregions.data[GFEDregions.data == 12.0] = 70.1
GFEDregions.data[GFEDregions.data == 13.0] = 86.2
GFEDregions.data[GFEDregions.data == 14.0] = 67.2


print(np.min(GFEDregions.data))
GFEDregions.data = ma.masked_where(GFEDregions.data <= 0.99, GFEDregions.data)
#GFEDregions.data[GFEDregions.data <= 0.01] = 50.0



cmap='RdBu' 
mesh = iplt.pcolormesh(GFEDregions, cmap=cmap, vmin=0, vmax=100)
plt.gca().coastlines() 
plt.colorbar(orientation='horizontal', label='Percentage overlap with observed burned area trend')
plt.show()
exit()




