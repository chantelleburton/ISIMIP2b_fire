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
print (GFEDregions)

cube = jules.load_cube('/scratch/cburton/scratch/ISIMIP_PAPER/HadGEM2/hadgem2-es_c20c.gen_ann_gb.1921.nc',iris.Constraint(cube_func=lambda x: x.var_name == 'tstar_gb'),missingdata=np.ma.masked).collapsed(['time'], iris.analysis.MEAN)
GFEDregions = GFEDregions.regrid(cube, iris.analysis.Linear())
print (GFEDregions)

colors = ('white','black','lightgrey','black','lightgrey','black','lightgrey','black','lightgrey','black','black','grey','lightgrey','black','lightgrey')
cmap=matplotlib.colors.ListedColormap(colors)  

#ticks on pcolormesh (for discrete colorbar)
levels = MaxNLocator(nbins=15).tick_values(0,15)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
mesh = iplt.pcolormesh(GFEDregions, cmap=cmap, norm=norm)
plt.gca().coastlines() 
plt.show()

exit()



