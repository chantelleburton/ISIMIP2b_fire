#module load scitools/default-current
#python3

# -*- coding: iso-8859-1 -*-
import os
import glob
import sys
sys.path.append('~kwilliam/fcm/python')
import iris
import cf_units
import re
import numpy as np
import jules
import iris.coord_categorisation
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from iris.plot import pcolormesh
import iris.quickplot as qplt
import iris.analysis.cartography
import iris.plot as iplt
import numpy.ma as ma

def main():


    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb') 
    control_cube = jules.load_cube('/MyScratchFolder/HadGEM2/RCP60_2p0/Years.nc', var_constraint)
    control_cube = control_cube.collapsed(['time'], iris.analysis.MEAN)*86400*365

    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb') 
    pop_cube = jules.load_cube('/MyScratchFolder/HadGEM2_PDpop/RCP6p0_2deg/Years.nc', var_constraint)
    pop_cube = pop_cube.collapsed(['time'], iris.analysis.MEAN)*86400*365

    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb')
    lu_cube = jules.load_cube('/MyScratchFolder/HadGEM2_PDlu/RCP6p0_2deg/Years.nc', var_constraint)
    lu_cube = lu_cube.collapsed(['time'], iris.analysis.MEAN)*86400*365



    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb') 
    control_PD = jules.load_cube('/MyScratchFolder/HadGEM2/2010-2019/Years.nc', var_constraint)
    control_PD = control_PD.collapsed(['time'], iris.analysis.MEAN)*86400*365

    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb') 
    pop_PD = jules.load_cube('/MyScratchFolder/HadGEM2_PDpop/2010-2019/Years.nc', var_constraint)
    pop_PD = pop_PD.collapsed(['time'], iris.analysis.MEAN)*86400*365

    var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'burnt_area_gb')
    lu_PD = jules.load_cube('/MyScratchFolder/HadGEM2_PDlu/2010-2019/Years.nc', var_constraint)
    lu_PD = lu_PD.collapsed(['time'], iris.analysis.MEAN)*86400*365


    control = control_cube - control_PD
    pop = pop_cube - pop_PD
    lu = lu_cube - lu_PD
   

    plt.subplot(1,3,1)
    iplt.pcolormesh(control, cmap='bwr', vmin=-0.1, vmax=0.1)
    plt.colorbar(orientation='horizontal')
    plt.gca().coastlines()
    plt.title('a) Transient')

    plt.subplot(1,3,2)
    iplt.pcolormesh(pop, cmap='bwr', vmin=-0.1, vmax=0.1)
    plt.colorbar(orientation='horizontal')
    plt.gca().coastlines()
    plt.title('b) Present day population')

    plt.subplot(1,3,3)
    iplt.pcolormesh(lu, cmap='bwr', vmin=-0.1, vmax=0.1)
    plt.colorbar(orientation='horizontal')
    plt.gca().coastlines()
    plt.title('c) Present day land-use')


    plt.suptitle('Change in Burnt Area')
    plt.show()


if __name__ == '__main__':
    main()


   



