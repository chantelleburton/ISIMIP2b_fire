

import numpy as np
import iris
 

#Find PI global mean temp 1860-1900

models = ('GFDL-ESM2M', 'IPSL-CM5A-LR', 'MIROC5')#'HADGEM2-ES',
for model in models:
    if model == 'HADGEM2-ES':
        filename = 'HadGEM2-ES'
    else:
        filename = model
    files = '/MyScratchFolder/DRIVE/tas_day_'+filename+'_historical_r1i1p1_EWEMBI*.nc'
    cubelist  = iris.load(files)
    iris.util.unify_time_units(cubelist)
    cube = cubelist.concatenate_cube()
    cube.convert_units('celsius')
    cube = cube.collapsed(('time'), iris.analysis.MEAN) 
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    weights = iris.analysis.cartography.area_weights(cube)
    PIcube = cube.collapsed(coords, iris.analysis.MEAN, weights = weights) 
    with open('RESULTS_HIST.txt', 'a') as f:
        f.write(model)
        f.write('\n')
        f.write('Historical')
        f.write('\n')
        f.write('PI ')
        f.write(str(PIcube.data))
        f.write('\n')

#Find GWT
    FutFiles = '/MyScratchFolder/DRIVE/'+model+'/Years.nc'
    FutCube  = iris.load_cube(FutFiles)
    print (FutCube)
    #iris.util.unify_time_units(FutCube)
    #FutCube = FutCube.concatenate_cube()
    FutCube.convert_units('celsius')
    years = np.arange(1860,2080)
    for startyear in years:
        endyear = startyear+20
        daterange = iris.Constraint(time=lambda cell: startyear <= cell.point.year < endyear)#start at 1870 (1860-1880)
        SliceCube = FutCube.extract(daterange)
        SliceCube = SliceCube.collapsed(('time'), iris.analysis.MEAN)
        coords = ('longitude', 'latitude')
        for coord in coords:
            if not SliceCube.coord(coord).has_bounds():
                SliceCube.coord(coord).guess_bounds()
        weights = iris.analysis.cartography.area_weights(SliceCube)
        SliceCube = SliceCube.collapsed(coords, iris.analysis.MEAN, weights = weights) 
        Anomaly = SliceCube - PIcube
        midpoint = startyear+10

        with open('RESULTS_HIST.txt', 'a') as f:
            f.write(str(midpoint))
            f.write(' ')
            f.write(str(Anomaly.data))
            f.write('\n')

       # print (midpoint)
       # print (Anomaly.data)


raise Exception ('Done')



