

import numpy as np
import iris
 

#Find PI global mean temp 1860-1900

models = ('GFDL-ESM2M', 'IPSL-CM5A-LR', 'MIROC5')#'HADGEM2-ES',
for model in models:
    if model == 'HADGEM2-ES':
        filename = 'HadGEM2-ES'
    else:
        filename = model
    #files = '/hpc/data/d00/hadea/isimip2b/historical/'+model+'/tas_day_'+filename+'_historical_r1i1p1_EWEMBI_*187*.nc'
    files = '/scratch/cburton/ISIMIP_PAPER/DRIVE/tas_day_'+filename+'_historical_r1i1p1_EWEMBI*.nc'
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
    #FutFiles = '/hpc/data/d00/hadea/isimip2b/historical/'+model+'/tas_day_'+filename+'_historical_r1i1p1_EWEMBI_*.nc'
    FutFiles = '/scratch/cburton/ISIMIP_PAPER/DRIVE/'+model+'/Years.nc'
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

'''
Pre-Industrial
RESULTS........................
HADGEM2-ES
13.73656565270101
RESULTS........................
GFDL-ESM2M
13.585987351477565
RESULTS........................
IPSL-CM5A-LR
13.242544160527304
RESULTS........................
MIROC5
13.726835381729321
'''



'''
models = ('HADGEM2-ES', 'GFDL-ESM2M', 'IPSL-CM5A-LR', 'MIROC5')
dates = ('18610101-18701231', '18710101-18801231')

for model in models:
    for date in dates:
        if model == 'HADGEM2-ES':
            filename == 'HadGEM2-ES'
        else:
            filename == model
        
        files = '/hpc/data/d00/hadea/isimip2b/rcp85/'+model+'tas_day_'+filename+'historical_r1i1p1_EWEMBI_'+date+'.nc'
        print (files)
        cube  = iris.load_cube(files)
        print (cube)

        if date == '18610101-18701231':
            year == 1861
        else:
           year == 1871
        for year in np.arange(10):
            date = iris.Constraint(time=lambda cell: cell.point.year == year)
            with iris.FUTURE.context(cell_datetime_objects=True):
                cube = cube.extract(date)   
                cube.convert_units('celsius')
                cube = cube.collapsed(('time'), iris.analysis.MEAN) 
coords = ('longitude', 'latitude')
for coord in coords:
    if not cube.coord(coord).has_bounds():
        cube.coord(coord).guess_bounds()
weights = iris.analysis.cartography.area_weights(cube)
cube = cube.collapsed(coords, iris.analysis.MEAN, weights = weights) 
print (cube.data)
##1850-1900 GLOBAL MEAN TEMP = 13.3265186241
raise Exception ('Stop')
'''
'''    
years = np.arange(2016,2099)
for year in years:

    temp = iris.load('/scratch/cburton/FutureUKESM/Data/be397a.py'+str(year)+'1201.pp')
    temp = temp[0]
    temp.convert_units('celsius')
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not temp.coord(coord).has_bounds():
            temp.coord(coord).guess_bounds()
    weights = iris.analysis.cartography.area_weights(temp)
    temp = temp.collapsed(coords, iris.analysis.MEAN, weights = weights) 
    Anomaly = temp-13.3265186241
    print(str(temp.coord('time'))), Anomaly.data

raise Exception ('Stop')
'''
