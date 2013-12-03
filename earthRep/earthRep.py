from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import numpy as np
from pseudoSpec1D import PeriodicGrid
#-----------------------------------------------------------

def grid2Sphere(v, grid, m, lon_0=-107, jetLat=50., amplMax=5.):
    nLons=grid.N
    lons = np.linspace(-180.+lon_0, 180+lon_0, nLons)
    lats = np.ones(nLons)*jetLat +v*amplMax
    x,y= m(lons, lats)
    return x,y

    

def ortho(v, grid,  
        lon_0=-107, lat_0=50, resDpi=200, projection='ortho', 
        resolution='l', jetLat=45., amplMax=5., coord=True, **kwargs):

    if not isinstance(grid, PeriodicGrid):
        raise TypeError()

    m = Basemap(projection=projection,lon_0=lon_0,lat_0=lat_0,
                resolution=resolution)

    x,y=grid2Sphere(v, grid, m, lon_0=lon_0, jetLat=45., amplMax=amplMax)

    m.plot(x,y, **kwargs)
    if coord:
        m.drawmeridian(np.linspace(0,350, 18))
        m.drawparallel(np.linspace(-80,80, 8))
    return m 



def overCanada(v, grid, amplMax=5., **kwargs):
    m = Basemap(width=12000000,height=9000000,projection='lcc',
            resolution='c',lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.)

    if not isinstance(grid, PeriodicGrid):
        raise TypeError()

    x,y=grid2Sphere(v, grid, m, lon_0=-107, jetLat=45., amplMax=amplMax)

    m.plot(x, y, **kwargs)
    return m
