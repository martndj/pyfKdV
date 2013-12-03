from mpl_toolkits.basemap import Basemap, maskoceans, interp
import matplotlib.cm as cm
import numpy as np
from pseudoSpec1D import PeriodicGrid, bumpSchwartz
#-----------------------------------------------------------

def earthRepresentation(v, grid,  
        lon_0=-100, lat_0=65, resDpi=200, projection='ortho', 
        resolution='l', jetLat=45., jetWidth=10., amplMax=10.,
        bluemarble=True):

    if not isinstance(grid, PeriodicGrid):
        raise TypeError()

    bump=np.vectorize(bumpSchwartz)

    m = Basemap(projection=projection,lon_0=lon_0,lat_0=lat_0,
                resolution=resolution)


    # make up some data on a regular lat/lon grid.
    nLons=grid.N
    nLats=nLons/2
    delta = 2.*np.pi/(nLons-1)

    lats = (0.5*np.pi-delta*np.indices((nLats,nLons))[0,:,:])
    lons = (delta*np.indices((nLats,nLons))[1,:,:])
    x, y = m(lons*180./np.pi, lats*180./np.pi)

    jetCoord=np.array([jetLat, jetWidth])/180.*np.pi 


    pert=v/grid.norm(v, metric=np.infty)*amplMax/180.*np.pi

    jetCoreMod=jetCoord[0]+pert
    longWave = bump((lats-jetCoreMod)/jetCoord[1])


    m.contour(x,y,longWave,15,linewidths=1.5, cmap=cm.autumn)
    if bluemarble : m.bluemarble()

