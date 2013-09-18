import numpy as np

from pseudoSpec1D import SpectralGrid, Trajectory  


class Launcher(object):
    """
    """
    class LauncherError(Exception):
        pass
    
    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, grid, dt):
        
        if not isinstance(grid, SpectralGrid):
            raise self.TLMLauncherError("grid <SpectralGrid>")
        self.grid=grid

        self.dt=dt

        self.propagator=None

    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------

    def integrate(self, ic, tInt):

        if not isinstance(ic, np.ndarray):
            raise self.LauncherError("ic <numpy.ndarray>")
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise self.LauncherError("ic.shape = (grid.N,)")
        self.tIntIn=tInt
        self.nDt=int(self.tIntIn/self.dt)
        self.tInt=self.nDt*self.dt
        if self.tInt<self.tIntIn:
            self.nDt+=1
            self.tInt=self.nDt*self.dt
        
        # Initialisation
        traj=Trajectory(self.grid)
        traj.initialize(ic, self.nDt, self.dt)

        # calling the propagator
        traj=self.propagator(ic, traj)

        return traj
