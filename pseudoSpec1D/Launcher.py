import numpy as np
from pseudoSpec1D import PeriodicGrid, Trajectory  


class Launcher(object):
    """
    Launcher master class

    Launcher class has one main method:
        * integrate(ic <numpy.ndarray>, tInt <float>)
    and one fundamental data which define the integral propagator
        * propagator <function>

    Its constructor can vary from one subclass another and should
    be explicitly overloaded.

    <!> This is a dummy master class (no propagator defined), only 
        defined subclasses should be instantiated. 
    """
    class LauncherError(Exception):
        pass
    
    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, grid, dt):
        
        if not isinstance(grid, PeriodicGrid):
            raise self.LauncherError("grid <PeriodicGrid>")
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
