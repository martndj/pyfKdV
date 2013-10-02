import numpy as np
from pseudoSpec1D import PeriodicGrid, Trajectory, specFilt 


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

    def integrate(self, ic, tInt, filtNtrc=True):

        if not isinstance(ic, np.ndarray):
            raise self.LauncherError("ic <numpy.ndarray>")
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise self.LauncherError("ic.shape = (grid.N,)")
        if ic.dtype<>'float64':
            raise LauncherError('Potential loss of precision')
        if filtNtrc:
            specFilt(ic, self.grid)

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
        if traj.getData().dtype<>'float64':
            raise LauncherError('Potential loss of precision')

        return traj

    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="====| Launcher |===================================\n"
        output+=self.grid.__str__()
        output+="\n| dt=%-23.15E"%self.dt
        output+="\n| propagator=%s"%self.propagator
        output+="\n===================================================\n"
        return output
