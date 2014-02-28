import numpy as np
from periodicGrid import PeriodicGrid
from trajectory import Trajectory
from spectralLib import specFilt 


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

    def integrate(self, ic, tInt, filtNtrc=True, t0=0.):
        """
        Call to the model propagator

            Launcher.integrate(ic, tInt, filtNtrc=True)

            ic  :   initial condition <numpy.ndarray>
            tInt:   integration time <float>

        """

        if not isinstance(ic, np.ndarray):
            raise self.LauncherError("ic <numpy.ndarray>")
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise self.LauncherError("ic.shape = (grid.N,)")
        if ic.dtype<>'float64':
            raise LauncherError('Potential loss of precision')
        if filtNtrc:
            specFilt(ic, self.grid)

        self.nDt=int(tInt/self.dt)
        
        if self.nDt*self.dt < tInt :
            self.nDt+=1

        self.tInt=self.nDt*self.dt
        
        # Initialisation
        traj=Trajectory(self.grid)
        traj.initialize(ic.copy(), self.nDt, self.dt)

        # calling the propagator
        traj=self.propagator(traj.ic, traj, t0=t0)
        if traj.getData().dtype<>'float64':
            raise LauncherError('Potential loss of precision')

        return traj

    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="####| Launcher |#######################################\n"
        output+=self.grid.__str__()
        output+="\n  dt=%-23.15E"%self.dt
        output+="\n\n  propagator\n  %s"%self.propagator
        output+="\n#######################################################\n"
        return output
