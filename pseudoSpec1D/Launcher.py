import numpy as np
from grid import Grid
from trajectory import Trajectory


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
    
    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, grid, dt):
        
        if not isinstance(grid, Grid):
            raise TypeError("grid <Grid>")
            
        self.grid=grid

        self.dt=dt

        self.propagator=None

    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------

    def integrate(self, ic, tInt, t0=0., propagator=None):
        """
        Call to the model propagator

            Launcher.integrate(ic, tInt)

            ic  :   initial condition <numpy.ndarray>
            tInt:   integration time <float>

        """

        if not isinstance(ic, np.ndarray):
            raise TypeError("ic <numpy.ndarray>")
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise ValueError("ic.shape = (grid.N,)")
        if ic.dtype<>'float64':
            raise ValueError('Potential loss of precision')

        self._nDt=int(tInt/self.dt)
        
        if self._nDt*self.dt < tInt :
            self._nDt+=1

        self._tInt=self._nDt*self.dt
        
        # Initialisation
        traj=Trajectory(self.grid)
        traj.initialize(ic.copy(), self._nDt, self.dt)

        # calling the propagator
        if propagator==None:
            traj=self.propagator(traj.ic, traj, t0=t0)
        else:
            traj=propagator(traj.ic, traj, t0=t0)
            
        if traj.getData().dtype<>'float64':
            raise RuntimeError('Potential loss of precision')

        return traj

    #-------------------------------------------------------


    def d_nDtInt(self, ic, nDtList, t0=0.):
        """ 
        Returns a dict with integrations at time steps requested
        (usefull for data assimilation model equivalent calculation)
        """
        if not isinstance(nDtList, (np.ndarray, list, set)):
            raise TypeError()
        nDtList=list(set(nDtList))
        nDtList.sort()
       
        d_x={}
        traj=self.integrate(ic, nDtList[-1]*self.dt, t0=t0)
        for i in nDtList:
            d_x[i]=traj[i]

        return d_x

    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="####| Launcher |#######################################\n"
        output+=self.grid.__str__()
        output+="\n  dt=%-23.15E"%self.dt
        output+="\n\n  propagator\n  %s"%self.propagator
        output+="\n#######################################################\n"
        return output
