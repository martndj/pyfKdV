import numpy as np

from periodicGrid import PeriodicGrid
from trajectory import Trajectory
from spectralLib import specFilt
from Launcher import Launcher



class TLMLauncher(object):
    """
    TLMLauncher master class

    TLMLauncher class has three main methods:
        * initialize(traj <Trajectory>)
        * integrate(ic <numpy.ndarray>, tInt <float>)
        * adjoint(fi <numpy.ndarray>, tInt <float>)
    and  fundamental data which define the integral propagator
        * propagator <function>
        * propagatorAdj <function>

    Its constructor can vary from one subclass another and should
    be explicitly overloaded.

    <!> Before integration (direct or adjoint), the TLMLauncher
        must be initialized with a reference trajectory.
    
    <!> This is a dummy master class (no propagator defined), only 
        defined subclasses should be instantiated. 
    """
    class TLMLauncherError(Exception):
        pass

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, grid, traj=None):

        if not isinstance(grid, PeriodicGrid):
            raise self.TLMLauncherError("grid <PeriodicGrid>")
        self.grid=grid

        self.isReferenced=False 
        if not traj==None: self.initialize(traj)

        # Status Attributes
        self.isIntegrated=False
        self.fullPertTraj=False
        self.tReal=0.

        self.propagator=None
        self.propagatorAdj=None

    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------
    
    def reference(self, traj):
        """
        Initialize the TLM by associating a reference trajectory
            (mandatory before integration)

            TLMLauncher.initialize(traj)

            traj    :   reference trajectory <Trajectory>
            
        """
        if not(isinstance(traj, Trajectory)):
            raise self.TLMLauncherError("traj <Trajectory>")
        if not (traj.grid==self.grid):
            raise PeriodicGridError("traj.grid <> grid")
        self.refTraj=traj
        self.isReferenced=True 
    
    def initialize(self, traj):
        """
        Alias of reference (for retro-compatibility)
        """
        return self.reference(traj)
    #------------------------------------------------------
    
    def incrmTReal(self, finished=False, tReal=None):
        if tReal<>None:
            self.tReal=tReal
        if finished:
            self.isIntegrated=True
        else:
            self.tReal+=self.refTraj.dt

    #-------------------------------------------------------

    def integrate(self, pert, tInt=None, t0=0., fullPertTraj=False,
                    filtNtrc=True):
        """
        Call to the TLM propagator

            TLMLauncher.integrate(pert, tInt, filtNtrc=True)

            pert    :   initial perturbation <numpy.ndarray>
            tInt    :   integration time <float>
            t0      :   initial time <float>
        """
        
        if not self.isReferenced:
            raise self.TLMLauncherError(
                        "Not initialized with a reference trajectory")
        if not isinstance(pert, np.ndarray):
            raise self.TLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise self.TLMLauncherError("pert.shape = (launcher.grid.N,)")
        if pert.dtype<>'float64':
            raise LauncherError('Potential loss of precision')

        if filtNtrc:
            specFilt(pert, self.grid)
        self._timeValidation(tInt, t0)

        if fullPertTraj:
            self.fullPertTraj=True
            self.pertTraj=Trajectory(self.grid)
            self.pertTraj.initialize(pert, self.nDt, self.dt)
        else:
            self.fullPertTraj=False

        fPert=self.propagator(pert.copy(), t0=t0)

        if fPert.dtype<>'float64':
            raise LauncherError('Potential loss of precision')
        return fPert
                                            

    #-------------------------------------------------------

    def adjoint(self, pert, tInt=None, t0=0., filtNtrc=True):         
        """
        Call to the TLM Adjoint retro-propagator

            TLMLauncher.integrate(pert, tInt, filtNtrc=True)

            pert    :   initial perturbation <numpy.ndarray>
            tInt    :   integration time <float>
            t0      :   initial time <float>
        """
        
        if not self.isReferenced:
            raise self.TLMLauncherError(
                        "Not initialized with a reference trajectory")
        if not isinstance(pert, np.ndarray):
            raise self.TLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise self.TLMLauncherError("pert.shape = (launcher.grid.N,)")
        if pert.dtype<>'float64':
            raise LauncherError('Potential loss of precision')

        if filtNtrc:
            specFilt(pert, self.grid)
        self._timeValidation(tInt, t0)

        adj= self.propagatorAdj(pert.copy())
        if adj.dtype<>'float64':
            raise LauncherError('Potential loss of precision')
        return adj
    

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------
    
    def _timeValidation(self, tInt, t0):

        # Time attributes
        if tInt==None : tInt=self.refTraj.tReal
        if t0<0.:
            raise self.TLMLauncherError("t0>=0.")
        elif t0<self.refTraj.t0:    
            raise self.TLMLauncherError("t0>=self.refTraj.t0")
        if t0+tInt>self.refTraj.t0+self.refTraj.tReal:
            raise self.TLMLauncherError("%.2f %.2f %.2f"%(t0,tInt,self.refTraj.tReal))
        self.dt=self.refTraj.dt
        self.nDt0=int((t0-self.refTraj.t0)/self.dt)
        self.nDt=int((tInt)/self.dt)
        # real times from step integers
        if self.nDt*self.dt<tInt:
            self.nDt+=1
        self.nDtFinal=self.nDt0+self.nDt
        self.t0=self.nDt0*self.dt
        self.tf=self.nDtFinal*self.dt+t0

    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="####| TLMLauncher |####################################\n"
        output+=self.grid.__str__()
        if not self.isReferenced:
            output+="\n  not initialized with a reference trajectory"
        else:
            output+="\n  reference trajectory:\n"
            output+=self.refTraj.__str__()
        output+="\n\n  propagator\n  %s"%self.propagator
        output+="\n\n  adjoint propagator\n  %s"%self.propagatorAdj
        output+="\n#######################################################\n"
        return output
