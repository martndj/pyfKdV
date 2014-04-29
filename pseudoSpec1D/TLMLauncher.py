import numpy as np

from grid import Grid
from trajectory import Trajectory
from Launcher import Launcher



class TLMLauncher(Launcher):
    """
    TLMLauncher sub-class

    TLMLauncher class has three main methods:
        * reference(traj <Trajectory>)
        * integrate(ic <numpy.ndarray>, tInt <float>)
        * adjoint(fi <numpy.ndarray>, tInt <float>)

    and two fundamental data which define the integral propagator and
    adjoint
        * propagator <function>
        * propagatorAdj <function>

    Its constructor can vary from one subclass another and should
    be explicitly overloaded.

    <!> Before integration (direct or adjoint), the TLMLauncher
        must be referenced with a trajectory.
    
    <!> This is a dummy master class (no propagator defined), only 
        defined subclasses should be instantiated. 
    """

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, grid, traj=None):

        if not isinstance(grid, Grid):
            raise ValueError()
        self.grid=grid

        self.isReferenced=False 
        if not traj==None: self.reference(traj)

        # Status Attributes
        self.isIntegrated=False
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

            TLMLauncher.reference(traj)

            traj    :   reference trajectory <Trajectory>
            
        """
        if not(isinstance(traj, Trajectory)):
            raise TypeError() 
        if not (traj.grid==self.grid):
            raise ValueError()
        self.refTraj=traj
        self.isReferenced=True 
    
    #-------------------------------------------------------

    def integrate(self, ic, tInt=None, t0=0., fullTraj=False):
        """
        Call to the TLM propagator

            TLMLauncher.integrate(pert, tInt)

            ic      :   initial perturbation <numpy.ndarray>
            tInt    :   integration time <float>
            t0      :   initial time <float>
        """
        
        if not self.isReferenced:
            raise RuntimeError(
                        "Not initialized with a reference trajectory")

        tInt=self._timeValidation(tInt, t0)
        traj=super(TLMLauncher, self).integrate(ic, tInt, 
                    propagator=self.propagator)
        if not fullTraj:
            return traj.final
        else:
            return traj

    #-------------------------------------------------------

    def adjoint(self, final, tInt=None, t0=0., fullTraj=False):         
        """
        Call to the TLM Adjoint retro-propagator

            TLMLauncher.integrate(pert, tInt)

            final   :   final perturbation <numpy.ndarray>
            tInt    :   integration time <float>
            t0      :   initial time <float>
        """
        
        if not self.isReferenced:
            raise RuntimeError(
                        "Not initialized with a reference trajectory")

        tInt=self._timeValidation(tInt, t0)
        traj=super(TLMLauncher, self).integrate(final, tInt, 
                    propagator=self.propagatorAdj)
        
        if not fullTraj:
            return traj.ic
        else:
            return traj
        # verifier que c'est bien .ic et non .final

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------
    
    def _timeValidation(self, tInt, t0):

        # Time attributes
        if tInt==None : tInt=self.refTraj.tReal
        if t0<0.:
            raise ValueError("t0>=0.")
        elif t0<self.refTraj.t0:    
            raise ValueError("t0>=self.refTraj.t0")
        if t0+tInt>self.refTraj.t0+self.refTraj.tReal:
            raise ValueError("%.2f %.2f %.2f"%(t0,tInt,self.refTraj.tReal))
        self.dt=self.refTraj.dt
        self.nDt0=int((t0-self.refTraj.t0)/self.dt)
        self.nDt=int((tInt)/self.dt)
        # real times from step integers
        if self.nDt*self.dt<tInt:
            self.nDt+=1
        self.nDtFinal=self.nDt0+self.nDt
        self.t0=self.nDt0*self.dt
        self.tf=self.nDtFinal*self.dt+t0
        return tInt

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
