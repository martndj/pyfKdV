import numpy as np

from grid import Grid
from trajectory import Trajectory
from Launcher import Launcher
from gradTest import gradientTest, gradTestString



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

    def integrate(self, ic, tInt=None, t0=0.):
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
        return traj

    #-------------------------------------------------------

    def adjoint(self, final, tInt=None, t0=0.):         
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
        
        return traj

    #-------------------------------------------------------

    def d_intTimesAdj(self, d_xIn, t0=0.):
        """
        Numerical adjoint of d_intTimes.()
        """
        d_x=d_xIn.copy()
        times=sorted(d_x.keys())
        nTimes=len(times)


        # Sn*
        for i in xrange(nTimes-1,0,-1):
            t_pre=times[i-1]
            t=times[i]
            d_x[t_pre]=(self.adjoint(d_x[t], t-t_pre, t0=t_pre).ic
                            + d_x[t_pre])
            d_x[t]=0.

        
        # I0*
        if times[0]==t0:
            adj=d_x[times[0]]
        else:
            adj=(self.adjoint(d_x[times[0]], times[0]-t0, t0=t0).ic)

        return adj
        
    


    #-------------------------------------------------------

    def gradTest(self, nlModel, tInt=None, t0=0., 
                    powRange=[-1,-14], euclidNorm=False,
                    output=True):
        if not isinstance(nlModel, Launcher):
            raise TypeError()
        if nlModel.grid<>self.grid:
            raise ValueError()
        if not self.isReferenced:
            raise RuntimeError(
                        "Not initialized with a reference trajectory")
        tInt=self._timeValidation(tInt, t0)
        gradJ0=self.adjoint(self.refTraj.final, tInt, t0=t0).ic
        if euclidNorm:
            J0=0.5*np.dot(self.refTraj.final, self.refTraj.final)
            n2GradJ0=np.dot(gradJ0, gradJ0)
        else:
            J0=0.5*self.grid.squareNorm(self.refTraj.final)
            n2GradJ0=self.grid.squareNorm(gradJ0)

        test={}
        for power in xrange(powRange[0],powRange[1], -1):
            eps=10.**(power)
            Mx=nlModel.integrate(self.refTraj.whereTime(t0)-eps*gradJ0, 
                                    tInt, t0=t0).final
            if euclidNorm:
                Jeps=0.5*np.dot(Mx, Mx)
            else:
                Jeps=0.5*self.grid.squareNorm(Mx)
            res=((J0-Jeps)/(eps*n2GradJ0))
            test[power]=[Jeps, res]

        if output:  print(gradTestString(J0, n2GradJ0, test))
        return (J0, n2GradJ0, test)

    def gradTest2(self, nlModel, tInt=None, t0=0., 
                    powRange=[-1,-14], 
                    output=True):
        if not isinstance(nlModel, Launcher):
            raise TypeError()
        if nlModel.grid<>self.grid:
            raise ValueError()
        if not self.isReferenced:
            raise RuntimeError(
                        "Not initialized with a reference trajectory")
        tInt=self._timeValidation(tInt, t0)
        def fct(x0):
            return 0.5*self.grid.squareNorm(
                        nlModel.integrate(x0, tInt, t0=t0).final)
        def gradFct(x0):
            return self.adjoint(nlModel.integrate(x0, tInt, t0=t0).final,
                                tInt, t0=t0).ic
        return gradientTest(self.refTraj.ic, fct, gradFct, 
                            powRange=powRange, output=output)
        

#    def gradTestString(self, J0, n2GradJ0, test):
#        s="----| Gradient test |------------------\n"
#        s+="  J0      =%+25.15e\n"%J0
#        s+=" |grad|^2 =%+25.15e\n"%n2GradJ0
#        for i in  (np.sort(test.keys())[::-1]):
#            s+="%4d %+25.15e  %+25.15e\n"%(i, test[i][0], test[i][1])
#        return s
    

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
        self._nDt0=int((t0-self.refTraj.t0)/self.dt)
        self._nDt=int((tInt)/self.dt)
        # real times from step integers
        if self._nDt*self.dt<tInt:
            self._nDt+=1
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
