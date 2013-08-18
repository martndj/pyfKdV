import numpy as np

from pseudoSpec1D import *
from kdvParam import *

import fKdV

class Launcher(object):
    """
    """
    class LauncherError(Exception):
        pass
    
    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, ic):

        if not (isinstance(param, Param)):
            raise LauncherError(
                  "param must be an instance of Param")

        self.grid=param.grid
        self.param=param
        self.ic=ic
        

        self.propagator=self.__kdvProp_Fortran

        
    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------

    def dtStable(self, maxA):
    
        maxK=2.0*np.pi*self.grid.Ntrc/self.grid.L
        denom=np.zeros(shape=self.grid.N)
        denom=np.sqrt((self.param[3]*maxK**3-self.param[1]*maxK
                       -self.param[2]*maxA*maxK)**2+self.param[4]**2)

        dt=1./denom.max()
        return dt


    #------------------------------------------------------

    def integrate(self, tInt, maxA, dtMod=0.7):
        dt=dtMod*self.dtStable(maxA)

        # Initialisation
        traj=Trajectory(self.grid)
        traj.initialize(self.ic, tInt, dt)

        # calling the propagator
        traj=self.propagator(traj)

        return traj
    

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------

    def __kdvProp_Fortran(self, traj):
        
        # Local variables names
        grid=self.grid
        param=self.param
        tReal=0.

        trajData=fKdV.fKdVPropagator(
                    grid.N, grid.Ntrc, grid.L, traj.dt, traj.nDt,
                    self.ic, param[1], param[2], param[3], param[4],
                    param[0]
                    )

        tReal=traj.nDt*traj.dt
        traj.putData(trajData)
        traj.incrmTReal(finished=True, tReal=tReal)

        return traj


#====================================================================
#--------------------------------------------------------------------
#====================================================================

if __name__=='__main__':
    
    grid=SpectralGrid(150,300.)
    tInt=50.
    maxA=2.

    def funcNulle(x):
        return 0.
    def funcNeg(x):
        return -1.
    def funcPos(x):
        return 1.
    def gauss(x):
        return -0.1*np.exp(-(x**2)/100.)
    def sinus(x):
        return 0.1*np.sin(2.*2*np.pi*x/150.)

    param=Param(grid, sinus, funcNulle, funcPos, funcNeg, gauss)

    def soliton(grid, x0, amp, alpha, beta, gamma):
        """ 
        Eigen solution of the  system
    
        """
        N=grid.N
        fct=np.zeros(N, dtype=np.float64)
    
        sgn=beta*gamma/np.abs(beta*gamma)
        amp=np.abs(amp)*sgn
        discr=np.sqrt(amp*beta/(12.*gamma))*(grid.x-x0)
        fct=amp*np.cosh(discr)**(-2)
    
        return fct 
    ic=soliton(grid, 0., 1., 0., 1., -1.)

    # NL model integration
    launcher=Launcher(param, ic)
    
    traj=launcher.integrate(tInt, maxA)
    axe=traj.waterfall()
    plt.show()
