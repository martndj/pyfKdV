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

    def __init__(self, tInt, param, maxA, dtMod=0.7):

        if not (isinstance(param, Param)):
            raise self.LauncherError(
                  "param must be an instance of Param")

        self.grid=param.grid
        self.param=param
        
        self.dtMod=dtMod
        self.dt=self.dtMod*self.dtStable(maxA)
        self.tIntIn=tInt
        self.nDt=int(self.tIntIn/self.dt)
        self.tInt=self.nDt*self.dt
        if self.tInt<self.tIntIn:
            self.nDt+=1
            self.tInt=self.nDt*self.dt
        self.maxA=maxA

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

    def integrate(self, ic):

        if not isinstance(ic, np.ndarray):
            raise self.LauncherError("ic <numpy.ndarray>")
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise self.LauncherError("ic.shape = (grid.N,)")
        
        # Initialisation
        traj=Trajectory(self.grid)
        traj.initialize(ic, self.nDt, self.dt)

        # calling the propagator
        traj=self.propagator(ic, traj)

        return traj
    

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------

    def __kdvProp_Fortran(self, ic, traj):
        
        # Local variables names
        grid=self.grid
        param=self.param
        tReal=0.

        trajData=fKdV.fKdVPropagator(
                    grid.N, grid.Ntrc, grid.L, self.dt, self.nDt,
                    ic, param[1], param[2], param[3], param[4],
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
    
    from kdvMisc import gauss, soliton
    grid=SpectralGrid(150,300.)
    tInt=50.
    maxA=2.
    
    def gaussNeg(x):
        x0=0.
        sig=5.
        return -0.1*gauss(x, x0, sig) 
    def sinus(x):
        return 0.1*np.sin(2.*2*np.pi*x/150.)

    param=Param(grid, beta=1., gamma=-1., rho=gaussNeg, forcing=sinus)
    ic=soliton(grid.x, 1., beta=1., gamma=-1. )

    # NL model integration
    launcher=Launcher(tInt, param, maxA)
    
    traj=launcher.integrate(ic)
    axe=traj.waterfall()
    plt.show()
