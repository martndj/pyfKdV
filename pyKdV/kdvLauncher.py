import numpy as np

from pseudoSpec1D import *
from kdvParam import *

import fKdV

class kdvLauncher(Launcher):
    """
    Launcher subclass for Augmented Korteweg-de Vries system

    kdvLauncher(param, maxA, dtMod=0.7)
        
        param   :   system local parametrisation <kdvParam>
        maxA    :   maximum signal amplitude <float>
                        for calculating stable time increment
        dtMod   :   explicit modificator to neutraly stable time
                        increment
    """
    class kdvLauncherError(Exception):
        pass
    
    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, maxA, dtMod=0.7):

        if not (isinstance(param, Param)):
            raise self.kdvLauncherError(
                  "param must be an instance of Param")

        self.grid=param.grid
        self.param=param
        
        self.dtMod=dtMod
        self.dt=self.dtMod*self.dtStable(maxA)
        self.maxA=maxA

        self.propagator=self.__kdvProp_Fortran

        
    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------

    def integrate(self, ic, tInt, filtNtrc=False):
    
        # the Fortran propagator filter implicitly
        return super(kdvLauncher, self).integrate(ic, tInt,
                                                filtNtrc)

    #------------------------------------------------------

    def dtStable(self, maxA):
    
        maxK=2.0*np.pi*self.grid.Ntrc/self.grid.L
        denom=np.zeros(shape=self.grid.N)
        denom=np.sqrt((self.param[3]*maxK**3-self.param[1]*maxK
                            -self.param[2]*maxA*maxK)**2
                       +self.param[4]**2)

        dt=1./denom.max()
        return dt

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
    grid=PeriodicGrid(150,300.)
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
    launcher=kdvLauncher(param, maxA)
    
    traj=launcher.integrate(ic, tInt)
    axe=traj.waterfall()
    plt.show()
