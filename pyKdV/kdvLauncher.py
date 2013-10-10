import numpy as np

from pseudoSpec1D import *
from kdvParam import *
from kdvMisc import dtStable

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

    def __init__(self, param, dt=None, maxA=3.):

        if not (isinstance(param, Param)):
            raise self.kdvLauncherError(
                  "param must be an instance of Param")

        self.grid=param.grid
        self.param=param
        self.isTimeDependant=self.param.isTimeDependant
        
        self.maxA=maxA
        if dt==None:
            self.dt=self.dtStable(self.maxA)
        else:
            self.dt=dt

        self.propagator=self.__kdvProp_Fortran

        
    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------

    def integrate(self, ic, tInt, filtNtrc=False):
        """
        Call to the model propagator

            Launcher.integrate(ic, tInt, filtNtrc=True)

            ic  :   initial condition <numpy.ndarray>
            tInt:   integration time <float>

        """
    
        # the Fortran propagator filter implicitly
        return super(kdvLauncher, self).integrate(ic, tInt,
                                                filtNtrc)

    #------------------------------------------------------

    def dtStable(self, maxA, dtMod=0.7):
        """
        Stable time incremement

            kdvLauncher.dtStable(maxA)

            maxA    :   expected maximum amplitude <float>
        """
        return dtStable(self.grid, self.param, maxA, dtMod=dtMod)

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------


    def __kdvProp_Fortran(self, ic, traj):
        
        # Local variables names
        grid=self.grid
        param=self.param
        tReal=0.

        # to be corrected
        if (not (param.nDt==0 and param.dt==0.)
            and (self.dt <> param.dt)):
            raise self.kdvLauncherError(
                    "incompatible parameter time increment (%f, %f)"%(
                                                        self.dt, param.dt))
        
        trajData=fKdV.fKdVPropagator(
                    grid.N, grid.Ntrc, grid.L, self.dt, self.nDt,
                    param.nDt,
                    ic, param[1].getData(), param[2].getData(),
                    param[3].getData(), param[4].getData(),
                    param[0].getData()
                    )

        tReal=traj.nDt*traj.dt
        traj.putData(trajData)
        traj.incrmTReal(finished=True, tReal=tReal)

        return traj

#====================================================================
#--------------------------------------------------------------------
#====================================================================

if __name__=='__main__':
    
    from kdvMisc import gauss, soliton, dtStable
    grid=PeriodicGrid(150,300.)
    tInt=50.
    maxA=2.
    
    def gaussNeg(x,t):
        x0=0.
        sig=5.
        return -0.1*gauss(x, x0, sig) 
    def sinus(x,t):
        return 0.1*np.sin(2.*2*np.pi*x/150.)

    def funcTimeDependant(x, t):
        return 0.1*np.sin(x/50.)*np.cos(t/10.)

    dt=dtStable(grid, Param(grid, beta=1.,gamma=-1., rho=0.1),
                    maxA, dtMod=0.7)
    nDtParam=int(tInt/dt)


    param=Param(grid, beta=1., gamma=-1., rho=gaussNeg,
                forcing=funcTimeDependant, 
                nDt=nDtParam, dt=dt)
    paramStatic=Param(grid, beta=1., gamma=-1.)

    ic=soliton(grid.x, 1., beta=1., gamma=-1. )

    # NL model integration
    launcher=kdvLauncher(param, maxA, dt=dt)
    
    traj=launcher.integrate(ic, tInt)
    axe=traj.waterfall()
    plt.show()
