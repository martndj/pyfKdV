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
    
    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, dt=None, maxA=3.):

        if not (isinstance(param, Param)):
            raise TypeError(
                  "param must be an instance of Param")

        self.grid=param.grid
        self.param=param
        self.isTimeDependant=self.param.isTimeDependant
        
        self.maxA=maxA
        if dt==None:
            self.dt=self.dtStable(self.maxA)
        else:
            self.dt=dt
            if dt>self.dtStable(self.maxA):
                raise RuntimeError(
                    "dt too small: potential numerical instability")

        self.propagator=self.__kdvProp_Fortran

        
    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------


    def dtStable(self, maxA, dtMod=1.):
        """
        Stable time incremement

            kdvLauncher.dtStable(maxA)

            maxA    :   expected maximum amplitude <float>
        """
        return dtStable(self.param, maxA, dtMod=dtMod)

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------


    def __kdvProp_Fortran(self, ic, traj, t0=0.):
        
        # Local variables names
        grid=self.grid

        # to be corrected (?: what is the problem?)
        if (not (self.param.nDt==0 and self.param.dt==0.)
            and (self.dt <> self.param.dt)):
            raise RuntimeError(
                    "incompatible parameter time increment (%f, %f)"%(
                                                        self.dt, param.dt))

        # if t0<>0 param must be adjusted before passed to propagator
        param=self.__t0AdjustParam(self.param, t0=t0)
                
        

        trajData=fKdV.fKdVPropagator(
                    grid.N, grid.Ntrc, grid.L, self.dt, self._nDt,
                    param.nDt,
                    ic, param[1].getData(), param[2].getData(),
                    param[3].getData(), param[4].getData(),
                    param.nu, param.nuN, param[0].getData()
                    )

        tReal=traj.nDt*traj.dt
        traj.putData(trajData)
        traj.incrmTReal(finished=True, tReal=tReal, t0=t0)

        return traj

    #------------------------------------------------------

    def __t0AdjustParam(self, param, t0=0., limit=False):
        
        if param.isTimeDependant:
            if t0==param.t0 :
                return param
            elif t0>param.t0:
                if t0>=param.tf:
                    if limit:
                        raise ValueError()
                    else:
                        return param.final
                else:
                    return param.cut(t0) 
        else:
            return param

#--------------------------------------------------------------------
#====================================================================
#--------------------------------------------------------------------
if __name__=='__main__':
    import matplotlib.pyplot as plt
    from kdvMisc import *
    
    Ntrc=42
    maxA=10.
    tInt=10.
    grid=PeriodicGrid(Ntrc)

    paramNu=Param(grid, nu=10., nuN=Ntrc/4)
    param=Param(grid)
    
    dt=0.01
    if dt > dtStable(param, maxA): raise RuntimeError()

    modelNu=kdvLauncher(paramNu, dt=dt)
    model=kdvLauncher(param, dt=dt)
    
    x0=soliton(grid.x,0.)
    trajNu=modelNu.integrate(x0, tInt)
    traj=model.integrate(x0, tInt)

    clf()
    ax1=plt.subplot(311)
    ax2=plt.subplot(312)
    ax3=plt.subplot(313)
    grid.plotAll(traj.ic, axeD=ax1, axeS=ax2, label='IC')
    grid.plotAll(trajNu.final, axeD=ax1, axeS=ax2, color='r', 
                 label=r'$\nu=%.3f\ \ N_\nu=%d$'%(paramNu.nu, paramNu.nuN))
    grid.plotAll(traj.final, axeD=ax1, axeS=ax2, color='g', label=r'$\nu=0$')

    grid.plotPSpec(trajNu.final-traj.final, axe=ax3, color='r',
                   label=r'$x_f^\nu-x_f^0$')
    ax2.legend()
    ax3.legend()
