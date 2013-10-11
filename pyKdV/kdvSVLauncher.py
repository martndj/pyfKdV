import numpy as np

from pseudoSpec1D import *
from kdvParam import *
from kdvLauncher import *
from kdvTLMLauncher import *

import fKdV

class kdvSVLauncher(object):
    """
    Singular vector calculator launcher class
    for Augmented Korteweg-de Vries system

    kdvSVLauncher(traj, param)

    To launch Lanczos calculation (lenghty process):

        kdvSVLauncher(traj, param).lanczos(Nev)

        traj    :   reference trajectory <Trajectory>
        param   :   parametrisation <kdvParam>
        Nev     :   number of singular vector calculated <int>
    """
    class kdvSVLauncherError(Exception):
        pass
    
    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, traj):

        if not (isinstance(traj, Trajectory)):
            raise self.kdvSVLauncherError("traj <Trajecotory>")
        if not traj.isIntegrated:
            raise self.kdvSVLauncherError("traj not integrated")

        self.traj=traj
        
        self.grid=param.grid
        self.param=param


    #------------------------------------------------------

    def lanczos(self, Nev, tInt=None):
        """
        Call to the Lanczos procedure to calculate singular vectors 

            kdvSVLauncher.lanczos(Nev, tInt=None)

            Nev     :   number of dominant singular vectors <int>
            tInt    :   integration time <float>
            
            if tInt==None, the full reference trajectory integration
            time is taken.
        """
        if tInt==None:
            self.tInt=self.traj.tInt
        elif isinstance(tInt, (int, float)):
            if tInt<=self.traj.tInt:
                self.tInt=tInt
            else:
                raise self.kdvSVLauncherError("tInt > traj.tInt")
        else:
            raise self.kdvSVLauncherError("tInt <None|int|float>")

        self.nDt=int(self.tInt/self.traj.dt)

        self.Nev=Nev
        grid=self.grid
        param=self.param
        traj=self.traj


        sVal, sVec=fKdV.fKdVLanczos(grid.N, grid.Ntrc, grid.L,
                                    traj.dt, self.nDt, param.nDt,
                                    traj.getData(), self.Nev,
                                    param[1].getData(), 
                                    param[2].getData(),
                                    param[3].getData(), 
                                    param[4].getData())

        self.sVal=sVal
        self.sVec=sVec
        return self.sVal

    #-------------------------------------------------------

    #def __getitem__(self, idx):
    #    if idx[0]>=self.Nev:
    #        raise IndexError()
#====================================================================
#--------------------------------------------------------------------
#====================================================================

if __name__=='__main__':
    
    import matplotlib.pyplot as plt
    from kdvMisc import gauss, soliton
    
    grid=PeriodicGrid(150,300.)
    tInt=2.
    maxA=2.
    Nev=2
    
    def gaussNeg(x,t):
        x0=0.
        sig=5.
        return -0.1*gauss(x, x0, sig) 
    def sinus(x,t):
        return 0.1*np.sin(2.*2*np.pi*x/150.)

    param=Param(grid, beta=1., gamma=-1., rho=gaussNeg, forcing=sinus)
    ic=soliton(grid.x, 1., beta=1., gamma=-1. )
    
    # NL model integration
    launcher=kdvLauncher(param, maxA=maxA)
    traj=launcher.integrate(ic, tInt)
    
    svLauncher=kdvSVLauncher(param, traj)
    sVal=svLauncher.lanczos(Nev, tInt=2.)
    
    print(sVal)
    for i in xrange(Nev):
        plt.plot(grid.x,svLauncher.sVec[i])
    plt.show()
