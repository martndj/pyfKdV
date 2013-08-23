import numpy as np

from pseudoSpec1D import *
from kdvParam import *
from kdvLauncher import *
from kdvTLMLauncher import *

import fKdV

class SVLauncher(Launcher):
    """
    """
    
    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, traj, Nev, tInt=None):

        if not (isinstance(traj, Trajectory)):
            raise self.LauncherError("traj <Trajecotory>")
        if not traj.isIntegrated:
            raise self.LauncherError("traj not integrated")

        self.traj=traj
        self.Nev=Nev
        
        if tInt==None:
            self.tInt=self.traj.tInt
        elif isinstance(tInt, (int, float)):
            if tInt<=self.traj.tInt:
                self.tInt=tInt
            else:
                raise self.LauncherError("tInt > traj.tInt")
        else:
            raise self.LauncherError("tInt <None|int|float>")

        super(SVLauncher, self).__init__(param, self.traj.ic)
        #Launcher.__init__(self, param, self.traj.ic)


    #------------------------------------------------------

    def lanczos(self):

        grid=self.grid
        param=self.param
        traj=self.traj

        self.nDt=int(self.tInt/traj.dt)

        sVal, sVec=fKdV.fKdVLanczos(grid.N, grid.Ntrc, grid.L,
                                    traj.dt, self.nDt, traj.getData(), 
                                    self.Nev,
                                    param[1], param[2], param[3], param[4])

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
    
    grid=SpectralGrid(150,300.)
    tInt=2.
    maxA=2.
    
    def gaussNeg(x):
        x0=0.
        sig=5.
        return -0.1*gauss(x, x0, sig) 
    def sinus(x):
        return 0.1*sin(2.*2*np.pi*x/150.)

    param=Param(grid, beta=1., gamma=-1., rho=gaussNeg, forcing=sinus)
    ic=soliton(grid.x, 1., beta=1., gamma=-1. )
    
    # NL model integration
    launcher=Launcher(param, ic)
    
    traj=launcher.integrate(tInt, maxA)
    
    svLauncher=SVLauncher(param, traj,2, tInt=2.3)
    sVal=svLauncher.lanczos()
    
    print(sVal)
