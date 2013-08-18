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

    def __init__(self, param, traj, Nev):

        if not (isinstance(traj, Trajectory)):
            raise self.LauncherError("traj <Trajecotory>")
        if not traj.isIntegrated:
            raise self.LauncherError("traj not integrated")

        self.traj=traj
        self.Nev=Nev

        super(SVLauncher, self).__init__(param, self.traj.ic)
        #Launcher.__init__(self, param, self.traj.ic)


    #------------------------------------------------------

    def lanczos(self):

        grid=self.grid
        param=self.param
        traj=self.traj

        sVal, sVec=fKdV.fKdVLanczos(grid.N, grid.Ntrc, grid.L,
                                    traj.dt, traj.nDt, traj.getData(), 
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
    
    grid=SpectralGrid(150,300.)
    tInt=2.
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
    
    svLauncher=SVLauncher(param, traj,4)
    sVal=svLauncher.lanczos()
    
    print(sVal)
