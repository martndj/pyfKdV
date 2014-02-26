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

        kdvSVLauncher(traj, param).lanczos(nSV)

        traj    :   reference trajectory <Trajectory>
        param   :   parametrisation <kdvParam>
        nSV     :   number of singular vector calculated <int>
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

        self.refTraj=traj
        
        self.grid=param.grid
        self.param=param
        self.isCalculated=False


    #------------------------------------------------------

    def lanczos(self, nSV, tInt=None):
        """
        Call to the Lanczos procedure to calculate singular vectors 

            kdvSVLauncher.lanczos(nSV, tInt=None)

            nSV     :   number of dominant singular vectors <int>
            tInt    :   integration time <float>
            
            if tInt==None, the full reference trajectory integration
            time is taken.
        """
        if tInt==None:
            self.tInt=self.refTraj.tInt
        elif isinstance(tInt, (int, float)):
            if tInt<=self.refTraj.tInt:
                self.tInt=tInt
            else:
                raise self.kdvSVLauncherError("tInt > traj.tInt")
        else:
            raise self.kdvSVLauncherError("tInt <None|int|float>")

        self.nDt=int(self.tInt/self.refTraj.dt)

        self.nSV=nSV
        self.Nev=self.nSV # retro compatibility
        grid=self.grid
        param=self.param
        traj=self.refTraj


        sVal, sVec=fKdV.fKdVLanczos(grid.N, grid.Ntrc, grid.L,
                                    traj.dt, self.nDt, param.nDt,
                                    traj.getData(), self.nSV,
                                    param[1].getData(), 
                                    param[2].getData(),
                                    param[3].getData(), 
                                    param[4].getData())

        self.sVal=sVal
        self.sVec=sVec
        self.isCalculated=True
        return self.sVal

    #------------------------------------------------------

    def straightenSV(self):

        if not self.isCalculted:
            raise self.kdvSVLauncherError("SV must be calculated first!")

        for i in self.nSV:
            if np.abs(np.min(self.sVec[i]))>np.max(self.sVec[i]):
                self.sVec[i]=-self.sVec[i]
            else:
                pass

    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="####| svLauncher |#####################################\n"
        output+="\n  parametrisation:\n"
        output+=self.param.__str__()
        output+="\n  reference trajectory:\n"
        output+=self.refTraj.__str__()
        if self.isCalculated:
            output+="\n %d SV obtained:\n"
            output+=str(self.sVal)
        output+=(
            "\n#######################################################\n")
        return output

    #-------------------------------------------------------

    #def __getitem__(self, idx):
    #    if idx[0]>=self.nSV:
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
    nSV=2
    
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
    sVal=svLauncher.lanczos(nSV, tInt=2.)
    
    print(sVal)
    for i in xrange(nSV):
        plt.plot(grid.x,svLauncher.sVec[i])
    plt.show()
