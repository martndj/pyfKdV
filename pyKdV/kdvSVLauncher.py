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

    def lanczos(self, nSV, tOpt=None, straighten=True):
        """
        Call to the Lanczos procedure to calculate singular vectors 

            kdvSVLauncher.lanczos(nSV, tOpt=None)

            nSV     :   number of dominant singular vectors <int>
            tOpt    :   optimization (integration) time <float>
            
            if tOpt==None, the full reference trajectory integration
            time is taken.
        """
        if tOpt==None:
            self.tOpt=self.refTraj.tReal
        elif isinstance(tOpt, (int, float)):
            if tOpt<=self.refTraj.tReal:
                self.tOpt=tOpt
            else:
                raise self.kdvSVLauncherError("tOpt > traj.tReal")
        else:
            raise self.kdvSVLauncherError("tOpt <None|int|float>")

        self.nDt=int(self.tOpt/self.refTraj.dt)

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
        if straighten: self.straightenSV()
        return self.sVal

    #------------------------------------------------------

    def straightenSV(self):

        if not self.isCalculated:
            raise self.kdvSVLauncherError("SV must be calculated first!")

        for i in xrange(self.nSV):
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

#if __name__=='__main__':
