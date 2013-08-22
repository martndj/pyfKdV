import numpy as np

from pseudoSpec1D import *
from kdvParam import *
#from kdvLauncher import *

import fKdV

class TLMLauncherError(Exception):
    pass

class TLMLauncher(object):
    """
    TLMLauncher class
    """
    doublePrecisionTolerance=1e-13

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, traj, pert):
        if not(isinstance(param, Param)):
            raise TLMLauncherError("param <Param>")

        self.param=param
        self.grid=param.grid
        
        if not(isinstance(traj, Trajectory)):
            raise TLMLauncherError("traj <Trajectory>")
        if not (traj.grid==self.grid):
            raise SpectralGridError("traj.grid <> grid")
        self.refTraj=traj

        if not isinstance(pert, np.ndarray):
            raise TLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise TLMLauncherError("pert.shape = (launcher.grid.N,)")

        # implicit filtering
        #self.pert=specFilt(pert, self.grid)
        self.pert=pert

        # Status Attributes
        self.isIntegrated=False
        self.fullPertTraj=False

        # Time attribute
        self.tInt=0.
        self.tReal=0.
    
    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------

    
    def incrmTReal(self, finished=False, tReal=None):
        if tReal<>None:
            self.tReal=tReal
        if finished:
            self.isIntegrated=True
        else:
            self.tReal+=self.refTraj.dt

    #-------------------------------------------------------

    def integrate(self, tInt, fullPertTraj=False):
        

        dt=self.refTraj.dt
        nDt=int(tInt/dt)

        if fullPertTraj:
            self.fullPertTraj=True
            self.pertTraj=Trajectory(self.grid)
            self.pertTraj.initialize(self.pert, tInt, dt)

        return self.__kdvTLMProp_Fortran(dt, nDt, fullPertTraj=fullPertTraj)
    
    #-------------------------------------------------------

    def adjoint(self, tInt, resetTReal=True): #, fullPertTraj=False
        

        dt=self.refTraj.dt
        nDt=int(tInt/dt)
        
        #if fullPertTraj:
        #    self.fullPertTraj=True
        #    self.pertTraj=Trajectory(self.grid)
        #    self.pertTraj.initialize(self.pert, tInt, dt)


        return self.__kdvTLMProp_Fortran_Adj(dt, nDt) 
        #                        fullPertTraj=fullPertTraj)

    
    #-------------------------------------------------------

    def singularOp(self, tInt):
        

        dt=self.refTraj.dt
        nDt=int(tInt/dt)

        return self.__kdvTLMSingularOp_Fortran(dt, nDt)

    #-------------------------------------------------------

    def gradTest(self):
        
        dt=self.refTraj.dt
        nDt=int(tInt/dt)
        grid=self.grid
        param=self.param
        fKdV.fKdVTestGradient(grid.N, grid.Ntrc, grid.L, 
                                dt, nDt, -10, self.pert,
                                self.refTraj.getData(),
                                param[1], param[2], param[3], param[4])
                                
    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------

    #----| Propagator (and adjoint) |-----------------------
    #------------------------------------------------------

    def __kdvTLMProp_Fortran(
            self, dt, nDt, fullPertTraj=True):

        # Local variables names
        grid=self.grid
        param=self.param

        if fullPertTraj:
            self.pertTraj.putData(fKdV.fKdVTLMPropagator(
                                grid.N, grid.Ntrc, grid.L, dt, nDt,
                                self.pert, self.refTraj.getData(),
                                param[1], param[2], param[3], param[4], 
                                fullTraj=True))

            tReal=nDt*dt
            self.pertTraj.incrmTReal(finished=True, tReal=tReal)
            fPert=self.pertTraj[nDt]
        else:
            fPert=fKdV.fKdVTLMPropagator(
                                grid.N, grid.Ntrc, grid.L, dt, nDt,
                                self.pert, self.refTraj.getData(),
                                param[1], param[2], param[3], param[4], 
                                fullTraj=False)
            tReal=nDt*dt

        self.incrmTReal(finished=True, tReal=tReal)
        return fPert


    #------------------------------------------------------

    def __kdvTLMProp_Fortran_Adj(
            self, dt, nDt, fullPertTraj=True):

        # Local variables names
        grid=self.grid
        param=self.param

        if fullPertTraj:
            self.pertTraj.putData(fKdV.fKdVTLMPropagatorAdj(
                                grid.N, grid.Ntrc, grid.L, dt, nDt,
                                self.pert, self.refTraj.getData(),
                                param[1], param[2], param[3], param[4], 
                                fullTraj=True))

            tReal=nDt*dt
            self.pertTraj.incrmTReal(finished=True, tReal=tReal)
            aPert=self.pertTraj[0]
        else:
            aPert=fKdV.fKdVTLMPropagator(
                                grid.N, grid.Ntrc, grid.L, dt, nDt,
                                self.pert, self.refTraj.getData(),
                                param[1], param[2], param[3], param[4], 
                                fullTraj=False)
            tReal=nDt*dt

        self.incrmTReal(finished=True, tReal=tReal)

        return aPert

    #------------------------------------------------------

    def __kdvTLMSingularOp_Fortran(self, dt, nDt):

        # Local variables names
        grid=self.grid
        param=self.param

        LAdjLx=fKdV.fKdVTLMSingularOp(
                                grid.N, grid.Ntrc, grid.L, dt, nDt,
                                self.pert, self.refTraj.getData(),
                                param[1], param[2], param[3], param[4])
        tReal=2.*nDt*dt

        self.incrmTReal(finished=True, tReal=tReal)
        return LAdjLx



#--------------------------------------------------------------------
#====================================================================
#--------------------------------------------------------------------



if __name__=='__main__':
    import matplotlib.pyplot as plt
    from kdvLauncher import *
    from kdvMisc import soliton, gauss

    grid=SpectralGrid(150,300.)
    tInt=3.
    maxA=2.


    param=Param(grid, beta=1., gamma=-1.)

    ic=soliton(param, 0., 1.)

    # NL model integration
    launcher=Launcher(param, ic)
    
    traj=launcher.integrate(tInt, maxA)


    pert=0.1*gauss(param, -10., grid.L/25. )

    tLauncher=TLMLauncher(param, traj, pert)
    tLauncher.gradTest()
    fPert=tLauncher.integrate(tInt, fullPertTraj=True)

    aPert=tLauncher.adjoint(tInt)

    sLauncher=TLMLauncher(param, traj, pert)
    sPert=sLauncher.singularOp(tInt)

    plt.plot(grid.x, fPert)
    plt.plot(grid.x, aPert)
    plt.plot(grid.x, sPert)
        
    plt.show()
