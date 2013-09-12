import numpy as np

from pseudoSpec1D import *
from kdvParam import *

import fKdV


class TLMLauncher(object):
    """
    TLMLauncher class
    """
    doublePrecisionTolerance=1e-13
    class TLMLauncherError(Exception):
        pass

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, traj, tInt=None, t0=0.):
        if not(isinstance(param, Param)):
            raise self.TLMLauncherError("param <Param>")

        self.param=param
        self.grid=param.grid
        
        if not(isinstance(traj, Trajectory)):
            raise self.TLMLauncherError("traj <Trajectory>")
        if not (traj.grid==self.grid):
            raise SpectralGridError("traj.grid <> grid")
        self.refTraj=traj

        # Time attributes
        self.dt=self.refTraj.dt
        if tInt==None : tInt=self.refTraj.tInt
        if t0<0:
            raise self.TLMLauncherError("t0>=0.")
        if t0+tInt>self.refTraj.tInt:
            raise self.TLMLauncherError("t0+tInt<=self.refTraj.tInt")
        self.nDt=int((tInt)/self.dt)
        self.nDt0=int(t0/self.dt)
        self.nDtFinal=self.nDt0+self.nDt
        # real times from step integers
        self.tInt=self.nDt*self.dt
        self.t0=self.nDt0*self.dt
        self.tFinal=self.nDtFinal*self.dt

        # Status Attributes
        self.isIntegrated=False
        self.fullPertTraj=False
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

    def integrate(self, pert, fullPertTraj=False):
        
        if not isinstance(pert, np.ndarray):
            raise self.TLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise self.TLMLauncherError("pert.shape = (launcher.grid.N,)")

        if fullPertTraj:
            self.fullPertTraj=True
            self.pertTraj=Trajectory(self.grid)
            self.pertTraj.initialize(pert, tInt, dt)
        else:
            self.fullPertTraj=False

        return self.__kdvTLMProp_Fortran(pert)
                                            

    #-------------------------------------------------------

    def adjoint(self, pert):         
        
        if not isinstance(pert, np.ndarray):
            raise self.TLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise self.TLMLauncherError("pert.shape = (launcher.grid.N,)")

        return self.__kdvTLMProp_Fortran_Adj(pert) 
    
    #-------------------------------------------------------

    def singularOp(self, pert):
        
        if not isinstance(pert, np.ndarray):
            raise self.TLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise self.TLMLauncherError("pert.shape = (launcher.grid.N,)")
        if tInt==None : tInt=self.refTraj.tInt

        return self.__kdvTLMSingularOp_Fortran(pert)

    #----| Diagnostics |------------------------------------
    #-------------------------------------------------------

    def gradTest(self, pert, maxPow=-10):
        if not isinstance(pert, np.ndarray):
            raise self.TLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise self.TLMLauncherError("pert.shape = (launcher.grid.N,)")
        
        grid=self.grid
        param=self.param
        fKdV.fKdVTestGradient(grid.N, grid.Ntrc, grid.L, 
                self.dt, self.nDt, maxPow, pert,
                self.refTraj.getData()[self.nDt0:self.nDt0+self.nDt+1],
                param[1], param[2], param[3], param[4])

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------

    #----| Propagator (and adjoint) |-----------------------
    #------------------------------------------------------

    def __kdvTLMProp_Fortran(self, pert, fullOutput=False):

        # Local variables names
        grid=self.grid
        param=self.param

        if self.fullPertTraj:
            self.pertTraj.putData(fKdV.fKdVTLMPropagator(
                grid.N, grid.Ntrc, grid.L, self.dt, self.nDt,pert, 
                self.refTraj.getData()[self.nDt0:self.nDt0+self.nDt+1],
                param[1], param[2], param[3], param[4],
                fullTraj=True))

            tReal=self.nDt*self.dt
            self.pertTraj.incrmTReal(finished=True, tReal=tReal)
            fPert=self.pertTraj[nDt]
        else:
            fPert=fKdV.fKdVTLMPropagator(
                grid.N, grid.Ntrc, grid.L, self.dt, self.nDt, pert, 
                self.refTraj.getData()[self.nDt0:self.nDt0+self.nDt+1],
                param[1], param[2], param[3], param[4], 
                fullTraj=False)

            tReal=self.nDt*self.dt

        self.incrmTReal(finished=True, tReal=tReal)
        return fPert


    #------------------------------------------------------

    def __kdvTLMProp_Fortran_Adj(self, pert):

        # Local variables names
        grid=self.grid
        param=self.param

        aPert=fKdV.fKdVTLMPropagator(
                grid.N, grid.Ntrc, grid.L, self.dt, self.nDt, pert,
                self.refTraj.getData()[self.nDt0:self.nDt0+self.nDt+1],
                param[1], param[2], param[3], param[4], 
                fullTraj=False)

        tReal=self.nDt*self.dt

        self.incrmTReal(finished=True, tReal=tReal)
        return aPert

    #------------------------------------------------------

    def __kdvTLMSingularOp_Fortran(self, pert):

        # Local variables names
        grid=self.grid
        param=self.param

        LAdjLx=fKdV.fKdVTLMSingularOp(
                grid.N, grid.Ntrc, grid.L, self.dt, self.nDt, pert,
                self.refTraj.getData()[self.nDt0:self.nDt0+selfnDt+1],
                param[1], param[2], param[3], param[4])

        tReal=2.*self.nDt*self.dt

        self.incrmTReal(finished=True, tReal=tReal)
        return LAdjLx



#--------------------------------------------------------------------
#====================================================================
#--------------------------------------------------------------------



if __name__=='__main__':
    import matplotlib.pyplot as plt
    from kdvLauncher import *
    from kdvMisc import *

    grid=SpectralGrid(150,300.)
    tInt=3.
    maxA=2.


    param=Param(grid, beta=1., gamma=-1.)


    
    # Adjoint testing
    print("Testting adjoint validity")
    dx=rndFiltVec(grid, amp=0.2)
    dy=rndFiltVec(grid, amp=0.2)
    u0=rndFiltVec(grid, Ntrc=grid.Ntrc/2,  amp=1.)

    M=Launcher(param, u0)
    u=M.integrate(tInt, 3.)
    L=TLMLauncher(param, u)

    Ldy=L.integrate(dy)
    LAdj_x=L.adjoint(dx)
    print("<dx, Ldy>-<L*dx, dt> = %+.15g"%(
                np.dot(dx, Ldy)-np.dot(LAdj_x, dy)))

    #----| partial trajectory time |----
    print("\nTestting adjoint validity for partial integration")
    L=TLMLauncher(param, u, tInt=tInt/3., t0=tInt/2.)
    Ldy=L.integrate(dy)
    print("%d(%.3f) --| L |---> %d(%.3f)"%(L.nDt0,  L.t0,
                                           L.nDtFinal, L.tFinal ))
    LAdj_x=L.adjoint(dx)
    print("<dx, Ldy>-<L*dx, dt> = %+.15g"%(
                np.dot(dx, Ldy)-np.dot(LAdj_x, dy)))


    #----| step integrations |----------
    print("\nTestting adjoint validity for successive step integrations")
    L1=TLMLauncher(param, u, tInt=tInt/3., t0=0.)
    L2=TLMLauncher(param, u, tInt=tInt/2., t0=L1.tFinal)

    L1dy=L1.integrate(dy)
    print("%d(%.3f) --| L1 |---> %d(%.3f)"%(L1.nDt0,  L1.t0,
                                            L1.nDtFinal, L1.tFinal ))
    L2L1dy=L2.integrate(L1dy)
    print("%d(%.3f) --| L2 |---> %d(%.3f)"%(L2.nDt0,  L2.t0,
                                            L2.nDtFinal, L2.tFinal ))

    A1dx=L1.adjoint(dx)
    A2dx=L2.adjoint(dx)
    A1A2dx=L1.adjoint(A2dx)


    print("<dx,L1dy>-<A1dx, dy>     = %+.15g"%(np.dot(dx, L1dy)\
                                            -np.dot(A1dx, dy)))
    print("<dx,L2L1dy>-<A2dx, L1dy> = %+.15g"%(np.dot(dx, L2L1dy)\
                                            -np.dot(A2dx, L1dy)))
    print("<dx,L2L1dy>-<A1A2dx, dy> = %+.15g"%(np.dot(dx, L2L1dy)\
                                            -np.dot(A1A2dx, dy)))
