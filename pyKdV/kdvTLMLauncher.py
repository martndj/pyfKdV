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

    def __init__(self, param, traj, pert):
        if not(isinstance(param, Param)):
            raise self.TLMLauncherError("param <Param>")

        self.param=param
        self.grid=param.grid
        
        if not(isinstance(traj, Trajectory)):
            raise self.TLMLauncherError("traj <Trajectory>")
        if not (traj.grid==self.grid):
            raise SpectralGridError("traj.grid <> grid")
        self.refTraj=traj

        if not isinstance(pert, np.ndarray):
            raise self.TLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise self.TLMLauncherError("pert.shape = (launcher.grid.N,)")

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

    def integrate(self, tInt=None, t0=0., fullPertTraj=False):
        
        if tInt==None : tInt=self.refTraj.tInt

        dt=self.refTraj.dt
        if t0<0:
            raise self.TLMLauncherError("t0>=0.")
        if t0+tInt>self.refTraj.tInt:
            raise self.TLMLauncherError("t0+tInt<=self.refTraj.tInt")
        nDt=int((tInt)/dt)
        nDt0=int(t0/dt)

        return  self.integrateNDt(nDt, nDt0, fullPertTraj=fullPertTraj)
   
    #-----------------------------------

    def integrateNDt(self, nDt, nDt0, fullPertTraj=False):
        self.nDt=nDt
        self.nDt0=nDt0

        if fullPertTraj:
            self.fullPertTraj=True
            self.pertTraj=Trajectory(self.grid)
            self.pertTraj.initialize(self.pert, tInt, dt)
        else:
            self.fullPertTraj=False

        return self.__kdvTLMProp_Fortran(self.refTraj.dt, nDt, nDt0, 
                                            fullOutput=self.fullPertTraj)

    #-------------------------------------------------------

    def adjoint(self, tInt=None, t0=0.):         

        if tInt==None : tInt=self.refTraj.tInt

        dt=self.refTraj.dt
        if t0<0:
            raise self.TLMLauncherError("t0>=0.")
        if t0+tInt>self.refTraj.tInt:
            raise self.TLMLauncherError("t0+tInt<=self.refTraj.tInt")
        nDt=int((tInt)/dt)
        nDt0=int(t0/dt)
        
        return self.adjointNDt(nDt, nDt0)

    #-----------------------------------
    def adjointNDt(self, nDt, nDt0):
        self.nDt=nDt
        self.nDt0=nDt0
        return self.__kdvTLMProp_Fortran_Adj(self.refTraj.dt, nDt, nDt0) 
    
    #-------------------------------------------------------

    def singularOp(self, tInt=None, t0=0.):
        
        if tInt==None : tInt=self.refTraj.tInt

        dt=self.refTraj.dt
        if t0<0:
            raise self.TLMLauncherError("t0>=0.")
        if t0+tInt>self.refTraj.tInt:
            raise self.TLMLauncherError("t0+tInt<=self.refTraj.tInt")
        nDt=int((tInt)/dt)
        nDt0=int(t0/dt)

        return self.__kdvTLMSingularOp_Fortran(dt, nDt, nDt0)

    #----| Diagnostics |------------------------------------
    #-------------------------------------------------------

    def gradTest(self, tInt=None, t0=0.):
        
        if tInt==None : tInt=self.refTraj.tInt

        dt=self.refTraj.dt
        if t0<0:
            raise self.TLMLauncherError("t0>=0.")
        if t0+tInt>self.refTraj.tInt:
            raise self.TLMLauncherError("t0+tInt<=self.refTraj.tInt")
        nDt=int((tInt)/dt)
        nDt0=int(t0/dt)

        grid=self.grid
        param=self.param
        fKdV.fKdVTestGradient(grid.N, grid.Ntrc, grid.L, 
                                dt, nDt, -10, self.pert,
                                self.refTraj.getData()[nDt0:nDt0+nDt+1],
                                param[1], param[2], param[3], param[4])

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------

    #----| Propagator (and adjoint) |-----------------------
    #------------------------------------------------------

    def __kdvTLMProp_Fortran(
            self, dt, nDt, nDt0=0, fullOutput=True):

        # Local variables names
        grid=self.grid
        param=self.param

        if fullOutput:
            self.pertTraj.putData(fKdV.fKdVTLMPropagator(
                                grid.N, grid.Ntrc, grid.L, dt, nDt,
                                self.pert, 
                                self.refTraj.getData()[nDt0:nDt0+nDt+1],
                                param[1], param[2], param[3], param[4], 
                                fullTraj=True))

            tReal=nDt*dt
            self.pertTraj.incrmTReal(finished=True, tReal=tReal)
            fPert=self.pertTraj[nDt]
        else:
            fPert=fKdV.fKdVTLMPropagator(
                                grid.N, grid.Ntrc, grid.L, dt, nDt,
                                self.pert, 
                                self.refTraj.getData()[nDt0:nDt0+nDt+1],
                                param[1], param[2], param[3], param[4], 
                                fullTraj=False)
            tReal=nDt*dt

        self.incrmTReal(finished=True, tReal=tReal)
        return fPert


    #------------------------------------------------------

    def __kdvTLMProp_Fortran_Adj(self, dt, nDt, nDt0=0):

        # Local variables names
        grid=self.grid
        param=self.param

        aPert=fKdV.fKdVTLMPropagator(
                            grid.N, grid.Ntrc, grid.L, dt, nDt,
                            self.pert,
                            self.refTraj.getData()[nDt0:nDt0+nDt+1],
                            param[1], param[2], param[3], param[4], 
                            fullTraj=False)
        tReal=nDt*dt

        self.incrmTReal(finished=True, tReal=tReal)

        return aPert

    #------------------------------------------------------

    def __kdvTLMSingularOp_Fortran(self, dt, nDt, nDt0=0):

        # Local variables names
        grid=self.grid
        param=self.param

        LAdjLx=fKdV.fKdVTLMSingularOp(
                                grid.N, grid.Ntrc, grid.L, dt, nDt,
                                self.pert,
                                self.refTraj.getData()[nDt0:nDt0+nDt+1],
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
    from kdvMisc import *

    grid=SpectralGrid(150,300.)
    tInt=3.
    maxA=2.


    param=Param(grid, beta=1., gamma=-1.)


    
    # Adjoint testing
    print("Testting adjoitn validity <dx, Ldy>==<L*dx, dy>")
    dx=rndFiltVec(grid, amp=0.2)
    dy=rndFiltVec(grid, amp=0.2)
    u0=rndFiltVec(grid, Ntrc=grid.Ntrc/2,  amp=1.)
    launcher=Launcher(param, u0)
    u=launcher.integrate(tInt, 3.)
    dyLauncher=TLMLauncher(param, u, dy)
    dxLauncher=TLMLauncher(param, u, dx)
    Ldy=dyLauncher.integrate()
    LAdj_x=dxLauncher.adjoint()
    print(np.dot(dx, Ldy)-np.dot(LAdj_x, dy))

    #----| partial trajectory time |----
    Ldy=dyLauncher.integrate(tInt/3., tInt/2.)
    LAdj_x=dxLauncher.adjoint(tInt/3., tInt/2.)
    print(np.dot(dx, Ldy)-np.dot(LAdj_x, dy))


    #----| step integrations |----------
    step=tInt/3.
    L1=TLMLauncher(param, u, dy)
    L1dy=L1.integrate(step, t0=0.)

    L2=TLMLauncher(param, u, L1dy)
    L2L1dy=L2.integrate(step, t0=L1.tReal)
    
    A2=TLMLauncher(param, u, L2L1dy)
    A2dx=A2.adjoint(step, t0=L1.tReal)

    A1=TLMLauncher(param, u, A2dx)
    A1A2dx=A1.adjoint(step, t0=0.)
    print(np.dot(dx, L2L1dy)-np.dot(A1A2dx, dy))

    #----| step integrations |----------
    step=tInt/3.
    L1=TLMLauncher(param, u, dy)
    L1dy=L1.integrate(step, t0=0.)
    nDt_1=L1.nDt
    nDt0_1=L1.nDt0
    print(nDt0_1, nDt_1)
    
    L2=TLMLauncher(param, u, L1dy)
    L2L1dy=L2.integrate(step, t0=L1.tReal)
    nDt_2=L2.nDt
    nDt0_2=L2.nDt0
    print(nDt0_2, nDt_2)

    A2=TLMLauncher(param, u, L2L1dy)
    A2dx=A2.adjointNDt(nDt_2, nDt0_2)
    A1=TLMLauncher(param, u, A2dx)
    A1A2dx=A2.adjointNDt(nDt_1, nDt0_1)

    A1=TLMLauncher(param, u, dx)
    A1dx=dxLauncher.adjointNDt(nDt_1, nDt0_1)
    print("<dx,L1dy>-<A1dx, dy>     = %+.15g"%(np.dot(dx, L1dy)\
                                            -np.dot(A1dx, dy)))
    print("<dx,L2L1dy>-<A2dx, L1dy> = %+.15g"%(np.dot(dx, L2L1dy)\
                                            -np.dot(A2dx, L1dy)))
    print("<dx,L2L1dy>-<A1A2dx, dy> = %+.15g"%(np.dot(dx, L2L1dy)\
                                            -np.dot(A1A2dx, dy)))

