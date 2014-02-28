import numpy as np

from pseudoSpec1D import *
from kdvParam import *

import fKdV


class kdvTLMLauncher(TLMLauncher):
    """
    TLMLauncher subclass for Augmented Korteweg-de Vries system

    kdvTLMLauncher(param)
        
        param   :   system local parametrisation <kdvParam>
    
    <!> Before integration (direct or adjoint), the TLMLauncher
        must be initialized with a reference trajectory.

    """
    class kdvTLMLauncherError(Exception):
        pass

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, traj=None):
        if not(isinstance(param, Param)):
            raise self.kdvTLMLauncherError("param <Param>")

        self.grid=param.grid
        self.param=param
        self.isTimeDependant=self.param.isTimeDependant

        self.isReferenced=False 
        if not traj==None: self.initialize(traj)


        # Status Attributes
        self.isIntegrated=False
        self.fullPertTraj=False
        self.tReal=0.

        self.propagator=self.__kdvTLMProp_Fortran
        self.propagatorAdj=self.__kdvTLMProp_Fortran_Adj
        self.propagatorSing=self.__kdvTLMSingularOp_Fortran
        self.propagatorGradTest=self.__kdvTLMGradTest_Fortran
    
    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------


    def singularOp(self, pert, tInt=None, t0=0., filtNtrc=True):
        """
        Call to the KdV singular operator 

            kdvTLMLauncher.singularOp(pert, tInt, filtNtrc=True)

            pert    :   initial perturbation <numpy.ndarray>
            tInt    :   integration time <float>
            t0      :   initial time <float>
        """
        
        if not isinstance(pert, np.ndarray):
            raise self.kdvTLMLauncherError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise self.kdvTLMLauncherError(
                                "pert.shape = (launcher.grid.N,)")
        if filtNtrc:
            specFilt(pert, self.grid)
        super(kdvTLMLauncher, self)._timeValidation(tInt, t0)

        return self.propagatorSing(pert)

    #----| Diagnostics |------------------------------------
    #-------------------------------------------------------

    def gradTest(self, ic, tInt=None, t0=0., maxPow=-10):
        """
        Gradient test

            Check consisten__kdvTLMGradTest_Fortran
ucy between the TLM adjoint and the
            nonlinear model.

            J(x)    =0.5|M(x)|^2
            gradJ(x)=L*M(x)


            kdvTLMLauncher.gradTest(ic, tInt=None, t0=0., maxPow=-10)

            ic      :   initial condition of the model <numpy.ndarray>
            tInt    :   integration time <float>
            t0      :   initial time <float>
            maxpower:   smallest power of 10 to test the difference
                        in the cost function


        """
        if not self.isReferenced:
            raise self.kdvTLMLauncherError(
                        "Not initialized with a reference trajectory")
        if not isinstance(ic, np.ndarray):
            raise self.kdvTLMLauncherError("ic <numpy.ndarray>")
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise self.kdvTLMLauncherError("ic.shape = (launcher.grid.N,)")
        super(kdvTLMLauncher, self)._timeValidation(tInt, t0)

        self.propagatorGradTest(ic, maxPow)


    #------------------------------------------------------
    #----| Private methods          |----------------------
    #----| Propagator (and adjoint) |-----------------------
    #------------------------------------------------------


    def __kdvTLMProp_Fortran(self, pert, t0=0.):

        # Local variables names
        grid=self.grid
        param=self.__t0AdjustParam(self.param, t0=t0)
        
        if self.fullPertTraj==True:
            self.pertTraj.putData(
                    fKdV.fKdVTLMPropagator(
                        grid.N, grid.Ntrc, grid.L, self.dt,
                        self.nDt, param.nDt, pert, 
                        self.refTraj.getData()[self.nDt0:
                                                self.nDt0+self.nDt+1],
                        param[1].getData(), param[2].getData(),
                        param[3].getData(),param[4].getData(), 
                        fullTraj=True))
            self.pertTraj.incrmTReal(finished=True, 
                                    tReal=self.nDt*self.dt, t0=t0)
            fPert=self.pertTraj.final

        else:
            fPert=fKdV.fKdVTLMPropagator(
                    grid.N, grid.Ntrc, grid.L, self.dt, self.nDt, param.nDt,
                    pert, 
                    self.refTraj.getData()[self.nDt0:self.nDt0+self.nDt+1],
                    param[1].getData(), param[2].getData(),
                    param[3].getData(),param[4].getData(), fullTraj=False)

        tReal=self.nDt*self.dt

        self.incrmTReal(finished=True, tReal=tReal)
        return fPert


    #------------------------------------------------------

    def __kdvTLMProp_Fortran_Adj(self, pert, t0=0.):

        # Local variables names
        grid=self.grid
        param=self.__t0AdjustParam(self.param, t0=t0)

        aPert=fKdV.fKdVTLMPropagatorAdj(
                grid.N, grid.Ntrc, grid.L, self.dt, self.nDt, param.nDt,
                pert,
                self.refTraj.getData()[self.nDt0:self.nDt0+self.nDt+1],
                param[1].getData(), param[2].getData(),
                param[3].getData(), param[4].getData())

        tReal=self.nDt*self.dt

        self.incrmTReal(finished=True, tReal=tReal)
        return aPert

    #------------------------------------------------------

    def __kdvTLMSingularOp_Fortran(self, pert, t0=0.):

        # Local variables names
        grid=self.grid
        param=self.__t0AdjustParam(self.param, t0=t0)

        LAdjLx=fKdV.fKdVTLMSingularOp(
                grid.N, grid.Ntrc, grid.L, self.dt, self.nDt, param.nDt,
                pert,
                self.refTraj.getData()[self.nDt0:self.nDt0+selfnDt+1],
                param[1].getData(), param[2].getData(), param[3].getData(),
                param[4].getData())

        tReal=2.*self.nDt*self.dt

        self.incrmTReal(finished=True, tReal=tReal)
        return LAdjLx

        
        
    #------------------------------------------------------

    def __kdvTLMGradTest_Fortran(self, ic, maxPow, t0=0.):
        
        grid=self.grid
        param=self.__t0AdjustParam(self.param, t0=t0)
        
        fKdV.fKdVTestGradient(grid.N, grid.Ntrc, grid.L, 
                self.dt, self.nDt, param.nDt, maxPow, ic,
                param[1].getData(), param[2].getData(), 
                param[3].getData(), param[4].getData(),
                param[0].getData())

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
    from kdvLauncher import *
    from kdvMisc import *

    testAdjoint=True
    tlmVsModel=False

    grid=PeriodicGrid(150,300.)
    tInt=3.
    maxA=2.

    param=Param(grid, beta=1., gamma=-1.)
    
    #----| Reference trajectory |-----------------
    u0=rndSpecVec(grid, Ntrc=grid.Ntrc/5,  amp=0.3, seed=0.1)
    M=kdvLauncher(param, maxA)
    u=M.integrate(u0, tInt)

    
    #----| Gradient test |------------------------
    L=kdvTLMLauncher(param, traj=u)
    L.gradTest(u0)

    #----| TLM vs NL model |----------------------
    if tlmVsModel:
        du=soliton(grid.x, 0. , amp=2., beta=1., gamma=-1. )
        L=kdvTLMLauncher(param, traj=u)
    
        u_pert=M.integrate(u0+du)
        pert=L.integrate(du, fullPertTraj=True)
        plt.plot(grid.x, u_pert.final)
        plt.plot(grid.x, u.final+pert)


    #----| Adjoint testing |----------------------
    if testAdjoint:
        Ntrc=grid.Ntrc
        print("Testting adjoint validity")
        L=kdvTLMLauncher(param, traj=u)
    
        dx=rndSpecVec(grid, Ntrc=Ntrc, amp=0.2, seed=0.2)
        dy=rndSpecVec(grid, Ntrc=Ntrc, amp=0.2, seed=0.3)
    
        Ldy=L.integrate(dy)
        print("dy     >>|%3d(%.3f) - L - %3d(%.3f)|>> Ldy"%(
                        L.nDt0,  L.t0, L.nDtFinal, L.tFinal ))
        Adx=L.adjoint(dx)
        print("Adx    <<|%3d(%.3f) - L*- %3d(%.3f)|>>  dx"%(
                        L.nDt0,  L.t0, L.nDtFinal, L.tFinal ))
        print("<dx, Ldy> = %+.15g"%np.dot(dx, Ldy))
        print("<Adx, dt> = %+.15g"%np.dot(Adx, dy))
        print("<dx, Ldy>-<Adx, dt> = %+.15g"%(
                    np.dot(dx, Ldy)-np.dot(Adx, dy)))


    #----| partial trajectory time |----
    print("\nTestting adjoint validity for partial integration")
    L=kdvTLMLauncher(param, traj=u)
    Ldy=L.integrate(dy, tInt=tInt/3., t0=tInt/2.)
    print("dy     >>|%3d(%.3f) - L - %3d(%.3f)|>> Ldy"%(
                    L.nDt0,  L.t0, L.nDtFinal, L.tFinal ))
    LAdj_x=L.adjoint(dx, tInt=tInt/3., t0=tInt/2.)
    print("L*dx   <<|%3d(%.3f) - L*- %3d(%.3f)|>>  dx"%(
                    L.nDt0,  L.t0, L.nDtFinal, L.tFinal ))
    print("<dx, Ldy>-<L*dx, dt> = %+.15g"%(
                np.dot(dx, Ldy)-np.dot(LAdj_x, dy)))


    #----| step integrations |----------
    print("\nTestting adjoint validity for successive step integrations")
    L1=kdvTLMLauncher(param, traj=u)
    L2=kdvTLMLauncher(param, traj=u)

    L1dy=L1.integrate(dy, tInt=tInt/3., t0=0.)
    print("dy     >>|%3d(%.3f) - L1 - %3d(%.3f)|>> L1dy"%(
                    L1.nDt0,  L1.t0, L1.nDtFinal, L1.tFinal ))
    L2L1dy=L2.integrate(L1dy, tInt=tInt/2., t0=L1.tFinal)
    print("L1dy   >>|%3d(%.3f) - L2 - %3d(%.3f)|>> L2L1dy"%(
                    L2.nDt0,  L2.t0, L2.nDtFinal, L2.tFinal ))

    A1dx=L1.adjoint(dx, tInt=tInt/3., t0=0.)
    print("A1dx   <<|%3d(%.3f) - L1*- %3d(%.3f)|<< dx"%(
                    L1.nDt0,  L1.t0, L1.nDtFinal, L1.tFinal ))
    A2dx=L2.adjoint(dx, tInt=tInt/2., t0=L1.tFinal)
    print("A2dx   <<|%3d(%.3f) - L2*- %3d(%.3f)|<< dx"%(
                    L2.nDt0,  L2.t0, L2.nDtFinal, L2.tFinal ))
    A1A2dx=L1.adjoint(A2dx, tInt=tInt/3., t0=0.)
    print("A1A2dx <<|%3d(%.3f) - L1*- %3d(%.3f)|<< A2dx"%(
                    L1.nDt0,  L1.t0, L1.nDtFinal, L1.tFinal ))


    print("<dx,L1dy>-<A1dx, dy>     = %+.15g"%(np.dot(dx, L1dy)\
                                            -np.dot(A1dx, dy)))
    print("<dx,L2L1dy>-<A2dx, L1dy> = %+.15g"%(np.dot(dx, L2L1dy)\
                                            -np.dot(A2dx, L1dy)))
    print("<dx,L2L1dy>-<A1A2dx, dy> = %+.15g"%(np.dot(dx, L2L1dy)\
                                            -np.dot(A1A2dx, dy)))
