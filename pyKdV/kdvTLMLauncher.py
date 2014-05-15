import numpy as np

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

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------

    def __init__(self, param, traj=None):
        if not(isinstance(param, Param)):
            raise TypeError("param <Param>")

        self.param=param
        self.isTimeDependant=self.param.isTimeDependant
        
        super(kdvTLMLauncher, self).__init__(param.grid, traj=traj)

        self.propagator=self.__kdvTLMProp_Fortran
        self.propagatorAdj=self.__kdvTLMProp_Fortran_Adj
        self.propagatorSing=self.__kdvTLMSingularOp_Fortran
        self.propagatorGradTest=self.__kdvTLMGradTest_Fortran
    
    #------------------------------------------------------
    #----| Public methods |--------------------------------
    #------------------------------------------------------


    def singularOp(self, pert, tInt=None, t0=0.):
        """
        Call to the KdV singular operator 

            kdvTLMLauncher.singularOp(pert, tInt, filtNtrc=True)

            pert    :   initial perturbation <numpy.ndarray>
            tInt    :   integration time <float>
            t0      :   initial time <float>
        """
        
        if not isinstance(pert, np.ndarray):
            raise TypeError("pert <numpy.ndarray>")
        if pert.ndim <> 1 or pert.size <> self.grid.N:
            raise TypeError(
                                "pert.shape = (launcher.grid.N,)")
        super(kdvTLMLauncher, self)._timeValidation(tInt, t0)

        return self.propagatorSing(pert)

    #----| Diagnostics |------------------------------------
    #-------------------------------------------------------

    def gradTestFortran(self, ic, dt, nDt, t0=0., maxPow=-10):
        """
        Gradient test

            Check consistency between the TLM adjoint and the
            nonlinear model.

            J(x)    =0.5|M(x)|^2
            gradJ(x)=L*M(x)


            kdvTLMLauncher.gradTest(ic, dt, nDt, t0=0., maxPow=-10)

            ic      :   initial condition of the model <numpy.ndarray>
            dt      :   integration step <float>
            nDt     :   number of integration steps <int>
            t0      :   initial time <float>
            maxpower:   smallest power of 10 to test the difference
                        in the cost function


        """
        if not isinstance(ic, np.ndarray):
            raise TypeError("ic <numpy.ndarray>")
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise TypeError("ic.shape = (launcher.grid.N,)")
        self.propagatorGradTest(ic, dt, nDt, maxPow, t0=t0)


    #------------------------------------------------------
    #----| Private methods          |----------------------
    #----| Propagator (and adjoint) |-----------------------
    #------------------------------------------------------


    def __kdvTLMProp_Fortran(self, pert, traj, t0=0.):

        # Local variables names
        grid=self.grid
        param=self.__t0AdjustParam(self.param, t0=t0)
        
        trajData=fKdV.fKdVTLMPropagator(
                    grid.N, grid.Ntrc, grid.L, self.dt,
                    self._nDt, param.nDt, pert, 
                    self.refTraj.getData()[self._nDt0:
                                            self._nDt0+self._nDt+1],
                    param[1].getData(), param[2].getData(),
                    param[3].getData(),param[4].getData(),
                    param.nu, param.nuN,
                    fullTraj=True)

        tReal=self._nDt*self.dt
        traj.putData(trajData)
        traj.incrmTReal(finished=True, tReal=tReal, t0=t0)

        return traj


    #------------------------------------------------------

    def __kdvTLMProp_Fortran_Adj(self, pert, traj, t0=0.):

        # Local variables names
        grid=self.grid
        param=self.__t0AdjustParam(self.param, t0=t0)

        trajData=fKdV.fKdVTLMPropagatorAdj(
                grid.N, grid.Ntrc, grid.L, self.dt, self._nDt, param.nDt,
                pert,
                self.refTraj.getData()[self._nDt0:self._nDt0+self._nDt+1],
                param[1].getData(), param[2].getData(),
                param[3].getData(), param[4].getData(), 
                param.nu, param.nuN, 
                fullTraj=True)

        tReal=self._nDt*self.dt
        traj.putData(trajData)
        traj.incrmTReal(finished=True, tReal=tReal, t0=t0)

        return traj

    #------------------------------------------------------

    def __kdvTLMSingularOp_Fortran(self, pert, t0=0.):

        # Local variables names
        grid=self.grid
        param=self.__t0AdjustParam(self.param, t0=t0)

        LAdjLx=fKdV.fKdVTLMSingularOp(
                grid.N, grid.Ntrc, grid.L, self.dt, self._nDt, param.nDt,
                pert,
                self.refTraj.getData()[self._nDt0:self._nDt0+selfnDt+1],
                param[1].getData(), param[2].getData(), param[3].getData(),
                param[4].getData(), param.nu, param.nuN)

        tReal=2.*self._nDt*self.dt

        self.incrmTReal(finished=True, tReal=tReal)
        return LAdjLx

        
        
    #------------------------------------------------------

    def __kdvTLMGradTest_Fortran(self, ic, dt, nDt, maxPow, t0=0.):
        
        grid=self.grid
        param=self.__t0AdjustParam(self.param, t0=t0)
   

        fKdV.fKdVTestGradient(grid.N, grid.Ntrc, grid.L, 
                dt, nDt, param.nDt, maxPow, ic,
                param[1].getData(), param[2].getData(), 
                param[3].getData(), param[4].getData(),
                param.nu, param.nuN,
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
    from pseudoSpec1D import gradientTest
    
    testGrad=True
    testAdjoint=True
    testNDtInt=True
    
    Ntrc=144
    maxA=10.
    nDt=100
    grid=PeriodicGrid(Ntrc)
    param=Param(grid)#, rho=-0.1*gauss(grid.x, 0., 10.))
    dt=0.001
    if dt > dtStable(param, maxA): raise RuntimeError()
    tInt=nDt*dt
    freq=7

    
    #----| Reference trajectory |-----------------
    u0=rndSpecVec(grid, seed=0)
    M=kdvLauncher(param, dt=dt)
    u=M.integrate(u0, tInt)

    L=kdvTLMLauncher(param, traj=u)

    if testGrad:
        #----| Gradient test |------------------------
        print("\nGradient test\n")
        L.gradTest(M)
        L.gradTestFortran(u0, dt, nDt)





    if testAdjoint:
        #----| Adjoint testing |----------------------
        print("Testing adjoint validity")
    
        dx=rndSpecVec(grid, amp=0.1, seed=1)
        dy=rndSpecVec(grid, amp=0.1, seed=2)
    
        Ldy=L.integrate(dy).final
        Adx=L.adjoint(dx).ic
        print("<dx, Ldy> = %+.15g"%np.dot(dx, Ldy))
        print("<Adx, dt> = %+.15g"%np.dot(Adx, dy))
        print("<dx, Ldy>-<Adx, dt> = %+.15g"%(
                    np.dot(dx, Ldy)-np.dot(Adx, dy)))


        #----| partial trajectory time |----
        print("\nTestting adjoint validity for partial integration")
        Ldy=L.integrate(dy, tInt=tInt/3., t0=tInt/2.).final
        LAdj_x=L.adjoint(dx, tInt=tInt/3., t0=tInt/2.).ic
        print("<dx, Ldy>-<L*dx, dt> = %+.15g"%(
                    np.dot(dx, Ldy)-np.dot(LAdj_x, dy)))
    
    
        #----| step integrations |----------
        print("\nTestting adjoint validity"
                +" for successive step integrations")
        L1=kdvTLMLauncher(param, traj=u)
        L2=kdvTLMLauncher(param, traj=u)
    
        L1dy=L1.integrate(dy, tInt=tInt/3., t0=0.).final
        L2L1dy=L2.integrate(L1dy, tInt=tInt/2., t0=tInt/3.).final
    
        A1dx=L1.adjoint(dx, tInt=tInt/3., t0=0.).ic
        A2dx=L2.adjoint(dx, tInt=tInt/2., t0=tInt/3.).ic
        A1A2dx=L1.adjoint(A2dx, tInt=tInt/3., t0=0.).ic
    
    
        print("<dx,L1dy>-<A1dx, dy>     = %+.15g"%(np.dot(dx, L1dy)\
                                                -np.dot(A1dx, dy)))
        print("<dx,L2L1dy>-<A2dx, L1dy> = %+.15g"%(np.dot(dx, L2L1dy)\
                                                -np.dot(A2dx, L1dy)))
        print("<dx,L2L1dy>-<A1A2dx, dy> = %+.15g"%(np.dot(dx, L2L1dy)\
                                                -np.dot(A1A2dx, dy)))


    if testNDtInt:
        dx=rndSpecVec(grid, amp=0.1, seed=1)
        nDtList=[i*u.nDt/freq for i in xrange(1,freq+1)]

        if testAdjoint:
            #----| Sequential integration adjoint test |--
            d_dy={}
            print("\nSequential adjoint test with d_nDtInt()")
            
            for j in nDtList:
                d_dy[j]=rndSpecVec(grid, amp=0.1, seed=j)

            d_Hx=L.d_nDtInt(dx, nDtList)
            Ay=L.d_nDtIntAdj(d_dy)
            
            Hx_y=0.
            for j in nDtList:
                Hx_y+=np.dot(d_Hx[j], d_dy[j])

            x_Ay=np.dot(dx, Ay)
            print(Hx_y, x_Ay, Hx_y-x_Ay)

        if testGrad:
            print("\n grad test between d_nDtInt()")
    
            def fct(x0, model, tlm, nDtL):
                d_x=model.d_nDtInt(x0, nDtL)
                J=0.
                for t in d_x.keys():
                    J+=0.5*np.dot(d_x[t], d_x[t])
                return J
            def gradFct(x0, model, tlm, nDtL):
                x=model.integrate(x0, dt*nDtL[-1])
                tlm.reference(x)
                d_x=model.d_nDtInt(x0, nDtL)
                return tlm.d_nDtIntAdj(d_x)
    
    
            gradientTest(dx, fct, gradFct, args=(M, L, nDtList))
