import numpy as np
import progressbar

from spectralGrid import *
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

    def __init__(self, grid, param, traj, pert):
        if not(isinstance(grid, SpectralGrid)):
            raise TLMLauncherError("grid <SpectralGrid>")
        self.grid=grid

        if not(isinstance(grid, SpectralGrid)):
            raise TLMLauncherError("param <Param>")

        self.param=param
        
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

        # Propagator
        self.propagator=self.__kdvTLMProp_Fortran
        self.propagator_Adj=self.__kdvTLMProp_Fortran_Adj
    
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

    def integrate(self, tInt, progress=False, fullPertTraj=False,
                  adjoint=False, resetTReal=True):
        

        dt=self.refTraj.dt
        nDt=int(tInt/dt)
        if resetTReal: self.tReal=0.

        if fullPertTraj:
            self.fullPertTraj=True
            self.pertTraj=Trajectory(self.grid)
            self.pertTraj.initialize(self.pert, tInt, dt)

        # calling the Tangent model
        #   defined in the herited classes (initial and final perturbations)

        if not adjoint:
            integratedPert=self.propagator(dt, nDt,
                progress=progress, fullPertTraj=fullPertTraj)
        else:
            integratedPert=self.propagator_Adj(dt, nDt, 
                progress=progress, fullPertTraj=fullPertTraj)

        return integratedPert
    
    #------------------------------------------------------

    
    def singular(self, tInt, progress=False):
        
        pert=self.integrate(tInt, progress=progress)
        pert=self.integrate(tInt, progress=progress,
                            resetTReal=False,adjoint=True)

        return pert

    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------

    #----| Propagator (and adjoint) |-----------------------
    #------------------------------------------------------

    def __kdvTLMProp_Fortran(
            self, dt, nDt, progress=False,fullPertTraj=True):

        # Local variables names
        grid=self.grid

        


        pertTrajArr=fKdV.fKdVTLMPropagator(
                    grid.N, grid.Ntrc, grid.L, dt, nDt,
                    self.pert, self.refTraj.getData(),
                    param[1], param[2], param[3], param[4], param[0])

        tReal=nDt*dt
        fPert=np.array(pertTrajArr[nDt])
        self.incrmTReal(finished=True, tReal=tReal)
        if True:
            self.pertTraj.putData(pertTrajArr)
            self.pertTraj.incrmTReal(finished=True, tReal=tReal)

        return fPert


    #------------------------------------------------------

    def __kdvTLMProp_Fortran_Adj(
            self, dt, nDt, progress=False,fullPertTraj=True):

        # Local variables names
        grid=self.grid

        


        adjTrajArr=fKdV.fKdVTLMPropagatorAdj(
                    grid.N, grid.Ntrc, grid.L, dt, nDt,
                    self.pert, self.refTraj.getData(),
                    param[1], param[2], param[3], param[4], param[0])

        tReal=nDt*dt
        aPert=np.array(adjTrajArr[nDt])
        self.incrmTReal(finished=True, tReal=tReal)
        if True:
            self.pertTraj.putData(adjTrajArr)
            self.pertTraj.incrmTReal(finished=True, tReal=tReal)

        return aPert

#--------------------------------------------------------------------
#====================================================================
#--------------------------------------------------------------------



if __name__=='__main__':
    import matplotlib.pyplot as plt
    from kdvLauncher import *

    grid=SpectralGrid(150,300.)
    tInt=3.
    maxA=2.

    def funcNulle(x):
        return 0.
    def funcNeg(x):
        return -1.
    def funcPos(x):
        return 1.

    param=Param(grid, funcNulle, funcNulle, funcPos, funcNeg, funcNulle)

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


    pert=0.1*np.exp(-(5.*np.pi*(grid.x-10)/grid.L)**2)

    tLauncher=TLMLauncher(grid, param, traj, pert)

#    print('Testing differential operator adjoint validity')
#    testAdj=tLauncher.testKdvPseudoSpecTLM_Adj(p4Term=True)
#    if testAdj[3]:
    if True:
#        print('Validating the test:'+str(testAdj[2]))
        fPert=tLauncher.integrate(tInt, fullPertTraj=True)
        aLauncher=TLMLauncher(grid, param, traj, fPert)
        aPert=aLauncher.integrate(tInt, fullPertTraj=True, adjoint=True)
        plt.plot(grid.x, fPert)
        plt.plot(grid.x, aPert)
    else:
        print('Test failed!:'+str(testAdj[2]))
        
    plt.show()