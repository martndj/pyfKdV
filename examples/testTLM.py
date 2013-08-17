import numpy as np
from pyfKdV import *
import matplotlib.pyplot as plt

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
fPert=tLauncher.integrate(tInt, fullPertTraj=True)

aLauncher=TLMLauncher(grid, param, traj, fPert)
aPert=aLauncher.adjoint(tInt, fullPertTraj=True)

sLauncher=TLMLauncher(grid, param, traj, pert)
sPert=sLauncher.singularOp(tInt)

plt.plot(grid.x, fPert)
plt.plot(grid.x, aPert)
plt.plot(grid.x, sPert)
    
plt.show()
