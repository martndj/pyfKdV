import numpy as np
from pyfKdV import *
import matplotlib.pyplot as plt

grid=SpectralGrid(150,300.)
tInt=50.
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
axe=traj.waterfall()
plt.show()
