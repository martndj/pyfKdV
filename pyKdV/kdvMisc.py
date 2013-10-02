import numpy as np
import random as rnd
from pseudoSpec1D import PeriodicGrid, specFilt


def rndSpecVec(g, Ntrc=None, amp=1., seed=0.848241945):
    """
    Pseudo random signal genrator
    (resolution independant : will generate the same signal
        independantly of g.Ntrc)
    """
    rnd.seed(seed)
    y=np.zeros(g.N, dtype='complex')
    if Ntrc==None:
        Ntrc=g.Ntrc
    y[0]=rnd.gauss(0., amp)
    for i in xrange(1,Ntrc):
        y[i]=rnd.gauss(0., 1.)+1j*rnd.gauss(0., 1.)
        y[-i]=y[i].conjugate()
    y=np.fft.ifft(y).real
    y=y/np.max(np.abs(y))*amp
    return y

def soliton(x, x0, amp=1., alpha=0., beta=1., gamma=-1.):
    """ 
    Eigen solution of the non augmented KdV system (rho==0)

    """
    sgn=beta*gamma/abs(beta*gamma)
    amp=abs(amp)*sgn
    discr=np.sqrt(amp*beta/(12.*gamma))*(x-x0)
    fct=amp*np.cosh(discr)**(-2)

    return fct 

def gauss(x, x0, sig):
    return np.exp(-((x-x0)**2)/(2*sig**2))


