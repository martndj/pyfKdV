import numpy as np
import random as rnd
from periodicGrid import PeriodicGrid
#-----------------------------------------------------------

def rndSpecVec(g, Ntrc=None, Nmin=0, amp=1., sig=1., seed=None):
    """
    Pseudo random signal genrator
        (resolution independant : will generate the same signal
            independantly of g.Ntrc)

        rndSpecVec(g, Ntrc=None, amp=1., sig=1., seed=0.848241945)

        g       :   grid <PeriodicGrid>
        Ntrc    :   truncature (high frequency) <int>
        Nmin    :   low frequency truncatrue <int>
        amp     :   amplitude <float>
        sig     :   standard deviation <float>
        seed    :   random generator seed <float>

    """
    if not isinstance(g, PeriodicGrid):
        raise Exception("g <PeriodicGrid>")
    rnd.seed(seed)
    y=np.zeros(g.N, dtype='complex')
    if Ntrc==None:
        Ntrc=g.Ntrc
    y[0]=rnd.gauss(0., sig)
    for i in xrange(Nmin+1,Ntrc):
        y[i]=rnd.gauss(0., sig)+1j*rnd.gauss(0., sig)
        y[-i]=y[i].conjugate()
    y=np.fft.ifft(y).real
    y=y/np.max(np.abs(y))*amp
    return y

#-----------------------------------------------------------

def gauss(x, x0, sig):
    """
    Au gaussian distribution
        
        gauss(x, x0, sig)

        x       :   grid space <numpy.ndarray>
        x0      :   center <float>
        sig     :   width <float>
    """
    return np.exp(-((x-x0)**2)/(2*sig**2))

#-----------------------------------------------------------

def bumpSchwartz(x):
    """
    L. Schwartz infinitly differentiable distribution

    <todo> vectorize it? (if isinstance(x, np.ndarray : return...)
    """
    if np.abs(x)>1.:
        return 0.
    else:
        return np.exp(1.)*np.exp(-1./(1.-x**2))

