import numpy as np
import random as rnd
from pseudoSpec1D import PeriodicGrid, specFilt


def rndSpecVec(g, Ntrc=None, amp=1., seed=0.848241945):
    """
    Pseudo random signal genrator
        (resolution independant : will generate the same signal
            independantly of g.Ntrc)

        rndSpecVec(g, Ntrc=None, amp=1., seed=0.848241945)

        g       :   grid <PeriodicGrid>
        Ntrc    :   truncature <int>
        amp     :   amplitude <float>
        seed    :   random generator seed <float>

    """
    if not isinstance(g, PeriodicGrid):
        raise Exception("g <PeriodicGrid>")
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
        
        soliton(x, x0, amp=1., alpha=0., beta=1., gamma=-1.)

        x       :   grid space <numpy.ndarray>
        x0      :   position <float>
        amp     :   amplitude <float>

        alpha, beta, gamma
                :   KdV parameters (constants) <float>

    """
    sgn=beta*gamma/abs(beta*gamma)
    amp=abs(amp)*sgn
    discr=np.sqrt(amp*beta/(12.*gamma))*(x-x0)
    fct=amp*np.cosh(discr)**(-2)

    return fct 

def cSoliton(amp=1., alpha=0., beta=1., gamma=-1.):
    """
    Speed of the KdV eigensolution
    
        cSoliton(amp=1., alpha=0., beta=1., gamma=-1.)
        
        amp     :   amplitude <float>

        alpha, beta, gamma
                :   KdV parameters (constants) <float>

    """
    return np.abs(beta)/3.*beta*gamma/(np.abs(beta*gamma))*amp

def gauss(x, x0, sig):
    """
    Au gaussian distribution
        
        gauss(x, x0, sig)

        x       :   grid space <numpy.ndarray>
        x0      :   center <float>
        sig     :   width <float>
    """
    return np.exp(-((x-x0)**2)/(2*sig**2))

def bumpSchwartz(x):
    """
    L. Schwartz infinitly differentiable distribution
    """
    if np.abs(x)>1.:
        return 0.
    else:
        return np.exp(-1./(1.-x**2))

def dtStable(grid, param, maxA, dtMod=0.7):
    """
    Stable time incremement

        dtStable(grid, param, maxA)

        grid    :   <Grid>
        maxA    :   expected maximum amplitude <float>
        param   :   KdV parameters <Param>
    """

    minRho=param[4].min()
    maxRho=param[4].max()
    if np.abs(minRho)>np.abs(maxRho):
        maxAbsRho=np.abs(minRho)
    else:
        maxAbsRho=np.abs(maxRho)

    maxK=2.0*np.pi*grid.Ntrc/grid.L
    denom=np.zeros(shape=grid.N)
    denom=np.sqrt((param[3].max()*maxK**3-param[1].min()*maxK
                        -param[2].min()*maxA*maxK)**2
                   +maxAbsRho**2)

    dt=1./denom
    return dtMod*dt
