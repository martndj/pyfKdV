import numpy as np
from pseudoSpec1D import PeriodicGrid, specFilt


#-----------------------------------------------------------

def soliton(x, x0=0., amp=1., alpha=0., beta=1., gamma=-1.):
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

#-----------------------------------------------------------

def cSoliton(amp=1., alpha=0., beta=1., gamma=-1.):
    """
    Speed of the KdV eigensolution
    
        cSoliton(amp=1., alpha=0., beta=1., gamma=-1.)
        
        amp     :   amplitude <float>

        alpha, beta, gamma
                :   KdV parameters (constants) <float>

    """
    return np.abs(beta)/3.*beta*gamma/(np.abs(beta*gamma))*amp+alpha

def ampSolStat(alpha=0., beta=1., gamma=-1.):
    return -3.*alpha/np.abs(beta)*beta*gamma/(np.abs(beta*gamma))

#-----------------------------------------------------------

def dtStable(grid, param, maxA, dtMod=0.7):
    """
    Stable time incremement

        dtStable(grid, param, maxA)

        grid    :   <Grid>
        maxA    :   expected maximum amplitude <float>
        param   :   KdV parameters <Param>
    """


    maxK=2.0*np.pi*grid.Ntrc/grid.L

    dt=None 
    discr=np.zeros(shape=grid.N)
    for t in xrange(param.nDt+1):
        rhoPos=np.zeros(shape=grid.N)
        for i in xrange(grid.N):
            if param.rho[t][i] > 0. :
                rhoPos[i]=param.rho[t][i]
        discr=(param.gamma[t]*maxK**2-np.abs(param.alpha[t])
                -np.abs(param.beta[t])*maxA)*maxK
        dtTmp=1./np.max(np.sqrt(discr**2+rhoPos**2))
        if dt==None or dtTmp<dt:
            dt=dtTmp

    return dtMod*dt
