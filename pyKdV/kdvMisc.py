import numpy as np

def soliton(x, x0, amp=1., alpha=0., beta=1., gamma=-1.):
    """ 
    Eigen solution of the  system

    """
    sgn=beta*gamma/abs(beta*gamma)
    amp=abs(amp)*sgn
    discr=np.sqrt(amp*beta/(12.*gamma))*(x-x0)
    fct=amp*np.cosh(discr)**(-2)

    return fct 

def gauss(x, x0, sig):
    return np.exp(-((x-x0)**2)/(2*sig**2))


