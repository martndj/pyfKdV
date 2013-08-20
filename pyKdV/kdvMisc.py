import numpy as np

def soliton(param, x0, amp):
    """ 
    Eigen solution of the  system

    """
    N=param.grid.N
    fct=np.zeros(N, dtype=np.float64)

    sgn=param[2]*param[3]/np.abs(param[2]*param[3])
    amp=np.abs(amp)*sgn
    discr=np.sqrt(amp*param[2]/(12.*param[3]))*(param.grid.x-x0)
    fct=amp*np.cosh(discr)**(-2)

    return fct 

def gauss(param, x0, sig):
    return np.exp(-((param.grid.x-x0)**2)/(2*sig**2))


