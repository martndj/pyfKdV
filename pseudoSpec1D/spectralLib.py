import numpy as np
import numpy.fft as fft
from periodicGrid import PeriodicGrid 

#----------------------------------------------------------

def conjT(M):
    """
    Transpose conjugate of an array
    """
    return M.conjugate().T

#----------------------------------------------------------

def fft2d(M):
    """
    Fast Fourier Transform of a 2nd order tensor
    """
    return conjT(np.fft.fft(conjT(np.fft.fft(M))))    


def ifft2d(M):
    """
    Inverse Fast Fourier Transform of a 2nd order tensor
    """
    return conjT(np.fft.ifft(conjT(np.fft.ifft(M)))) 

#----------------------------------------------------------

def fftOrder(grid):
    """
    Return the FFT order of frequency index
        
        fftOrder(grid)

        grid    :   <PeriodicGrid>
    """
    if not isinstance(grid, PeriodicGrid):
        raise TypeError("grid <PeriodicGrid>")
    n=np.zeros(grid.N)
    for i in xrange(grid.N):
        if i <= (grid.N-1)/2:
            n[i]=i
        else:
            n[i]=i-grid.N
    return n

#----------------------------------------------------------

def diffMatrixDiag(order, grid):
    """
    Return the diagonal of the spectral differentiation matrix

        diffMatrixDiag(order, grid)

        order   :   order of differentiation <int>
        grid    :   <PeriodicGrid>
    """
    if not isinstance(grid, PeriodicGrid):
        raise TypeError("grid <PeriodicGrid>")
    D=np.zeros(grid.N, dtype=complex)
    n=fftOrder(grid)
    for i in xrange(grid.N):
        # truncature
        if np.abs(n[i])<=grid.Ntrc:
            D[n[i]]=(1.j*n[i]*2.*np.pi/grid.L)**order
    return D

#----------------------------------------------------------

def specDiff(f,grid, order=1):
    """
    Spectral differentiation

        specDiff(f,grid, order=1)

        f       :   input <numpy.ndarray>
        grid    :   <PeriodicGrid>
        order   :   order of differentiation <int>
    """
    if not isinstance(grid, PeriodicGrid):
        raise TypeError("grid <PeriodicGrid>")
    if not isinstance(f, np.ndarray):
        raise TypeError("f <numpy.ndarray>")
    if not (f.ndim==1 and len(f)==grid.N):
        raise TypeError("f.shape=(grid.N)")
    tf=fft.fft(f)
    tfdf=diffMatrixDiag(order, grid)*tf
    df=(fft.ifft(tfdf)).real
    return df

def specDiff_Adj(f, grid, order=1):
    """
    Spectral differentiation adjoint

        specDiff_Adj(f,grid, order=1)

        f       :   input <numpy.ndarray>
        grid    :   <PeriodicGrid>
        order   :   order of differentiation <int>
    """
    if not isinstance(grid, PeriodicGrid):
        raise TypeError("grid <PeriodicGrid>")
    if not isinstance(f, np.ndarray):
        raise TypeError("f <numpy.ndarray>")
    if not (f.ndim==1 and len(f)==grid.N):
        raise ValueError("f.shape=(grid.N)")
    tf=fft.fft(f)
    tfdf=(diffMatrixDiag(order, grid)).conjugate()*tf
    df=(fft.ifft(tfdf)).real.copy(order='C')
    return df
    
#----------------------------------------------------------

def specFilt(f,grid, Ntrc=None, Nlow=0.):
    """ 
    Spectral truncature
        <!> Auto-adjoint

        specFilt(f,grid,Ntrc=None)
        
        f       :   input <numpy.ndarray>
        grid    :   <PeriodicGrid>
        Ntrc    :   truncature <int>
    """
    if not isinstance(grid, PeriodicGrid):
        raise TypeError("grid <PeriodicGrid>")
    if not isinstance(f, np.ndarray):
        raise TypeError("f <numpy.ndarray>")
    if not (f.ndim==1 and len(f)==grid.N):
        raise TypeError("f.shape=(grid.N)")
    if Ntrc==None:
        Ntrc=grid.Ntrc
    
    tf=fft.fft(f)
    n=fftOrder(grid)
    for i in xrange(grid.N):
        # truncature
        if np.abs(n[i])>=Ntrc or np.abs(n[i])<Nlow:
            tf[i]=0.
    f=(fft.ifft(tf)).real.copy(order='C')
    return f
    
#----------------------------------------------------------

def specTrunc(w, grid, Ntrc=None):
    """ 
    Spectral truncature

        specTrunc(w,grid,Ntrc=None)
        
        w       :   spectral input <numpy.ndarray>
        grid    :   <PeriodicGrid>
        Ntrc    :   truncature <int>
    """
    if not isinstance(grid, PeriodicGrid):
        raise TypeError("grid <PeriodicGrid>")
    if not isinstance(w, np.ndarray):
        raise TypeError("w <numpy.ndarray>")
    if not (w.ndim==1 and len(w)==grid.N):
        raise ValueError("w.shape=(grid.N)")
    if Ntrc==None:
        Ntrc=grid.Ntrc
    n=fftOrder(grid)
    for i in xrange(grid.N):
        # truncature
        if np.abs(n[i])>=Ntrc:
            w[i]=0.
    return w


#----------------------------------------------------------

def specShift(f, grid, dx):
    """
    Shift a periodic signal.
    """

    if not isinstance(grid, PeriodicGrid):
        raise TypeError("grid <PeriodicGrid>")

    kSpace=fftOrder(grid)
    tf=fft.fft(f)
    w=tf*np.exp(-2.j*np.pi*dx*kSpace/grid.L)
    return fft.ifft(w).real
