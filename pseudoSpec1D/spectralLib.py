import numpy as np
import numpy.fft as fft


class spectralLibError(Exception):
    pass
#----------------------------------------------------------

def conjT(M):
    return M.conjugate().T

#----------------------------------------------------------

def fft2d(M):
    return conjT(np.fft.fft(conjT(np.fft.fft(M))))    


def ifft2d(M):
    #return g.N*np.fft.ifft2(M)
    return conjT(np.fft.ifft(conjT(np.fft.ifft(M)))) 

#----------------------------------------------------------

def fftOrder(grid):
    n=np.zeros(grid.N)
    for i in xrange(grid.N):
        if i <= (grid.N-1)/2:
            n[i]=i
        else:
            n[i]=i-grid.N
    return n

#----------------------------------------------------------

def diffMatrixDiag(order, grid):
    D=np.zeros(grid.N, dtype=complex)
    n=fftOrder(grid)
    for i in xrange(grid.N):
        # truncature
        if np.abs(n[i])<=grid.Ntrc:
            D[n[i]]=(1.j*n[i]*2.*np.pi/grid.L)**order
    return D

#----------------------------------------------------------

def specDiff(f,grid, order=1):
    tf=fft.fft(f)
    tfdf=diffMatrixDiag(order, grid)*tf
    df=(fft.ifft(tfdf)).real
    return df

def specDiff_Adj(f, grid, order=1):
    tf=fft.fft(f)
    tfdf=(diffMatrixDiag(order, grid)).conjugate()*tf
    df=(fft.ifft(tfdf)).real.copy(order='C')
    return df
    
#----------------------------------------------------------

def specFilt(f,grid,Ntrc=None):
    """ 
    In-place spectral truncature
    
        <!> Auto-adjoint
    """
    if Ntrc==None:
        Ntrc=grid.Ntrc
    tf=fft.fft(f)
    n=fftOrder(grid)
    for i in xrange(grid.N):
        # truncature
        if np.abs(n[i])>Ntrc:
            tf[i]=0.
    f=(fft.ifft(tf)).real.copy(order='C')
    return f
    
#----------------------------------------------------------

def specTrunc(w, grid, Ntrc=None):
    if Ntrc==None:
        Ntrc=grid.Ntrc
    n=fftOrder(grid)
    for i in xrange(grid.N):
        # truncature
        if np.abs(n[i])>Ntrc:
            w[i]=0.
    return w
