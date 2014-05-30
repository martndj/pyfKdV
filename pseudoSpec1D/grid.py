import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec, SubplotSpec

class Grid(object):
    """
    Grid class
    1+1D partial differential system

    Grid(N, L, centered=True)

        N       :   number of grid points <int>
        L       :   domain lenght <float>
        centered:   centered grid: x in (-L/2, L/2)
                    not centered grid: x in (0., L)
    """
    
    
    def __init__(self, N, L, centered=True):
        if  not(type(N) is int) or not(isinstance(L, (float,int))):
            raise TypeError("Grid(N <int>| L <float|int>)")
        
    
        self.N=N
        self.L=L
        self.centered=centered

        self.__setGridPoints()
    
    #-------------------------------------------------------
    #----| Public methods |---------------------------------
    #-------------------------------------------------------
    def max(self):
        return np.max(self.x)
       
    #-------------------------------------------------------

    def min(self):
        return np.min(self.x)

    #-------------------------------------------------------

    def pos2Idx(self, pos):
        """
        Convert space position to grid index

            Grid.pos2Idx(pos)

            pos :   array of positions <numpy.ndarray>
        """

        def findIdx(posValue):
            whereVec=np.where(self.x>=posValue)
            if len(whereVec[0])==0 or self.min()>posValue:
                raise ValueError(
                    "Position outside of the grid:\n %f not in [%f,%f]"%(
                                posValue, self.min(), self.max()))
            else:
                return np.min(whereVec)

        if isinstance(pos, (list, np.ndarray)):
            if isinstance(pos, np.ndarray):
                if pos.ndim<>1:
                    raise ValueError("pos.ndim=1")
            N=len(pos)
            idx=np.empty(N, dtype=int)
            for i in xrange(N):
                idx[i]=findIdx(pos[i])
        elif isinstance(pos, (float, int)):
            idx=np.empty(1, dtype=int)
            idx[0]=findIdx(pos)
        else:
            raise TypeError("pos <list|numpy.ndarray|float>")
        return idx 

    def where(self, pos):
        return self.pos2Idx(pos)

    #-------------------------------------------------------

    def squareNorm(self, field, metric=None, metricArgs=()):
        if not isinstance(field, np.ndarray):
            raise TypeError("field <numpy.ndarray>")
        if field.ndim<>1 or field.shape[-1]<>self.N:
            raise ValueError("field icompatible dimensions")

        if metric==None:
            return np.dot(field,field)*self.dx
        elif metric==np.infty:
            return np.abs(field).max()**2
        elif isinstance(metric, (float, int)):
            return np.dot(field, field)*metric*self.dx
        elif isinstance(metric, np.ndarray):
            if ((not (metric.ndim==1 or metric.ndim==2))
                    or metric.shape[-1]<>self.N):
                raise ValueError("metric icompatible dimensions") 
            return np.dot(field, np.dot(metric, field))*self.dx
        elif callable(metric):
            return metric(field, *metricArgs)
        else:
                raise TypeError("metric <numpy.ndarray>")
    
    #-------------------------------------------------------

    def norm(self, field, metric=None, metricArgs=()):
        return np.sqrt(self.squareNorm(field, metric=metric, metricArgs=metricArgs))
    
    #-------------------------------------------------------

    def zeros(self):
        return np.zeros(self.N)
    
    #-------------------------------------------------------

    def plot(self, field,  axe=None, xlabel=r'$x$',  **kwargs):
        axe=self._checkAxe(axe)
        axe.plot(self.x, field,  **kwargs)
        axe.set_xlabel(xlabel)
        return axe


    #-------------------------------------------------------

    def gradient(self, field):
        return np.gradient(field, self.dx)
    

    #-------------------------------------------------------
    #----| Private methods |--------------------------------
    #-------------------------------------------------------

    def __setGridPoints(self):
        self.dx=self.L/(self.N)
        if self.centered:
            self.x=np.linspace(-self.L/2.,self.L/2-self.dx,self.N)
        else:
            self.x=np.linspace(0.,self.L-self.dx,self.N)
            

    #-------------------------------------------------------

    def _checkAxe(self, axe):
        if axe==None:
            axe=plt.subplot(111)
        if isinstance(axe, int):
            if len(str(axe))<>3:
                raise ValueError(
                    "Single argument to subplot must be a 3-digit integer")
            axe=plt.subplot(axe)
        elif isinstance(axe,SubplotSpec):
            axe=plt.subplot(axe)
        elif isinstance(axe,Axes):
            pass
        else:
            raise TypeError(
            "axe < matplotlib.axes.Axes | matplotlib.gridspec.GridSpec >")
        return axe

    #-------------------------------------------------------
    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="----| Grid |---------------------------\n"
        output+="| N=%d\n"%(self.N)
        output+="| dx=%-23.15E\n"%(self.dx)
        output+="| L=%-23.15E\n"%(self.L)
        output+="---------------------------------------"
        return output


    def __eq__(self, grid2):
        value=(self.N==grid2.N and self.L==grid2.L)
        return value 
#--------------------------------------------------------------------
#####################################################################
#--------------------------------------------------------------------

if __name__=="__main__":
    g=Grid(200,100.)
    print(g) 

