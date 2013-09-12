import numpy as np
import copy

import matplotlib.pyplot as plt
from matplotlib import collections

from spectralGrid import SpectralGrid

class Trajectory(object):
    """
    Trajectory for 1+1D partial differential system

        Light version for pyKdV5
            * no I\O methods
            * fully encapsulated
    """

    __version__='Light version for pyKdV5'
    __author__='Martin Deshaies-Jacques <deshaies.martin@sca.uqam.ca>'

    class TrajectoryError(Exception):
        pass

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------


    def __init__(self, grid):
        """
        Empty Trajectory instance constructor
        """
        # Grid
        if not (isinstance(grid, SpectralGrid)):
            raise self.TrajectoryError(
                  "grid <SpectralGrid>")
        self.grid=grid

        # Time attributes
        self.tInt=0.
        self.tReal=0.
        self.nDt=0
        self.dt=0.
        
        # Status Attributes
        self.isIntegrated=False
        self.isInitialised=False

        # Data attributes
        self.__data=None
        self.ic=None
        self.A=None
        self.A2=None
        self.time=None

    #------------------------------------------------------
    #----| general public methods |------------------------
    #------------------------------------------------------

    
    def copy(self):
        return copy.deepcopy(self)

        
    #-------------------------------------------------------
    
    def initialize(self, ic, nDt, dt):
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise self.TrajectoryError("ic.shape<>(grid.N,)")
        self.__allocate(nDt, dt)
        self.ic=ic
        self.isInitialised=True

    
    #-------------------------------------------------------
    
    def incrmTReal(self, finished=False, tReal=None):
        if tReal<>None:
            self.tReal=tReal
        if finished:
            self.isIntegrated=True
            self.time=np.linspace(0.,self.tReal, self.nDt+1)
        else:
            self.tReal+=self.dt

    #-------------------------------------------------------

    def final(self):
        if self.isIntegrated:
            return self.__data[self.nDt]
        else:
            raise self.TrajectoryError("Trajectory not integrated")

    #-------------------------------------------------------

    def norm2(self, ret=True):
        if not self.isIntegrated:
            raise self.TrajectoryError("Trajectory not integrated")

        self.A2=np.zeros(self.nDt+1)
        for i in xrange(self.nDt+1):
            self.A2[i]=self.__SquareNorm(self.__data[i])
        if ret:
            return self.A2

    #-------------------------------------------------------

    def norm(self, ret=True):
        if not self.isIntegrated:
            raise self.TrajectoryError("Trajectory not integrated")

        self.A=np.zeros(self.nDt+1)
        if not (self.A2==None):
            self.A=np.sqrt(self.A2)
        else:
            for i in xrange(self.nDt+1):
                self.A[i]=np.sqrt(self.__SquareNorm(self.__data[i]))
        if ret:
            return self.A

    #-------------------------------------------------------

    def getData(self):
        return self.__data

    #-------------------------------------------------------

    def putData(self, data):
        if not self.isInitialised :
            raise self.TrajectoryError("Trajectory not initialised")
        if (data.shape<>(self.nDt+1, self.grid.N)):
            raise self.TrajectoryError("Incompatible data affectation")
        self.__data=data
        self.shape=data.shape


    #-------------------------------------------------------
    #----| Fortran compatible I/O Methods |-----------------
    #-------------------------------------------------------

    def write(self, trajFile):

        fich=open(trajFile, 'w')
        
        # Header
        fich.write('%5i%5i\n'%(self.grid.Ntrc, self.grid.N))
        fich.write(' %-23.15E%-23.15E\n'%(self.grid.dx, self.grid.L))
        fich.write(' %-23.15E%-23.15E%-23.15E\n'%(self.dt, self.tInt,
                                                 self.tReal))
        fich.write('%10i\n'%(self.nDt))
    
        # Data
        if self.nDt==0:
            for j in xrange(self.grid.N):
                fich.write('%23.15E\n'%self.data[j])
        elif self.nDt>0:
            for i in xrange(self.nDt+1):
                for j in xrange(self.grid.N):
                    fich.write('%23.15E\n'%self.data[i,j])
    
        fich.close()
        print("  >> ")+trajFile



    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------

    def __SquareNorm(self, z):
        return np.dot(z,z)*self.grid.dx

    #-------------------------------------------------------

    def __allocate(self, nDt, dt):
        self.nDt=nDt
        self.dt=dt
        self.tInt=self.nDt*self.dt
        
        self.__data=np.zeros(shape=(self.nDt+1,self.grid.N))
        self.ic=np.zeros(self.grid.N)



    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="====| Trajectory |=================================\n"
        output+=self.grid.__str__()
        #if not self.isIntegrated:
        #    output+="\n| Trajectory not integrated"
        #else:
        output+="\n| nDt=%d"%self.nDt
        output+="\n| tInt=%-23.15E"%self.tInt
        output+="\n| dt=%-23.15E"%self.dt
        output+="\n| tReal=%-23.15E"%self.tReal
        output+="\n===================================================\n"
        return output



    #-------------------------------------------------------

    def __getitem__(self, idx):
        if idx>self.nDt:
            raise IndexError()
        return self.__data[idx]


    #-------------------------------------------------------
    
    def __setitem__(self, idx, state):
        self.__data[idx]=state

    #-------------------------------------------------------

    def __sub__(self, traj2):
        if (not isinstance(traj2, Trajectory)):
            raise self.TrajectoryError("Operation on Trajectory objects")

        if (self.dt != traj2.dt):
            raise self.TrajectoryError("Incompatible time increment")
        elif (not self.grid == traj2.grid):
            print("%d %d"%(self.grid.Ntrc,traj2.grid.Ntrc))
            print("%g %g"%(self.grid.dx,traj2.grid.dx))
            print("%g %g"%(self.grid.L,traj2.grid.L))
            print(self.grid.Ntrc==traj2.grid.Ntrc and 
                  self.grid.dx==traj2.grid.dx and 
                  self.grid.L==traj2.grid.L)
            raise self.TrajectoryError("Incompatible grids")

        trajSub=__data2Traj(self.grid, self.__data-traj2.__data,self.dt)
        return trajSub

    #-------------------------------------------------------

    def __add__(self, traj2):
        if (not isinstance(traj2, Trajectory)):
            raise self.TrajectoryError("Operation on Trajectory objects")

        if (self.dt != traj2.dt):
            raise self.TrajectoryError("Incompatible time increment")
        elif (not self.grid == traj2.grid):
            print("%d %d"%(self.grid.Ntrc,traj2.grid.Ntrc))
            print("%g %g"%(self.grid.dx,traj2.grid.dx))
            print("%g %g"%(self.grid.L,traj2.grid.L))
            print(self.grid.Ntrc==traj2.grid.Ntrc and 
                  self.grid.dx==traj2.grid.dx and 
                  self.grid.L==traj2.grid.L)
            raise self.TrajectoryError("Incompatible grids")

        trajSub=__data2Traj(self.grid, self.__data+traj2.__data,self.dt)
        return trajSub

    #-------------------------------------------------------

    def __mul__(self, scalar):
        trajMult=self.copy()
        trajMult.__data*=scalar
        return trajMult

    #-------------------------------------------------------
    #----| Public plotting methods |------------------------
    #-------------------------------------------------------

    class TrajectoryPlotError(Exception):
        pass

    #-------------------------------------------------------

    def waterfall(self, xlim=None, nbLines=50, title=None, 
                  offset=None, ampl=None, color='b', axe=None):
        """
            @TODO: ajout d'une echelle de grandeur
        """

        axe=self.__checkAxe(axe)
    
        if self.nDt<nbLines:
            nbLines=self.nDt

        freq=int((self.nDt+1.)/nbLines)
        if freq==1:
            nbLines=self.nDt

        if offset==None:
            offset=self.tInt/nbLines
        if ampl==None:
            ampl=8./(np.max(np.abs(self.__data)))*\
                    (self.tInt/nbLines)
        lines=[]
        for i in xrange(nbLines):
            j=i*freq
            curve=list(zip(self.grid.x, ampl*self.__data[j]
                           +self.time[j]-i*offset))
            lines.append(curve)

    
        col=collections.LineCollection(lines, offsets=(0.0,offset))
        col.set_color(color)
    
        axe.add_collection(col, autolim=True)
        axe.autoscale_view()
        axe.set_xlabel(r'$x$')
        axe.set_ylabel(r'$\tau$')
        axe.set_xlim(xlim)
        
        if title!=None:
           axe.set_title(title)

        return axe

    
    
    #-------------------------------------------------------

    def plotA(self, title=None, axe=None, **kwargs):

        axe=self.__checkAxe(axe)

        self.normA(ret=False)
        axe.plot(self.time, self.A, kwargs)
        
        if title!=None:
            axe.set_title(title)

        return axe

    #-------------------------------------------------------

    def plotA2(self, title=None, **kwargs):

        axe=self.__checkAxe(axe)

        self.normA2(ret=False)
        axe.plot(self.time, self.A2, **kwargs)
        
        if title!=None:
            axe.set_title(title)

        return axe

    #-------------------------------------------------------
    #----| Private plotting methods |-----------------------
    #-------------------------------------------------------

    def __checkAxe(self, axe):
        if not self.isIntegrated:
            raise self.TrajectoryPlotError(
                "Trajectory not integrated: nothing to plot")

        if axe==None:
            axe=plt.subplot(111)
        elif not (isinstance(axe,[matplotlib.axes.AxesSubplot,
                                  matplotlib.gridspec.GridSpec])):
            raise self.TrajectoryPlotError(
                "axe < matplotlib.axes.AxesSubplot | matplotlib.gridspec.GridSpec >")
        return axe
#--------------------------------------------------------------------
#####################################################################
#--------------------------------------------------------------------

if __name__=="__main__":

    grid1=SpectralGrid(100,100.)
    noTrj1=Trajectory(grid1)
    grid2=SpectralGrid(102,100.)
    noTrj2=Trajectory(grid2)
    print(noTrj1)
    print(noTrj1==noTrj2)
    #print(noTrj1+noTrj2)
