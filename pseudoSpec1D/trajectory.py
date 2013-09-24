import numpy as np
import copy

import matplotlib.pyplot as plt
from matplotlib import collections, axes, gridspec
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec

from periodicGrid import PeriodicGrid, Grid

class Trajectory(object):
    """
    Trajectory class
    for periodic 1+1D partial differential system

    """
    
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
        if not (isinstance(grid, Grid)):
            raise self.TrajectoryError(
                  "grid <Grid>")
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
            self.final=self.__data[self.nDt]

        else:
            self.tReal+=self.dt

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

    def whereTimeIdx(self, time):
        return np.where(self.time>=time)[0].min()

    def whereTime(self, time):
        return self.__data[self.whereTimeIdx(time)]



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
            print("%d %d"%(self.grid.N,traj2.grid.N))
            print("%g %g"%(self.grid.dx,traj2.grid.dx))
            print("%g %g"%(self.grid.L,traj2.grid.L))
            print(self.grid.N==traj2.grid.N and 
                  self.grid.dx==traj2.grid.dx and 
                  self.grid.L==traj2.grid.L)
            raise self.TrajectoryError("Incompatible grids")

        trajSub=putData(self.__data-traj2.__data,self.dt)
        return trajSub

    #-------------------------------------------------------

    def __add__(self, traj2):
        if (not isinstance(traj2, Trajectory)):
            raise self.TrajectoryError("Operation on Trajectory objects")

        if (self.dt != traj2.dt):
            raise self.TrajectoryError("Incompatible time increment")
        elif (not self.grid == traj2.grid):
            print("%d %d"%(self.grid.N,traj2.grid.N))
            print("%g %g"%(self.grid.dx,traj2.grid.dx))
            print("%g %g"%(self.grid.L,traj2.grid.L))
            print(self.grid.N==traj2.grid.N and 
                  self.grid.dx==traj2.grid.dx and 
                  self.grid.L==traj2.grid.L)
            raise self.TrajectoryError("Incompatible grids")

        trajSub=putData(self.__data+traj2.__data,self.dt)
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
    
    def fftTraj(self):
        nDemi=int(self.grid.N-1)/2
        data=np.zeros(shape=(self.nDt+1, nDemi))
        for i in xrange(len(self.time)):
            data[i]=np.abs(np.fft.fft(self[i])[0:nDemi])

        k=Grid(nDemi,nDemi, centered=False)
        fftTraj=SpectralTrajectory(k, Ntrc=self.grid.Ntrc)
        fftTraj.initialize(data[0],self.nDt, self.dt)
        fftTraj.putData(data)
        fftTraj.incrmTReal(finished=True, tReal=self.tReal)
        return fftTraj
    
    #-------------------------------------------------------
    
    def waterfall(self, xlim=None, nbLines=50, title=None, 
                  offset=None, ampl=None, color='b', axe=None):
        """
            @TODO: ajout d'une echelle de grandeur
        """

        axe=self._checkAxe(axe)
    
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

        axe=self._checkAxe(axe)

        self.normA(ret=False)
        axe.plot(self.time, self.A, kwargs)
        
        if title!=None:
            axe.set_title(title)

        return axe

    #-------------------------------------------------------

    def plotA2(self, title=None, **kwargs):

        axe=self._checkAxe(axe)

        self.normA2(ret=False)
        axe.plot(self.time, self.A2, **kwargs)
        
        if title!=None:
            axe.set_title(title)

        return axe

    #-------------------------------------------------------
    #----| Private plotting methods |-----------------------
    #-------------------------------------------------------

    def _checkAxe(self, axe):
        if not self.isIntegrated:
            raise self.TrajectoryPlotError(
                "Trajectory not integrated: nothing to plot")

        if axe==None:
            axe=plt.subplot(111)
        elif not (isinstance(axe,(Axes, GridSpec))):
            raise self.TrajectoryPlotError(
                "axe < matplotlib.axes.Axes | matplotlib.gridspec.GridSpec >")
        return axe

#--------------------------------------------------------------------
#====================================================================  
#--------------------------------------------------------------------
class SpectralTrajectory(Trajectory):

    class SpectralTrajectoryError(Exception):
        pass

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------


    def __init__(self, grid, Ntrc=None):
        super(SpectralTrajectory, self).__init__(grid)
        if Ntrc==None:
            self.Ntrc=grid.N
        else:
            self.Ntrc=Ntrc

    def waterfall(self, xlim=None, nbLines=50, title=None, 
                  offset=None, ampl=None, color='b', axe=None):
        axe=self._checkAxe(axe)
        super(SpectralTrajectory, self).waterfall(xlim, nbLines, title,
                                                    offset, ampl, color,
                                                    axe)
        axe.axvline(x=self.Ntrc, color='k', linewidth=2)
    


#--------------------------------------------------------------------
#####################################################################
#--------------------------------------------------------------------

if __name__=="__main__":

    grid1=PeriodicGrid(100,100.)
    noTrj1=Trajectory(grid1)
    grid2=PeriodicGrid(102,100.)
    noTrj2=Trajectory(grid2)
    print(noTrj1)
    print(noTrj1==noTrj2)
    #print(noTrj1+noTrj2)
