import numpy as np
import random as rnd
import copy

import matplotlib.pyplot as plt
from matplotlib import collections, axes, gridspec
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec, SubplotSpec

from periodicGrid import PeriodicGrid, Grid

class Trajectory(object):
    """
    Trajectory class
    for periodic 1+1D partial differential system

        Trajectory(grid)

        grid    :   <Grid>
    """
    
    class TrajectoryError(Exception):
        pass

    #------------------------------------------------------
    #----| Init |------------------------------------------
    #------------------------------------------------------


    def __init__(self, grid, label=None):
        """
        Empty Trajectory instance constructor
        """
        # Grid
        if not (isinstance(grid, Grid)):
            raise self.TrajectoryError(
                  "grid <Grid>")
        self.grid=grid

        # Time attributes
        self.tReal=0.
        self.nDt=0
        self.dt=0.
        
        # Status Attributes
        self.isIntegrated=False
        self.isInitialised=False
        self.isConcatenated=False
        self.isTrimmed=False

        # Data attributes
        self.__data=None
        self.ic=None
        self.A=None
        self.A2=None
        self.time=None

        self.label=label

    #------------------------------------------------------
    #----| general public methods |------------------------
    #------------------------------------------------------

    def setLabel(self, label):
        self.label=label
    
    #-------------------------------------------------------

    def copy(self):
        return copy.deepcopy(self)

        
    #-------------------------------------------------------
    
    def initialize(self, ic, nDt, dt):
        """
        Trajectory initialization
            (memory allocation and IC association)

            Trajectory.initialize(ic, nDt, dt)

            ic  :   initial condition <numpy.ndarray>
            nDt :   time steps <int>
            dt  :   time step increment <float>
        """
        if not (isinstance(ic, np.ndarray)):
            raise self.TrajectoryError(
                  "ic <numpy.ndarray>")
        if ic.ndim <> 1 or ic.size <> self.grid.N:
            raise self.TrajectoryError("ic.shape<>(grid.N,)")
        self.__allocate(nDt)
        self.dt=dt
        self.ic=ic
        self.isInitialised=True

    #-------------------------------------------------------
    
    def zeros(self, nDt, dt=0.):
        """
        Trajectory initialization
            (memory allocation)

            Trajectory.zeros(nDt)

            nDt :   time steps <int>
        """
        self.__allocate(nDt)
        self.dt=dt
        self.ic=np.zeros(self.grid.N)
        self.isInitialised=True
    
    #-------------------------------------------------------
 
    def null(self, nDt, dt):
        """ 
        Nul trajectory

        """

        self.initialize(np.zeros(self.grid.N), nDt, dt)
        self.putData(np.zeros(shape=(self.nDt+1, self.grid.N)))
        self.incrmTReal(finished=True, tReal=self.nDt*self.dt)
    
    #-------------------------------------------------------
    
    def incrmTReal(self, finished=False, tReal=None, t0=0.):
        if tReal<>None:
            self.tReal=tReal
        if finished:
            self.isIntegrated=True
            self.t0=t0
            self.tf=self.tReal+self.t0
            self.time=np.linspace(self.t0,self.tf, self.nDt+1)
            self.final=self.__data[self.nDt]

        else:
            self.tReal+=self.dt

    #-------------------------------------------------------

    def norm2(self, ret=True, metric=None):
        """
        Euclidian Square norm evolution
        """
        if not self.isIntegrated:
            raise self.TrajectoryError("Trajectory not integrated")
        self.A2=np.zeros(self.nDt+1)
        for i in xrange(self.nDt+1):
            self.A2[i]=self.__SquareNorm(self.__data[i], metric=metric)
        if ret:
            return self.A2

    #-------------------------------------------------------

    def norm(self, ret=True, metric=None):
        """
        Euclidian Norm evolution
        """
        if not self.isIntegrated:
            raise self.TrajectoryError("Trajectory not integrated")
        self.A=np.zeros(self.nDt+1)
        if (self.A2==None):
            self.norm2(ret=False, metric=metric)
        self.A=np.sqrt(self.A2)
        if ret:
            return self.A

    #-------------------------------------------------------

    def getData(self):
        """
        Trajectory values access method (output)
        """
        return self.__data

    #-------------------------------------------------------

    def putData(self, data, tReal=None, t0=0.):
        """
        Trajectory values access method (input)
        """
        if not self.isInitialised :
            raise self.TrajectoryError("Trajectory not initialised")
        if not isinstance(data, np.ndarray):
            raise self.TrajectoryError("data <numpy.ndarray>")
        if (data.shape<>(self.nDt+1, self.grid.N)):
            raise self.TrajectoryError("Incompatible data affectation")
        self.__data=data
        self.shape=data.shape
        self.ic=self[0]

        if tReal<>None:
            self.incrmTReal(finished=True, tReal=tReal, t0=t0)

    #-------------------------------------------------------

    def whereTimeIdx(self, time):
        return np.where(self.time>=time)[0].min()

    def whereTime(self, time):
        """
        Return the instantaneous field at specified time

            Trajectory.whereTime(time)

            time    :   desired time <float>

            <!> will find the nearest egal or superior time
                in the discretized time serie
        """
        return self.__data[self.whereTimeIdx(time)]


    #-------------------------------------------------------
    
    def max(self):
        return self.getData().max()


    def min(self):
        return self.getData().min()
    
    def abs(self):
        return np.abs(self.getData())

    #-------------------------------------------------------

    def concatenate(self, traj):
        """
        <TODO>
        Concatenate two trajectories which are consecutive

        Trajectories must share spacial and temporal caracteristics
        """
        if not isinstance(traj, Trajectory):
            raise self.TrajectoryError("traj <Trajectory>")
        if traj.dt<>self.dt:
            raise self.TrajectoryError(
                "incompatible time increments: %f, %f"%(self.dt, traj.dt))
        if not traj.grid==self.grid:
            raise self.TrajectoryError("incompatible grid")
            
        trajComp=Trajectory(self.grid)
        trajComp.initialize(self.ic, self.nDt+traj.nDt,self.dt)
        data=np.empty(shape=(trajComp.nDt+1, self.grid.N))
        for t in xrange(self.nDt+1):
            data[t]=self.__data[t]
        trajData=traj.getData()
        for i in xrange(1,traj.nDt+1):
            t=i+self.nDt
            data[t]=trajData[i]

        trajComp.putData(data)
        trajComp.incrmTReal(finished=True, tReal=trajComp.nDt*trajComp.dt)
        trajComp.isConcatenated=True
        trajComp.isTrimmed=self.isTrimmed
        return trajComp
        
    #-------------------------------------------------------

    def trim(self, freq):


        if not isinstance(freq, int):
            raise self.TrajectoryError("freq <int>")
        if freq<1 or freq>self.nDt:
            raise self.TrajectoryError("1 =< freq < self.nDt")

        nDtTrim=(self.nDt+1)/freq-1
        trajTrim=Trajectory(self.grid)
        trajTrim.initialize(self.ic, nDtTrim+1,self.dt*freq)
        trajTrim[0]=self.ic
        for i in xrange(trajTrim.nDt):
            trajTrim[i+1]=self.__data[i*freq]
        trajTrim.incrmTReal(finished=True, tReal=trajTrim.nDt*trajTrim.dt,
                            t0=self.t0)
        trajTrim.isTrimmed=True
        trajTrim.isConcatenated=self.isConcatenated
        return trajTrim
    
    #-------------------------------------------------------
    
    def cut(self, t0=None, tf=None):


        if t0==None:
            idx0=0
            t0=self.t0
        elif t0<self.tReal:
            t0=float(t0)
            idx0=self.whereTimeIdx(t0)
        else: raise ValueError()
        if tf==None or tf==self.tReal:
            idxF=self.nDt
        elif tf<self.tReal:
            tf=float(tf)
            idxF=self.whereTimeIdx(tf)
            if idxF>self.nDt: idxF=self.nDt
        else: raise ValueError()

        cutNDt=idxF-idx0-1
        cutTraj=Trajectory(self.grid)
        cutTraj.initialize(self[idx0], cutNDt+1, self.dt)
        for i in xrange(cutTraj.nDt+1):
            cutTraj[i]=self.__data[idx0+i]
        cutTraj.incrmTReal(finished=True, tReal=(cutNDt+1)*self.dt, t0=t0)
        return cutTraj

    #------------------------------------------------------

    def degrad(self, mu, sigma, seed=None):
        '''
        Gaussian noise trajectory degradation
    
        degrad(u,mu,sigma,seed=...)
    
        mu      :  noise mean (gaussian mean)
        sigma   :  noise variance
        '''
        if not isinstance(self, Trajectory):
            raise TypeError("self <Trajectory>")
        rnd.seed(seed)
        ic_degrad=self.ic.copy()
        for i in xrange(self.grid.N):
            ic_degrad[i]=self.ic[i]+rnd.gauss(mu, sigma)
        
        degrad=Trajectory(self.grid)
        degrad.initialize(ic_degrad, self.nDt, self.dt)
        degrad[0]=ic_degrad
        for i in xrange(1,self.nDt+1):
            for j in xrange(self.grid.N):
                degrad[i][j]=self[i][j]+rnd.gauss(mu, sigma)
        degrad.incrmTReal(finished=True, tReal=self.tReal, t0=self.t0)
        return degrad

    #------------------------------------------------------

    def gradient(self, n=1):

        delXTraj=Trajectory(self.grid)
        delXTraj.initialize(self.grid.zeros(), self.nDt, self.dt)
        data=np.empty(shape=self.shape)
        for t in xrange(self.nDt+1):
            data[t]=self[t]
            for i in xrange(n):
                data[t]=self.grid.gradient(data[t])
        delXTraj.putData(data)
        delXTraj.incrmTReal(finished=True, tReal=self.tReal, t0=self.t0)
        return delXTraj

    #------------------------------------------------------

    def delT(self):

        if not self.isIntegrated:
            raise self.TrajectoryError("Trajectory not integrated") 
        delTTraj=Trajectory(self.grid)
        if self.nDt==0:
            delTTraj.null(self.nDt, self.dt)
        else:
            delTTraj.initialize(self.grid.zeros(), self.nDt, self.dt)
            data=np.empty(shape=self.shape)
            data=np.gradient(self.getData(), self.dt)[0]
            delTTraj.putData(data)
            delTTraj.incrmTReal(finished=True, tReal=self.tReal, t0=self.t0)
        return delTTraj
    
    #------------------------------------------------------

    def pointMult(self, traj, filter=False):
        # add coherency check (nDt, t0, ...)
        prodTraj=Trajectory(self.grid)
        prodTraj.initialize(self.grid.zeros(), self.nDt, self.dt)
        data=np.empty(shape=self.shape)
        for t in xrange(self.nDt+1):
            data[t]=self[t]*traj[t]
        if filter:
            pass
        prodTraj.putData(data)
        prodTraj.incrmTReal(finished=True, tReal=self.tReal)
        return prodTraj


    #------------------------------------------------------
    #----| Private methods |-------------------------------
    #------------------------------------------------------

    def __SquareNorm(self, z, metric=None):
        return self.grid.squareNorm(z, metric=metric)

    #-------------------------------------------------------

    def __allocate(self, nDt):
        self.nDt=nDt
        self.__data=np.zeros(shape=(self.nDt+1,self.grid.N))
        self.shape=self.__data.shape

    #-------------------------------------------------------

    def compatibility(self, traj2):
        if (not isinstance(traj2, Trajectory)):
            raise self.TrajectoryError("Operation on Trajectory objects")

        if (self.dt != traj2.dt):
            #raise self.TrajectoryError("Incompatible time increment")
            return False
        if (self.nDt != traj2.nDt):
            #raise self.TrajectoryError("Incompatible time step number")
            return False
        if (self.t0 != traj2.t0):
            return False
        if (not self.grid == traj2.grid):
            #print("%d %d"%(self.grid.N,traj2.grid.N))
            #print("%g %g"%(self.grid.dx,traj2.grid.dx))
            #print("%g %g"%(self.grid.L,traj2.grid.L))
            #print(self.grid.N==traj2.grid.N and 
            #      self.grid.dx==traj2.grid.dx and 
            #      self.grid.L==traj2.grid.L)
            #raise self.TrajectoryError("Incompatible grids")
            return False
        return True

    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="====| Trajectory |=================================\n"
        output+=self.grid.__str__()
        output+="\n| nDt=%d"%self.nDt
        output+="\n| dt=%-23.15E"%self.dt
        output+="\n| tReal=%-23.15E"%self.tReal
        output+="\n| t0=%-23.15E"%self.t0
        if self.isConcatenated:
            output+="\n is concatenated"
        if self.isTrimmed:
            output+="\n is trimmed"
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

        if (not self.compatibility(traj2)):
           raise self.TrajectoryError("Incompatible trajectories")

        trajSub=self.copy()
        trajSub.putData(self.__data-traj2.__data, tReal=self.tReal,
                        t0=self.t0)
        return trajSub

    #-------------------------------------------------------

    def __add__(self, traj2):
        if (not isinstance(traj2, Trajectory)):
            raise self.TrajectoryError("Operation on Trajectory objects")

        if (not self.compatibility(traj2)):
           raise self.TrajectoryError("Incompatible trajectories")

        trajAdd=self.copy()
        trajAdd.putData(self.__data+traj2.__data, tReal=self.tReal,
                            t0=self.t0)
        return trajAdd

    #-------------------------------------------------------

    def __mul__(self, scalar):
        trajMult=self.copy()
        trajMult.__data*=scalar
        trajMult.ic*=scalar
        return trajMult

    #-------------------------------------------------------
    #----| Public plotting methods |------------------------
    #-------------------------------------------------------

    class TrajectoryPlotError(Exception):
        pass

    #-------------------------------------------------------
    
    def fftTraj(self):
        """
        Return the SpectralTrajectory associated
        """
        nDemi=int(self.grid.N-1)/2
        data=np.zeros(shape=(self.nDt+1, nDemi))
        for i in xrange(len(self.time)):
            data[i]=np.abs(np.fft.fft(self[i])[0:nDemi])

        k=Grid(nDemi,nDemi, centered=False)
        fftTraj=SpectralTrajectory(k, Ntrc=self.grid.Ntrc)
        fftTraj.initialize(data[0],self.nDt, self.dt)
        fftTraj.putData(data)
        fftTraj.incrmTReal(finished=True, tReal=self.tReal, t0=self.t0)
        return fftTraj
    
    #-------------------------------------------------------
    
    def waterfall(self, xlim=None, ylim=None, nbLines=50, title=None, 
                  offset=None, ampl=None, color='b', axe=None):
        """
        Make a normalized waterfall plot of the Trajectory

            Trajectory.waterfall(xlim=None, nbLines=50, title=None, 
                  offset=None, ampl=None, color='b', axe=None)

                xlim    :   x axis limits <tuple>
                nbLines :   number of lines in the plot <int>
                title   :   title <str>
                color   :   color <str>
                axe     :   subplot object <Axes | GridSpec>

            @TODO: a scale
        """

        if isinstance(self, SpectralTrajectory):
            xlabel="$k$"
        else:
            xlabel="$x$"

        axe=self._checkAxe(axe)
    
        if not self.isIntegrated or self.nDt==0 :
            axe.plot(self.grid.x, self[0], color)
            axe.set_xlim(xlim)
            axe.set_xlabel(xlabel)
            axe.set_ylabel(r'$A$')

            if title!=None:
                axe.set_title(title)
            return axe

        if ylim==None:
            tReal=self.tReal
            nDt=self.nDt
            data=self.__data
            time=self.time
        elif isinstance(ylim, tuple):
            tReal=ylim[1]-ylim[0]
            lIdx=self.whereTimeIdx(ylim[0])
            hIdx=self.whereTimeIdx(ylim[1])
            if len(ylim)<>2 or hIdx<lIdx:
                raise self.TrajectoryPlotError("ylim=(tMin,tMax)")
            nDt=hIdx-lIdx
            data=self.__data[lIdx:hIdx]
            time=self.time[lIdx:hIdx]
        else:
                raise self.TrajectoryPlotError("ylim=(tMin,tMax)")

        if nDt<nbLines:
            nbLines=nDt

        freq=int((nDt+1.)/nbLines)
        if freq==1:
            nbLines=nDt

        if offset==None:
            offset=tReal/nbLines
        if ampl==None:
            denom=(np.max(np.abs(data)))
            if denom==0.: denom=1.
            ampl=8./denom*\
                    (tReal/nbLines)
        lines=[]
        for i in xrange(nbLines):
            j=i*freq
            curve=list(zip(self.grid.x, ampl*data[j]
                           +time[j]-i*offset))
            lines.append(curve)

    
        col=collections.LineCollection(lines, offsets=(0.0,offset))
        col.set_color(color)
    
        axe.add_collection(col, autolim=True)
        axe.autoscale_view()
        axe.set_xlabel(xlabel)
        axe.set_ylabel(r'$t$')
        axe.set_xlim(xlim)
        
        if title!=None:
           axe.set_title(title)

        return axe

    
    
    #-------------------------------------------------------

    def plotA(self, title=None, axe=None,  metric=None, **kwargs):
        """
        Amplitude evolution plot

                Trajectory.plotA(title=None, axe=None)

                title       :   title <str>
                axe         :   subplot object <Axes | GridSpec>
        """
        axe=self._checkAxe(axe)

        self.norm(ret=False, metric=metric)
        axe.plot(self.time, self.A, **kwargs)
        axe.set_xlabel("$t$")
        axe.set_ylabel(r"$\Vert A\Vert$")
        
        if title!=None:
            axe.set_title(title)

        return axe

    #-------------------------------------------------------

    def plotA2(self, title=None, axe=None, metric=None,  **kwargs):
        """
        Square Amplitude evolution plot

                Trajectory.plotA(title=None, axe=None, **kwargs )

                title       :   title <str>
                axe         :   subplot object <Axes | GridSpec>
        """

        axe=self._checkAxe(axe)

        self.norm2(ret=False, metric=metric)
        axe.plot(self.time, self.A2, **kwargs)
        axe.set_xlabel("$t$")
        axe.set_ylabel(r"$\Vert A\Vert^2$")
 
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
            raise self.TrajectoryPlotError(
                "axe < matplotlib.axes.Axes | matplotlib.gridspec.GridSpec >")
        return axe

#--------------------------------------------------------------------
#====================================================================  
#--------------------------------------------------------------------
class SpectralTrajectory(Trajectory):
    """
    Spectral Trajectory class
    for periodic 1+1D partial differential system

        Trajectory(grid, Ntrc=None)

        grid    :   <Grid>
        Ntrc    :   truncature <int>
    """

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

    def waterfall(self, xlim=None, ylim=None, nbLines=50, title=None, 
                  offset=None, ampl=None, color='b', axe=None):
        """
        Make a normalized waterfall plot of the Spectral Trajectory

            Trajectory.waterfall(xlim=None, nbLines=50, title=None, 
                  offset=None, ampl=None, color='b', axe=None)

                xlim    :   x axis limits <tuple>
                nbLines :   number of lines in the plot <int>
                title   :   title <str>
                color   :   color <str>
                axe     :   subplot object <Axes | GridSpec>

            @TODO: a scale
        """
        axe=self._checkAxe(axe)
        if xlim==None:
            xlim=(0,self.Ntrc)
        super(SpectralTrajectory, self).waterfall(xlim, ylim, nbLines,
                                                    title, offset, ampl,
                                                    color, axe)
        if xlim[1]>self.Ntrc: 
            axe.axvline(x=self.Ntrc, color='k', linewidth=2)
    
        return axe


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
