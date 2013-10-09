import numpy as np
from pseudoSpec1D import *

class Param(object):
    r"""
     Augmented KdV differential system parameters

        \partial_t A(x,t)= forcing(x) - \alpha(x) A 
                           - \beta(x)A \partial_x A 
                           - \gamma(x) \partial_x^3 A - \rho(x) A 

        Param(grid, alpha=None, beta=None, gamma=None, 
                    rho=None, forcing=None)

        Parameters can be specified either as constant or function
        In the later case, function must be of one variable
        (the space variable 'x')

        <TODO>  parameters defined as array  

    """
    class ParamError(Exception):
        pass


    def __init__(self, grid, forcing=None, alpha=None, beta=None,
                    gamma=None,  rho=None, nDt=0, dt=0.):

        if not (isinstance(grid, PeriodicGrid)):
            raise self.ParamError(
                "grid must be an instance of PeriodicGrid")
        self.grid=grid

        if (isinstance(forcing,Trajectory) or
                isinstance(alpha,Trajectory) or
                isinstance(beta,Trajectory) or
                isinstance(gamma,Trajectory) or
                isinstance(rho,Trajectory)):
            self.inputAsTrajectory=True
        else:
            self.inputAsTrajectory=False

        if (callable(forcing) or
                callable(alpha) or
                callable(beta) or
                callable(gamma) or
                callable(rho)):
            self.inputAsFunction=True
        else:
            self.inputAsFunction=False

        if self.inputAsFunction and dt==0. and nDt>0:
            raise self.ParamError(
            "if one parameter is <function(x,t)> and nDt>0 => dt>0")
        self.nDt=nDt
        self.dt=dt

        self.__initTraj(forcing, alpha, beta, gamma, rho)
        self.putDt(self.dt)
    #----------------------------------------------------------------
    #----| Public methods |------------------------------------------
    #----------------------------------------------------------------

    def max(self):
        val=0.
        for i in xrange(5):
            tryVal=np.max(self[i])
            if tryVal>val: val=tryVal
        return val

    #----------------------------------------------------------------

    def min(self):
        val=0.
        for i in xrange(5):
            tryVal=np.min(self[i])
            if tryVal<val: val=tryVal
        return val
    #----------------------------------------------------------------

    def putDt(self, dt):
        self.dt=dt
        self.tReal=self.dt*self.nDt
        self.alpha.incrmTReal(finished=True, tReal=self.tReal)
        self.beta.incrmTReal(finished=True, tReal=self.tReal)
        self.gamma.incrmTReal(finished=True, tReal=self.tReal)
        self.rho.incrmTReal(finished=True, tReal=self.tReal)
        self.forcing.incrmTReal(finished=True, 
                                    tReal=self.tReal)
        

    #----------------------------------------------------------------
    #----| Private methods |-----------------------------------------
    #----------------------------------------------------------------

    def __initTraj(self, forcing, alpha, beta, gamma, rho):
        params=(forcing, alpha, beta, gamma, rho)
        if self.inputAsTrajectory:
            for p in params:
                if isinstance(p, Trajectory):
                    if self.dt==0.:
                        self.dt=p.dt
                    else:
                        if p.dt<>self.dt:
                            raise ParamError(
                        "parameter trajectories must have the same dt")

                    if p.nDt>self.nDt:
                        self.nDt=p.nDt

        self.forcing=self.__setTraj(forcing)
        self.alpha=self.__setTraj(alpha)
        self.beta=self.__setTraj(beta)
        self.gamma=self.__setTraj(gamma)
        self.rho=self.__setTraj(rho)
        
        if self.nDt>0:
            self.isTimeDependant=True
        else:   
            self.isTimeDependant=False
        self.shape=(5, self.nDt ,self.grid.N)
    #----------------------------------------------------------------

    def __setTraj(self, param):

        if isinstance(param, Trajectory):
            if param.nDt<self.nDt:
                data=params[i].getData()
                newdata=np.empty(shape=(self.nDt+1, self.grid.N))
                newdata[0:params[i].nDt+1,:]=data
                newdata[params[i].nDt+1:, :]=data[params[i].nDt+1, :]
                traj=Trajectory(self.grid)
                traj.zeros(self.nDt)
                traj.putData(newdata)
                 
            if param.nDt==self.nDt:
                traj=param
    
        elif param==None:
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt)
            traj.putData(np.zeros(shape=(self.nDt+1, self.grid.N)))

        elif isinstance(param, (int, float)):
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt)
            traj.putData(np.ones(shape=(self.nDt+1, self.grid.N))\
                            *param)

        elif callable(param):
            data=self.__matricize(param)
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt)
            traj.putData(data)

        else:
            raise self.ParamError(
                    "<None|float|function(x,t)|Trajectory>")
        
        return traj
            
    #----------------------------------------------------------------

    def __getitem__(self,i):
        if i==0:
            return self.forcing
        elif i==1:
            return self.alpha
        elif i==2:
            return self.beta
        elif i==3:
            return self.gamma
        elif i==4:
            return self.rho
        else:
            raise IndexError('[0:forcing(x,t), 1:alpha(x,t),'+
                            '2:beta(x,t), 3:gamma(x,t), 4:rho(x,t)]')

    
    #----------------------------------------------------------------

    def __matricize(self, funcXT):
        time=np.linspace(0., self.nDt*self.dt, self.nDt+1)
        matrix=np.empty(shape=(self.nDt+1, self.grid.N))
        for i in xrange(self.grid.N):
            for j in xrange(self.nDt+1):
                matrix[j][i]=funcXT(self.grid.x[i], time[j])

        return matrix
        

#==========================================================
#----------------------------------------------------------
#==========================================================
if __name__=='__main__':
    import matplotlib.pyplot as plt
    
    grid=PeriodicGrid(142, 5.)
    
    def func1(x,t):
        return  -1.
    
    def func2(x,t):
        return 0.
    
    def func3(x,t):
        return np.sin(x)

    def funcTimeDependant(x, t):
        return np.sin(x-t)*np.cos(t)
    nDtParam=100
    dt=0.1

    p=Param(grid, funcTimeDependant, func2, func1, func1, func3,
                nDt=nDtParam, dt=dt)
    if p.isTimeDependant:
        p[0].waterfall()
    else:
        plt.plot(grid.x, p[0].getData()[0,:])
    
    plt.show() 


