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
                    gamma=None,  rho=None):

        if not (isinstance(grid, PeriodicGrid)):
            raise self.ParamError("grid must be an instance of PeriodicGrid")
        self.grid=grid
        
        if (isinstance(forcing, Trajectory) 
            or isinstance(alpha, Trajectory) 
            or isinstance(beta, Trajectory)
            or isinstance(gamma, Trajectory)
            or isinstance(rho, Trajectory)):
            self.__initTraj(forcing, alpha, beta, gamma, rho)

            self.isTimeDependant=True
            self.shape=(5, self.nDt ,grid.N)

        else:
            self.forcing=self.__setFunc(forcing)
            self.alpha=self.__setFunc(alpha)
            self.beta=self.__setFunc(beta)
            self.gamma=self.__setFunc(gamma)
            self.rho=self.__setFunc(rho)
            
            self.isTimeDependant=False
            self.shape=(5,grid.N)

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
        if not self.isTimeDependant:
            raise self.ParamError(
                    "not a time dependant parametrisation")
        else:
            self.dt=dt
            self.alpha.incrmTReal(finished=True, tReal=self.dt*self.nDt)
            self.beta.incrmTReal(finished=True, tReal=self.dt*self.nDt)
            self.gamma.incrmTReal(finished=True, tReal=self.dt*self.nDt)
            self.rho.incrmTReal(finished=True, tReal=self.dt*self.nDt)
            self.forcing.incrmTReal(finished=True, tReal=self.dt*self.nDt)

    #----------------------------------------------------------------
    #----| Private methods |-----------------------------------------
    #----------------------------------------------------------------

    def __initTraj(self, forcing, alpha, beta, gamma, rho):
        params=(forcing, alpha, beta, gamma, rho)
        nDtParam=0
        issetDt=False
        for p in params:
            if isinstance(p, Trajectory):
                if not issetDt:
                    self.dt=p.dt
                    issetDt=True
                else:
                    if p.dt<>self.dt:
                        raise ParamError(
                        "parameter trajectories must have the same dt")

                if p.nDt>nDtParam:
                    nDtParam=p.nDt
            self.nDt=nDtParam

        self.forcing=self.__setTraj(forcing)
        self.alpha=self.__setTraj(alpha)
        self.beta=self.__setTraj(beta)
        self.gamma=self.__setTraj(gamma)
        self.rho=self.__setTraj(rho)
            

    #----------------------------------------------------------------

    def __setTraj(self, param):

        if isinstance(param, Trajectory):
            if param.nDt<self.nDt:
                data=params[i].getData()
                newdata=np.empty(shape=(self.nDt, self.grid.N))
                newdata[0:params[i].nDt,:]=data
                newdata[params[i].nDt:, :]=data[params[i].nDt, :]
                traj=Trajectory(self.grid)
                traj.zeros(self.nDt-1)
                traj.putData(newdata)
                 
            if param.nDt==self.nDt:
                traj=param
    
        elif param==None:
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt-1)
            traj.putData(np.zeros(shape=(self.nDt, self.grid.N)))
        elif isinstance(param, (int, float)):
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt-1)
            traj.putData(np.ones(shape=(self.nDt, self.grid.N))\
                            *param)

        elif callable(param):
            constData=np.vectorize(param)(self.grid.x)
            data=np.empty(shape=(self.nDt, self.grid.N))
            data[:]=constData
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt-1)
            traj.putData(data)

        else:
            raise self.ParamError("<None|float|Trajectory>")
        
        return traj
            
    #----------------------------------------------------------------

    def __setFunc(self, param):
        def wrapper(c):
            def funcConst(x):
                return float(c)
            return funcConst
        
        if param==None:
            return wrapper(0.)
        elif isinstance(param, (int, float)):
            return wrapper(param)
        elif callable(param):
            return param
        else:
            raise self.ParamError("<None|float|function>")
            


    #----------------------------------------------------------------

    def __getitem__(self,i):
        if self.isTimeDependant:
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
        else:
            if i==0:
                return np.vectorize(self.forcing)(self.grid.x)
            elif i==1:
                return np.vectorize(self.alpha)(self.grid.x)
            elif i==2:
                return np.vectorize(self.beta)(self.grid.x)
            elif i==3:
                return np.vectorize(self.gamma)(self.grid.x)
            elif i==4:
                return np.vectorize(self.rho)(self.grid.x)
            else:
                raise IndexError('[0:forcing(x), 1:alpha(x),'+
                                '2:beta(x), 3:gamma(x), 4:rho(x)]')
        


#==========================================================
#----------------------------------------------------------
#==========================================================
if __name__=='__main__':
    import matplotlib.pyplot as plt
    
    grid=PeriodicGrid(142, 5.)
    
    def func1(x):
        return  -1.
    
    def func2(x):
        return 0.
    
    def func3(x):
        return np.sin(x)

    def funcTimeDependant(x, t):
        return np.sin(x-t)*np.cos(t)
    nDtParam=100
    dt=0.1
    data=np.zeros(shape=(nDtParam+1, grid.N))
    for n in xrange(nDtParam+1):
        data[n,:]=funcTimeDependant(grid.x, n*dt)
    traj=Trajectory(grid)
    traj.initialize(data[0,:], nDtParam, dt)
    traj.putData(data)

    p=Param(grid, traj, func2, func1, func1, func3)
    if p.isTimeDependant:
        p[0].waterfall()
    else:
        plt.plot(grid.x, p[4])
    
    plt.show() 


