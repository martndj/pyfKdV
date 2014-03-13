import numpy as np
from pseudoSpec1D import *

class Param(object):
    """
     Augmented KdV differential system parameters

        \partial_t A(x,t)= forcing(x) - \alpha(x) A 
                           - \beta(x)A \partial_x A 
                           - \gamma(x) \partial_x^3 A - \rho(x) A 

        Param(grid, alpha=0., beta=0., gamma=0., rho=0., forcing=0.,
                nDt=0, dt=0., tInt=0.)

        Parameters can be specified either as constant, function
        or Trajectory object.
        None will
        In the function case, function must be of two variables
        (the space variable 'x' and time variable 'y')

            def func(x,t):
                ...
                return <some expression of x and t>

            <@> although both must be in the function declaration
                it is not necessary for anyone of them to actually
                be used in the return expression.

    """
    class ParamError(Exception):
        pass


    def __init__(self, grid, forcing=0., alpha=0., beta=1.,
                    gamma=-1.,  rho=0., nDt=0, dt=0., tInt=0., t0=0.):

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


        if nDt>0. or tInt>0.:
            self.isTimeDependant=True
        else:
            self.isTimeDependant=False
            
        if self.isTimeDependant:
            if self.inputAsFunction and dt==0.: 
                raise self.ParamError("if one parameter is <function(x,t)> and (nDt>0 or tInt>0) => dt>0")
            if nDt==0.:
                self.nDt=int(tInt/dt)
            else:
                self.nDt=nDt
        else:
            self.nDt=nDt
        self.dt=dt
        self.t0=t0

        self.__initTraj(forcing, alpha, beta, gamma, rho, t0=t0)
        self.__finalize()
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

    def whereTime(self, time):
        return Param(self.grid, 
                forcing=self.forcing.whereTime(time),
                alpha=self.alpha.whereTime(time),
                beta=self.beta.whereTime(time),
                gamma=self.gamma.whereTime(time),
                rho=self.rho.whereTime(time), t0=time)

    #----------------------------------------------------------------

    def cut(self, t0=None, tf=None):
   
        return Param(self.grid, 
                forcing=self.forcing.cut(t0=t0,tf=tf),
                alpha=self.alpha.cut(t0=t0, tf=tf),
                beta=self.beta.cut(t0=t0, tf=tf),
                gamma=self.gamma.cut(t0=t0, tf=tf),
                rho=self.rho.cut(t0=t0, tf=tf),
                dt=self.dt, nDt=self.nDt)

    #----------------------------------------------------------------

    def __incrmTReal(self, dt, t0=0.):
        self.tReal=self.dt*self.nDt
        self.alpha.incrmTReal(finished=True, tReal=self.tReal, t0=t0)
        self.beta.incrmTReal(finished=True, tReal=self.tReal, t0=t0)
        self.gamma.incrmTReal(finished=True, tReal=self.tReal, t0=t0)
        self.rho.incrmTReal(finished=True, tReal=self.tReal, t0=t0)
        self.forcing.incrmTReal(finished=True, 
                                    tReal=self.tReal, t0=t0)
        
    #----------------------------------------------------------------

    def __finalize(self):
        self.time=self.forcing.time
        self.tf=self.time.max()
        if self.isTimeDependant:
            self.final=Param(self.grid, 
                            forcing=self.forcing.final,
                            alpha=self.alpha.final,
                            beta=self.beta.final,
                            gamma=self.gamma.final,
                            rho=self.rho.final, t0=self.tf)
               
        self.fluxR=self.alpha.gradient()-self.rho
        self.fluxR.setLabel(r'$r$')


    #----------------------------------------------------------------
    #----| Private methods |-----------------------------------------
    #----------------------------------------------------------------

    def __initTraj(self, forcing, alpha, beta, gamma, rho, t0=0.):
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
                        self.isTimeDependant=True

        self.forcing=self._setTraj(forcing, t0=t0, label=r'$f$')
        self.alpha=self._setTraj(alpha, t0=t0, label=r'$\alpha$')
        self.beta=self._setTraj(beta, t0=t0, label=r'$\beta$')
        self.gamma=self._setTraj(gamma, t0=t0, label=r'$\gamma$')
        self.rho=self._setTraj(rho, t0=t0, label=r'$\rho$')
        
        self.shape=(5, self.nDt ,self.grid.N)
        self.__incrmTReal(self.dt, t0=t0)
    #----------------------------------------------------------------

    def _setTraj(self, param, t0=0., label=None):

        if isinstance(param, Trajectory):
            if param.nDt<self.nDt:
                data=param.getData()
                newdata=np.empty(shape=(self.nDt+1, self.grid.N))
                newdata[0:param.nDt+1,:]=data
                newdata[param.nDt+1:, :]=data[param.nDt]
                traj=Trajectory(self.grid)
                traj.zeros(self.nDt, self.dt)
                traj.putData(newdata)
                 
            if param.nDt==self.nDt:
                traj=param
    
        elif isinstance(param, (int, float)):
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt, self.dt)
            traj.putData(np.ones(shape=(self.nDt+1, self.grid.N))\
                            *param)

        elif callable(param):
            data=self._matricize(param)
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt, self.dt)
            traj.putData(data)

        elif isinstance(param, np.ndarray):
            if param.shape[-1] <> self.grid.N:
                raise self.ParamError("param.shape[-1]=self.grid.N")
                
            if param.ndim==1:
                data=np.outer(np.ones(self.nDt+1), param)
            elif param.ndim==2:
                if param.shape[0] <> self.nDt+1:
                    raise self.ParamError("param.shape[0]=self.nDt+1")
                data=param
            else:   
                raise self.ParamError("param.ndim in [1,2]")
            traj=Trajectory(self.grid)
            traj.zeros(self.nDt, self.dt)
            traj.putData(data)
                

        else:
            raise self.ParamError(
                    "<float|function(x,t)|numpy.ndarray|Trajectory>")
        
        traj.ic=traj[0]
        traj.incrmTReal(finished=True, tReal=self.dt*self.nDt, t0=t0)
        if label<>None:traj.setLabel(label)
        return traj
            
    #----------------------------------------------------------------

    def _matricize(self, funcXT):
        time=np.linspace(0., self.nDt*self.dt, self.nDt+1)
        matrix=np.empty(shape=(self.nDt+1, self.grid.N))
        for i in xrange(self.grid.N):
            for j in xrange(self.nDt+1):
                matrix[j][i]=funcXT(self.grid.x[i], time[j])

        return matrix

    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="~~~~| Param |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        output+=self.grid.__str__()
        output+="\n nDt=%d"%self.nDt
        output+="\n dt=%-23.15E"%self.dt
        output+="\n tReal=%-23.15E"%self.tReal
        paramList=("forcing", "alpha", "beta", "gamma", "rho")
        for i in xrange(len(paramList)):
            output+="\n  :: "+paramList[i]+" ::"
            output+="\n  min= %f"%self[i].min()
            output+="\n  max= %f\n"%self[i].max()

        output+="\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        return output

    #----------------------------------------------------------------

    def __getitem__(self,i, j=None, k=None):
        if i==0:
            param=self.forcing
        elif i==1:
            param=self.alpha
        elif i==2:
            param=self.beta
        elif i==3:
            param=self.gamma
        elif i==4:
            param=self.rho
        elif i==5:
            param=self.fluxR
        else:
            raise IndexError('[0:forcing(x,t), 1:alpha(x,t),'+
                            '2:beta(x,t), 3:gamma(x,t), 4:rho(x,t), 5:fluxR(x,t)]')
        if j<>None:
            if k<>None:
                return (param.getData())[j][k]
            else:
                return (param.getData())[j]
        else:
            return param

    
        

#==========================================================
#----------------------------------------------------------
#==========================================================
#if __name__=='__main__':
