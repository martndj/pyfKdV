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

        <TODO>  For now parameters (including forcing) are static
                although localized (Param.py)
                It would be interresting for parameters to be
                Trajectory object, although, the Fortran model
                would have to be profundly modified.

        <TODO>  parameters defined as array  

    """
    class ParamError(Exception):
        pass


    def __init__(self, grid, forcing=None, alpha=None, beta=None,
                    gamma=None,  rho=None):

        if not (isinstance(grid, PeriodicGrid)):
            raise self.ParamError("grid must be an instance of PeriodicGrid")

        self.grid=grid
        self.forcing=self.__setFunc(forcing)
        self.alpha=self.__setFunc(alpha)
        self.beta=self.__setFunc(beta)
        self.gamma=self.__setFunc(gamma)
        self.rho=self.__setFunc(rho)

        self.shape=(5,grid.N)

    #----------------------------------------------------------------

    def __setFunc(self, fct):
        def wrapper(c):
            def funcConst(x):
                return float(c)
            return funcConst
        
        if fct==None:
            return wrapper(0.)
        elif isinstance(fct, (int, float)):
            return wrapper(fct)
        elif callable(fct):
            return fct
        else:
            raise self.ParamError("<None|function>")
            


    #----------------------------------------------------------------

    def __getitem__(self,i):
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

    p=Param(grid, func2, func2, func1, func1, func3)
    plt.plot(grid.x, p[4])
    
    plt.show() 
