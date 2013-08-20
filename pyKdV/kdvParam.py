import numpy as np
from pseudoSpec1D import *

class Param(object):
    """
     augmented dufferential system parameters

        \partial_t A(x,t)= forcing(x,t) - \alpha(x,t) A 
                           - \beta(x,t)A \partial_x A 
                           - \gamma(x,t) \partial_x^3 A - \rho(x,t) A 

        <TODO>  For now parameters (including forcing) are static
                although localized (Param.py)

        <TODO>  forcing, alpha, ...
                could either be:
                    * constant
                    * vectors (localized static parameters)
                    * Trajectory object
                    * user defined function param(x,t)
                        <!> x and t must be explicit
                            so that the type be conserved
                        def func(x,t):
                            return 0.*x + 0.1*t**2
                the later being the more general, we'll use this one
                for start

        <TODO>  if param=None return null function
    """
    class ParamError(Exception):
        pass


    def __init__(self, grid, alpha=None, beta=None, gamma=None, rho=None,
                    forcing=None):

        if not (isinstance(grid, SpectralGrid)):
            raise self.ParamError("grid must be an instance of SpectralGrid")

        self.grid=grid
        self.forcing=self.__setFunc(forcing)
        self.alpha=self.__setFunc(alpha)
        self.beta=self.__setFunc(beta)
        self.gamma=self.__setFunc(gamma)
        self.rho=self.__setFunc(rho)

        self.shape=(5,grid.N)

    #----------------------------------------------------------------

    def __setFunc(self, fct):
        def funcNulle(x):
            return 0.
        if fct==None:
            return funcNulle
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


#==========================================================
#----------------------------------------------------------
#==========================================================
if __name__=='__main__':
    import matplotlib.pyplot as plt
    
    grid=SpectralGrid(142, 5.)
    
    def func1(x):
        return  -1.
    
    def func2(x):
        return 0.
    
    def func3(x):
        return np.sin(x)

    p=Param(grid, func2, func2, func1, func1, func3)
    plt.plot(grid.x, p[4])
    
    plt.show() 
