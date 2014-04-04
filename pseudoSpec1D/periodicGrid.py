import numpy as np
from grid import Grid

class PeriodicGrid(Grid):
    """
    Periodic Grid class
    for periodic 1+1D partial differential system

    PeriodicGrid(Ntrc, L, aliasing=3)

        Ntrc    :   truncature <int>
        L       :   periodic domain lenght <float>
        aliasing:   parameter to prevent aliasing <int>
                        (system dependent) 
                        N=aliasing*Ntrc+1 
                        reference: Durran, D. R. Numerical Methods
                            for Wave Equations in Geophysical Fluid
                            Dynamics, Springer-Verlag, 1999
    """
    
    class PeriodicGridError(Exception):
        pass
    
    #-------------------------------------------------------
    #----| Init |-------------------------------------------
    #-------------------------------------------------------
    def __init__(self, Ntrc, L=360., aliasing=3):
        if  not(type(Ntrc) is int) or not(isinstance(L, (float,int))):
            raise self.PeriodicGridError("PeriodicGrid(Ntrc <int>| L <float>)")
        
    
        self.Ntrc=Ntrc
        self.L=L

        if not(type(aliasing) is int):
            raise self.PeriodicGridError("PeriodicGrid(aliasing= <int>)")
        if aliasing<1:
            raise self.PeriodicGridError("aliasing > 1")
        else:
            self.aliasing=aliasing

        self.__setGridPoints()

        self.centered=True

    #-------------------------------------------------------
    #----| Public methods |---------------------------------
    #-------------------------------------------------------

    def mod(self, x):
        if x>self.max():
            diff=np.abs(self.max()-x)
            return self.x.min()+np.mod(diff, self.L)
        elif x<self.min():
            diff=np.abs(self.min()-x)
            return self.max()-np.mod(diff, self.L)
        else:
            pass
        return x

    #-------------------------------------------------------

    def plotPSpec(self, field,  axe=None, trunc=True, **kwargs):
        axe=self._checkAxe(axe)
        nDemi=int(self.N-1)/2
        data=np.zeros(nDemi)
        data=np.abs(np.fft.fft(field)[0:nDemi])
        k=Grid(nDemi,nDemi, centered=False)
        
        axe=k.plot(data, axe=axe, **kwargs)
        if trunc:
            axe.axvline(x=self.Ntrc, color='k', linestyle=':')
        return axe
    #-------------------------------------------------------
    #----| Private methods |--------------------------------
    #-------------------------------------------------------

    def __setGridPoints(self):
        self.N=self.aliasing*self.Ntrc+1
        if not self.N%2:
            raise self.PeriodicGridError(
                "N must be odd (add/substract 1 to Ntrc)")
        self.dx=self.L/(self.N)
        self.x=np.linspace(-self.L/2.,self.L/2-self.dx,self.N)


    #-------------------------------------------------------
    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="----| Grid |---------------------------\n"
        output+="| Ntrc=%d   N=%d\n"%(self.Ntrc, self.N)
        output+="| dx=%-23.15E\n"%(self.dx)
        output+="| L=%-23.15E\n"%(self.L)
        output+="---------------------------------------"
        return output


    def __eq__(self, grid2):
        value=(self.Ntrc==grid2.Ntrc and self.L==grid2.L)
        return value 
#--------------------------------------------------------------------
#####################################################################
#--------------------------------------------------------------------

if __name__=="__main__":
    g=PeriodicGrid(200,100.)
    print(g) 

