import numpy as np

class SpectralGrid:
    """
    Spectral Grid class
    for 1+1D partial differential system

        Light version for pyKdV5
            * no I\O methods


        <TODO>  include possibility to change N(Ntrc) relation
                which is problem dependant
                __setGridPoints()
    """

    __version__='Light version for pyKdV5'
    __author__='Martin Deshaies-Jacques <deshaies.martin@sca.uqam.ca>'
    
    class SpectralGridError(Exception):
        pass
    
    def __init__(self, Ntrc, L):
        """
        Grid default constructor
        """
        if  not(type(Ntrc) is int) or not(type(L)is float):
            raise self.SpectralGridError("SpectralGrid(Ntrc <int>| L <float>)")
        
    
        self.Ntrc=Ntrc
        self.L=L

        self.__setGridPoints()
        self.dx=self.L/(self.N)

        self.x=np.linspace(-self.L/2.,self.L/2,self.N)
    
    #-------------------------------------------------------
    #----| Private methods |--------------------------------
    #-------------------------------------------------------

    def __setGridPoints(self):
        self.N=3*self.Ntrc+1
        if not self.N%2:
            raise self.SpectralGridError(
                "N must be odd (add/substract 1 to Ntrc)")


    #-------------------------------------------------------
    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        """
        """
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
    g=SpectralGrid(200,100.)
    print(g) 

