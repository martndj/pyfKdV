import numpy as np

class Grid(object):
    """
    Grid class
    1+1D partial differential system

    Grid(N, L, centered=True)

        N       :   number of grid points <int>
        L       :   domain lenght <float>
        centered:   centered grid: x in (-L/2, L/2)
                    not centered grid: x in (0., L)
    """
    
    class GridError(Exception):
        pass
    
    def __init__(self, N, L, centered=True):
        if  not(type(N) is int) or not(isinstance(L, (float,int))):
            raise self.GridError("Grid(N <int>| L <float|int>)")
        
    
        self.N=N
        self.L=L
        self.centered=centered

        self.__setGridPoints()
    
    #-------------------------------------------------------
    #----| Public methods |---------------------------------
    #-------------------------------------------------------
    def pos2Idx(self, pos):
        """
        Convert space position to grid index

            PeriodicGrid.pos2Idx(pos)

            pos :   array of positions <numpy.ndarray>
        """
        if isinstance(pos, np.ndarray):
            if pos.ndim<>1:
                raise self.PeriodicGridError("pos.ndim=1")
            N=len(pos)
            idx=np.empty(N, dtype=int)
            for i in xrange(N):
                idx[i]=np.min(np.where(self.x>=pos[i]))
        elif isinstance(pos, (float, int)):
            idx=np.empty(1, dtype=int)
            idx[0]=np.min(np.where(self.x>=pos))
        else:
            raise self.PeriodicGridError("pos <numpy.ndarray>")
        return idx 
    #-------------------------------------------------------
    #----| Private methods |--------------------------------
    #-------------------------------------------------------

    def __setGridPoints(self):
        self.dx=self.L/(self.N)
        if self.centered:
            self.x=np.linspace(-self.L/2.,self.L/2-self.dx,self.N)
        else:
            self.x=np.linspace(0.,self.L-self.dx,self.N)
            


    #-------------------------------------------------------
    #----| Classical overloads |----------------------------
    #-------------------------------------------------------

    def __str__(self):
        output="----| Grid |---------------------------\n"
        output+="| N=%d\n"%(self.N)
        output+="| dx=%-23.15E\n"%(self.dx)
        output+="| L=%-23.15E\n"%(self.L)
        output+="---------------------------------------"
        return output


    def __eq__(self, grid2):
        value=(self.N==grid2.N and self.L==grid2.L)
        return value 
#--------------------------------------------------------------------
#####################################################################
#--------------------------------------------------------------------

if __name__=="__main__":
    g=Grid(200,100.)
    print(g) 

