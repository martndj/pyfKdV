import numpy as np
from grid import Grid
from Launcher import Launcher
from TLMLauncher import TLMLauncher

class IdLauncher(Launcher):
    """
    Dummy Launcher (identity)

    Does nothing!  Its purpose is to tests data assimilation schemes.
    """

    def __init__(self, grid, dt=1.):
        if not isinstance(grid, Grid):
            raise TypeError("grid <Grid>")

        self.grid=grid
        self.dt=dt
        self.propagator=self._identity

    def _identity(self, ic, traj, t0=0.):
        trajData=np.zeros(shape=(traj.nDt+1, self.grid.N))
        for i in xrange(traj.nDt+1):
            trajData[i]=ic
        
        tReal=traj.nDt*traj.dt
        traj.putData(trajData)
        traj.incrmTReal(finished=True, tReal=tReal, t0=t0)
        return traj

class IdTLM(IdLauncher, TLMLauncher):
    def __init__(self, grid, traj=None):
        if not isinstance(grid, Grid):
            raise TypeError("grid <Grid>")

        self.grid=grid

        self.isReferenced=False 
        if not traj==None: self.reference(traj)

        self.propagator=self._identity
        self.propagatorAdj=self._identity
