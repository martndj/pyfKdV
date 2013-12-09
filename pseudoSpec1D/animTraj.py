from trajectory import Trajectory
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

def animTraj(traj, dt, tf=None, interval=25, fig=None, axes=None,
                filename=None, ylim=None):

    if not isinstance(traj, Trajectory):
        raise TypeError()

    if fig==None or axes==None:
        fig, ax = plt.subplots()

    dN=traj.whereTimeIdx(dt)
    if tf==None:
        IdxF=traj.shape[0]
    else:
        IdxF=traj.whereTimeIdx(tf)

    if ylim==None:
        ax.set_ylim(traj.min(), traj.max())
    else:   
        ax.set_ylim(ylim)

    line, = ax.plot(traj.grid.x, traj[0]) 
    def animate(i):
        line.set_ydata(traj[i])
        return line,

    def init():
        line.set_ydata(np.ma.array(traj.grid.x, mask=True))
        return line,

    ani = animation.FuncAnimation(fig, animate, 
            np.arange(0, IdxF, dN), interval= 25, 
            blit=True, init_func=init)
    if filename<>None:
        ani.save(filename, fps=15)
    
    return ani

