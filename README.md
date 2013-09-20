Korteweg de-Vries Pseudospectral Integrator (pyfKdV)
====================================================

Martin Deshaies-Jacques

[deshaies.martin@sca.uqam.ca](mailto:deshaies.martin@sca.uqam.ca)

[www.science.martn.info](http://www.science.martn.info)

pyfKdV is a Fortran 95 pseudospectral propagator wrapped with python 2.7.

Program released in this bundle are licenced under the terms of the opensource GPL licence version 3 (http://www.gnu.org/licenses/gpl.html).
Detailed can be found in ./doc/gpl.txt.

Use, share, remix!
...and criticize


See also:

 * Data Assimilation Lab: [dVar](https://github.com/martndj/dVar)

About Korteweg-de-Vrie _augmented_ system and it's use in atmospheric dynamics
---------------------------------------------------------------------------

 * [Warn, T. and Brasnett, B. The amplification and capture of atmospheric solitons by topography:38, 1982.](http://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1983\)040%3C0028%3ATAACOA%3E2.0.CO%3B2)
 * [Hodyss, D. and Nathan, T. R. Solitary rossby waves in zonally varying jet ows. Geophys. Astrophys., 262, 2002.](http://www.tandfonline.com/doi/abs/10.1080/03091920290011012#.Ug1egSHPTMU)


Compiling informations
----------------------
### Fortran Dependencies
 * [FFTW3](http://www.fftw.org/): fortran 77 fast-fourier transform librairy (http://www.fftw.org/)
 * [ARPACK](http://www.caam.rice.edu/software/ARPACK/) is a collection of Fortran77 subroutines designed to solve large scale eigenvalue problems.

### Python Dependencies
 * Numpy
 * Matplotlib
 * Cython

### Instructions
 * ./ refer to the root of pyfKdV installation;
 * [something] means optional arguments;
 * \<something\> means you must replace it with what is appropriate to your environment.


 1. On Linux OS Debian/Ubuntu, you can install it running

        [sudo] apt-get install libfftw3* libarpack2*
        [sudo] apt-get install python python-matplotlib python-numpy

    The version of cython shipping with Debian at the moment (0.15.1) is not cutting edge enough (missing memoryviews component.)
    You should [download the latest](http://cython.org/#download).
    (I used [0.19.1](http://cython.org/release/Cython-0.19.1.tar.gz) personnaly).
    Installing is straightforward, instruction are included.
    You'll need python header files:
    
        [sudo] apt-get install python-dev
        
    (there was a ascII error and I had to include a --no-compile-something to the setup.py call but it work perfectly afterward.)


 2. To compile, you'll have to modify ./src/Makefile, to specify your libraries directories, and fortran compiler:
    
        libDir=/usr/lib
        userLibDir=/home/<user>/lib
        compiler=gfortran

   Be sure that ${userLibDir} is in your ${LD_LIBRARY_PATH}:

        echo $LD_LIBRARY_PATH

   if not:

        export LD_LIBRARY_PATH=".:/home/<user>/lib/"

 3. Build the ./src/ directory

        make [tests] pyKdV [clean]


 4. Add the path to the python modules pyKdV and pseudoSpec1D to your PYTHONPATH environment variable, in Linux bash environment you can do it running.
 
        export PYTHONPATH=${PYTHONPATH}:<path to ./>

    Adding this export lines to your startup script (.bashrc or .profile) is a way to do it.

 4. Have fun!



Python Wrapper
--------------
Configuring and launching an integration
----------------------------------------

 1. Write a launcher python script:
  
        import numpy as np
        import pyKdV as kdv
        import matplotlib.pyplot as plt 
            
        #----| Grid configuration |-------------------
        grid=kdv.PeriodicGrid(150,300.)
        tInt=30.
        maxA=4.
            
        #----| KdV parameters arguments |-------------
        def gauss(x):
            x0=0.
            sig=5.
            return -0.1*np.exp(-((x-x0)**2)/(2*sig**2))
                
        param=kdv.Param(grid, beta=1., gamma=-1, rho=gauss)
            
        #----| Initial condition |--------------------
        ic=kdv.rndFiltVec(grid, Ntrc=grid.Ntrc/4)\
                    +kdv.soliton(grid.x, 0., amp=3.)
            
        #----| Launching the integration |------------
        launcher=kdv.kdvLauncher(param, maxA)
        traj=launcher.integrate(ic, tInt)
            
        
        #----| Perturbating the trajectory |----------
        pert=0.5*kdv.gauss(grid.x, 0., 20.)
        #-- linear perturbation
        tlmLauncher=kdv.kdvTLMLauncher(param)
        tlmLauncher.initialize(traj)
        fLinearPert=tlmLauncher.integrate(pert, tInt)
        #-- nonlinear perturbation
        fNLPert=launcher.integrate(ic+pert, tInt).final()-traj.final()
        
        #----| Plotting the result |------------------
        subplt1=plt.subplot(311)
        subplt2=plt.subplot(312)
        subplt3=plt.subplot(313)
        traj.waterfall(axe=subplt1)
        subplt1.legend([r"$x(t)$"], loc="lower right")
        subplt2.plot(grid.x, fLinearPert)
        subplt2.plot(grid.x, fNLPert)
        subplt2.legend(["$\mathbf{L}\delta x$", 
                        "$\mathcal{M}(x+\delta x)-\mathcal{M}(x)$"],
                        loc='lower right')
        traj.fftTraj().waterfall(axe=subplt3, color='r')
        subplt3.legend([ "$N_{trc}$", r"$\mathcal{F}[x(t)]$"], loc="lower right")
        
        plt.show()
  

 2. Singular vector calculation work similarly, but are way longer to obtain:

        svLauncher=kdv.kdvSVLauncher(traj, param)
        sVal=svLauncher.lanczos(3, tInt=1.)
        print(sVal)
        plt.figure()
        for sv in svLauncher.sVec:
            plt.plot(grid.x, sv) 


Acknoledgement
--------------
A special thank to [Mauro Werder](http://www.sfu.ca/~mawerder/), through the blog of whom I could understand [how to link python and fortran](http://www.sfu.ca/~mawerder/notes/calling_fortran_from_python.html)!
