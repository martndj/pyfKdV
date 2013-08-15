import numpy as np

cdef extern:
    void c_kdvpropagator(int N, int Ntrc, double L,
                            double dt, int nDt, double tReal,
                            double* ic, double* traj, 
                            double* alph, double* beta, double* gamm,
                            double* rho, double* forc)

def fKdVPropagator(int N, int Ntrc, double L,
                    double dt, int nDt, double tReal,
                    double[::1] ic not None,
                    double[::1] alph not None,
                    double[::1] beta not None,
                    double[::1] gamm not None,
                    double[::1] rho not None,
                    double[::1] forc not None):

    traj=np.empty(shape=(nDt+1, N))
    cdef double[:,::1] c_traj = traj

    c_kdvpropagator(N, Ntrc, L, dt, nDt, tReal, &ic[0], &c_traj[0,0], 
                    &alph[0], &beta[0], &gamm[0], &rho[0], &forc[0])

    return c_traj
