import numpy as np

cdef extern:
    void c_kdvpropagator(int N, int Ntrc, double L,
                            double dt, int nDt,
                            double* ic, double* traj, 
                            double* alph, double* beta, double* gamm,
                            double* rho, double* forc)

def fKdVPropagator(int N, int Ntrc, double L,
                    double dt, int nDt,
                    double[::1] ic not None,
                    double[::1] alph not None,
                    double[::1] beta not None,
                    double[::1] gamm not None,
                    double[::1] rho not None,
                    double[::1] forc not None):

    traj=np.empty(shape=(nDt+1, N))
    cdef double[:,::1] c_traj = traj

    c_kdvpropagator(N, Ntrc, L, dt, nDt, &ic[0], &c_traj[0,0], 
                    &alph[0], &beta[0], &gamm[0], &rho[0], &forc[0])

    return c_traj

#--------------------------------------------------------------------

cdef extern:
    void c_kdvtlmpropagator(int N, int Ntrc, double L,
                            double dt, int nDt,
                            double* u, double* p0, double* pTraj, 
                            double* alph, double* beta, double* gamm,
                            double* rho)

def fKdVTLMPropagator(int N, int Ntrc, double L,
                    double dt, int nDt,
                    double[::1] p0 not None,
                    double[:,::1] u not None,
                    double[::1] alph not None,
                    double[::1] beta not None,
                    double[::1] gamm not None,
                    double[::1] rho not None,
                    double[::1] forc not None):

    pTraj=np.empty(shape=(nDt+1, N))
    cdef double[:,::1] c_pTraj = pTraj

    c_kdvtlmpropagator(N, Ntrc, L, dt, nDt, &u[0,0], &p0[0], &c_pTraj[0,0], 
                        &alph[0], &beta[0], &gamm[0], &rho[0])

    return c_pTraj

#--------------------------------------------------------------------

cdef extern:
    void c_kdvtlmpropagatoradj(int N, int Ntrc, double L,
                                double dt, int nDt,
                                double* u, double* pf, double* aTraj, 
                                double* alph, double* beta, double* gamm,
                                double* rho)

def fKdVTLMPropagatorAdj(int N, int Ntrc, double L,
                            double dt, int nDt,
                            double[::1] pf not None,
                            double[:,::1] u not None,
                            double[::1] alph not None,
                            double[::1] beta not None,
                            double[::1] gamm not None,
                            double[::1] rho not None,
                            double[::1] forc not None):

    aTraj=np.empty(shape=(nDt+1, N))
    cdef double[:,::1] c_aTraj = aTraj

    c_kdvtlmpropagatoradj(N, Ntrc, L, dt, nDt, &u[0,0], &pf[0], 
                            &c_aTraj[0,0], 
                            &alph[0], &beta[0], &gamm[0], &rho[0])

    return c_aTraj

