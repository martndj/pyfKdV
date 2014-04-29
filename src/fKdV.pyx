import numpy as np

cdef extern:
    void c_kdvpropagator(int N, int Ntrc, double L,
                            double dt, int nDt, int nDtParam, 
                            double* ic, double* traj, 
                            double* alph, double* beta, double* gamm,
                            double* rho, double* forc)

def fKdVPropagator(int N, int Ntrc, double L,
                    double dt, int nDt, int nDtParam, 
                    double[::1] ic not None,
                    double[:,::1] alph not None,
                    double[:,::1] beta not None,
                    double[:,::1] gamm not None,
                    double[:,::1] rho not None,
                    double[:,::1] forc not None):

    traj=np.empty(shape=(nDt+1, N))
    cdef double[:,::1] c_traj = traj

    c_kdvpropagator(N, Ntrc, L, dt, nDt, nDtParam,
                    &ic[0], &c_traj[0,0], 
                    &alph[0,0], &beta[0,0], &gamm[0,0], &rho[0,0],
                    &forc[0,0])

    return np.array(c_traj)

#--------------------------------------------------------------------
#--------------------------------------------------------------------

cdef extern:
    void c_kdvtlmpropagator(int N, int Ntrc, double L,
                            double dt, int nDt, int nDtParam, 
                            double* u,
                            double* p0, double* pf, 
                            double* alph, double* beta, double* gamm,
                            double* rho)

cdef extern:
    void c_kdvtlmpropagatorfulltraj(int N, int Ntrc, double L,
                                    double dt, int nDt, int nDtParam, 
                                    double* u, double* p0, double* pf, 
                                    double* pTraj, 
                                    double* alph, double* beta,
                                    double* gamm, double* rho)

def fKdVTLMPropagator(int N, int Ntrc, double L,
                    double dt, int nDt, int nDtParam, 
                    double[::1] p0 not None,
                    double[:,::1] u not None,
                    double[:,::1] alph not None,
                    double[:,::1] beta not None,
                    double[:,::1] gamm not None,
                    double[:,::1] rho not None,
                    fullTraj=False):

    pTraj=np.empty(shape=(nDt+1, N))
    cdef double[:,::1] c_pTraj = pTraj
    cdef double[::1] pf = np.empty(N)

    if fullTraj:
        c_kdvtlmpropagatorfulltraj(N, Ntrc, L, dt, nDt, nDtParam,
                                    &u[0,0], &p0[0], &pf[0], 
                                    &c_pTraj[0,0], &alph[0,0], 
                                    &beta[0,0], &gamm[0,0], &rho[0,0])

        return np.array(c_pTraj)
    else:
        c_kdvtlmpropagator(N, Ntrc, L, dt, nDt, nDtParam, &u[0,0], 
                            &p0[0], &pf[0], &alph[0,0], &beta[0,0],
                            &gamm[0,0], &rho[0,0])
        return np.array(pf)



#--------------------------------------------------------------------
#--------------------------------------------------------------------

cdef extern:
    void c_kdvtlmpropagatoradj(int N, int Ntrc, double L,
                            double dt, int nDt, int nDtParam, double* u,
                            double* pf, double* adj, 
                            double* alph, double* beta, double* gamm,
                            double* rho)

cdef extern:
    void c_kdvtlmpropagatoradjfulltraj(int N, int Ntrc, double L,
                                    double dt, int nDt, int nDtParam, 
                                    double* u, double* pf, double* adj, 
                                    double* aTraj, 
                                    double* alph, double* beta,
                                    double* gamm, double* rho)

def fKdVTLMPropagatorAdj(int N, int Ntrc, double L,
                    double dt, int nDt, int nDtParam, 
                    double[::1] pf not None,
                    double[:,::1] u not None,
                    double[:,::1] alph not None,
                    double[:,::1] beta not None,
                    double[:,::1] gamm not None,
                    double[:,::1] rho not None,
                    fullTraj=False):

    aTraj=np.empty(shape=(nDt+1, N))
    cdef double[:,::1] c_aTraj = aTraj
    cdef double[::1] adj= np.empty(N)

    if fullTraj:
        c_kdvtlmpropagatoradjfulltraj(N, Ntrc, L, dt, nDt, nDtParam,
                                    &u[0,0], &pf[0], &adj[0], 
                                    &c_aTraj[0,0], &alph[0,0], 
                                    &beta[0,0], &gamm[0,0], &rho[0,0])

        return np.array(c_aTraj)
    else:
        c_kdvtlmpropagatoradj(N, Ntrc, L, dt, nDt, nDtParam, &u[0,0], 
                        &pf[0], &adj[0], &alph[0,0], &beta[0,0],
                        &gamm[0,0], &rho[0,0])
        return np.array(adj)

#--------------------------------------------------------------------
#--------------------------------------------------------------------


cdef extern:
    void c_kdvtlmsingularop(int N, int Ntrc, double L,
                            double dt, int nDt, int nDtParam, double* u,
                            double* x, double* y, 
                            double* alph, double* beta, double* gamm,
                            double* rho)

def fKdVTLMSingularOp(int N, int Ntrc, double L,
                    double dt, int nDt, int nDtParam,
                    double[::1] x not None,
                    double[:,::1] u not None,
                    double[:,::1] alph not None,
                    double[:,::1] beta not None,
                    double[:,::1] gamm not None,
                    double[:,::1] rho not None):

    cdef double[::1] y = np.empty(N)

    c_kdvtlmsingularop(N, Ntrc, L, dt, nDt, nDtParam, &u[0,0], 
                        &x[0], &y[0], &alph[0,0], &beta[0,0], 
                        &gamm[0,0], &rho[0,0])
    return np.array(y)

#--------------------------------------------------------------------
#--------------------------------------------------------------------

cdef extern:
    void c_kdvlanczos(int N, int Ntrc, double L,
                            double dt, int nDt, int nDtParam, double* u,
                            int Nev, double* V, double* sv, 
                            double* alph, double* beta, double* gamm,
                            double* rho)

def fKdVLanczos(int N, int Ntrc, double L,
                    double dt, int nDt, int nDtParam,
                    double[:,::1] u not None,
                    int Nev,
                    double[:,::1] alph not None,
                    double[:,::1] beta not None,
                    double[:,::1] gamm not None,
                    double[:,::1] rho not None):

    cdef double[:,::1] V = np.empty(shape=(N, Nev))
    cdef double[::1] sv = np.empty(Nev)

    c_kdvlanczos(N, Ntrc, L, dt, nDt, nDtParam, &u[0,0], 
                    Nev, &V[0,0], &sv[0], 
                    &alph[0,0], &beta[0,0], &gamm[0,0], &rho[0,0])

    # descending sv order
    return np.array(sv[::-1]), np.array(V[:,::-1]).transpose()

#--------------------------------------------------------------------
#--------------------------------------------------------------------

cdef extern:
    void c_testgradient(int N, int Ntrc, double L,
                            double dt, int nDt, int nDtParam, 
                            int maxPower,  double* p, 
                            double* alph, double* beta, double* gamm,
                            double* rho, double* forc)

def fKdVTestGradient(int N, int Ntrc, double L,
                    double dt, int nDt, int nDtParam, int maxPower,
                    double[::1] p not None,
                    double[:,::1] alph not None,
                    double[:,::1] beta not None,
                    double[:,::1] gamm not None,
                    double[:,::1] rho not None,
                    double[:,::1] forc not None):


    c_testgradient(N, Ntrc, L, dt, nDt, nDtParam, maxPower, 
                    &p[0], &alph[0,0], &beta[0,0], &gamm[0,0],
                    &rho[0,0], &forc[0,0])
