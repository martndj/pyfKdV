module kdvWrapC

use iso_c_binding, only: c_double, c_int
use kdvProp
use kdvTLMProp

contains
!====================================================================

subroutine c_kdvPropagator(N, Ntrc, L, dt, nDt, ic, traj, &
                             alph, beta, gamm, rho, forc) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    forc
    real(c_double), dimension(N) ::  ic

    ! note that in C the indices will be reversed!:
    real(c_double), intent(out), dimension(N, nDt+1)    ::  traj
    
    ! ...hence the transpose:
    traj=transpose(kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, &
                                    alph, beta, gamm, rho, forc))

end subroutine c_kdvPropagator


!--------------------------------------------------------------------


subroutine c_kdvTLMPropagator(N, Ntrc, L, dt, nDt, u, p0, pTraj, &
                             alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    p0

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    real(c_double), intent(out), dimension(N, nDt+1)    ::  pTraj
    
    ! ...hence the transpose:
    pTraj=transpose(kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal,&
                                    transpose(u), p0, &
                                    alph, beta, gamm, rho))

end subroutine c_kdvTLMPropagator


!--------------------------------------------------------------------


subroutine c_kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, u, pf, aTraj, &
                                    alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    pf

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    real(c_double), intent(out), dimension(N, nDt+1)    ::  aTraj
    
    ! ...hence the transpose:
    aTraj=transpose(kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, tReal,&
                                    transpose(u), pf, &
                                    alph, beta, gamm, rho))

end subroutine c_kdvTLMPropagatorAdj


!====================================================================
end module kdvWrapC
