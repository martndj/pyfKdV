module kdvWrapC

use iso_c_binding, only: c_double, c_int
use kdvProp

contains
!====================================================================

subroutine c_kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, traj, &
                             alph, beta, gamm, rho, forc) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double), intent(inout)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    forc
    real(c_double), intent(inout), dimension(N) ::  ic

    ! note that in C the indices will be reversed!:
    real(c_double), intent(out), dimension(N, nDt+1)    ::  traj
    
    ! ...hence the transpose:
    traj=transpose(kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, &
                                    alph, beta, gamm, rho, forc))

end subroutine c_kdvPropagator

!====================================================================
end module kdvWrapC
