module kdvWrapC

use iso_c_binding, only: c_double, c_int, c_bool
use kdvProp
use kdvTLMProp
use kdvTLMTest
use kdvLanczos

contains
!====================================================================


subroutine c_kdvPropagator(N, Ntrc, L, dt, nDt, nDtParam, ic, traj, &
                                alph, beta, gamm, rho, forc) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt,&
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), dimension(N) ::  ic

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam+1)&
                                        ::  alph, beta, gamm, rho, forc
    real(c_double), intent(out), dimension(N, nDt+1)    ::  traj
    
    ! ...hence the transpose:
    traj=transpose(kdvPropagator(N, Ntrc, L, dt, nDt, nDtParam, &
                                    tReal, ic, &
                                    transpose(alph), transpose(beta),&
                                    transpose(gamm), transpose(rho),&
                                    transpose(forc))&
                                    )

end subroutine c_kdvPropagator


!--------------------------------------------------------------------
!--------------------------------------------------------------------


subroutine c_kdvTLMPropagator(N, Ntrc, L, dt, nDt, nDtParam,&
                                    u, p0, pf, &
                                    alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt,&
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

                                                   
    real(c_double), intent(in), dimension(N)   ::  p0
    real(c_double), intent(out), dimension(N)   ::  pf

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam+1)&
                                            ::  alph, beta, gamm, rho 
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    pf=kdvTLMPropagator(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                        transpose(u), p0, transpose(alph), &
                        transpose(beta), transpose(gamm), &
                        transpose(rho))

end subroutine c_kdvTLMPropagator

!--------------------------------------------------------------------

subroutine c_kdvTLMPropagatorFullTraj(N, Ntrc, L, dt, nDt, nDtParam,&
                                         u, p0, pf, pTraj, alph, &
                                         beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, &
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)   ::  p0
    real(c_double), intent(out), dimension(N)   ::  pf

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam+1)&
                                            ::  alph, beta, gamm, rho 
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    real(c_double), intent(out), dimension(N, nDt+1)    ::  pTraj
    real(c_double), dimension(nDt+1, N)    ::  f_pTraj
    
    ! ...hence the transpose:
    pf=kdvTLMPropagator(N, Ntrc, L, dt, nDt, nDtParam, tReal, &
                        transpose(u), p0, transpose(alph), &
                        transpose(beta), transpose(gamm), &
                        transpose(rho), f_pTraj)

    pTraj=transpose(f_pTraj)

end subroutine c_kdvTLMPropagatorFullTraj

!--------------------------------------------------------------------
!--------------------------------------------------------------------


subroutine c_kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, nDtParam, &
                                    u, pf, adj, &
                                    alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, &
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  pf
    real(c_double), intent(out), dimension(N)   ::  adj

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam+1)&
                                            ::  alph, beta, gamm, rho
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    adj=kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                        transpose(u), pf, transpose(alph), &
                        transpose(beta), transpose(gamm), &
                        transpose(rho))

end subroutine c_kdvTLMPropagatorAdj

!--------------------------------------------------------------------

subroutine c_kdvTLMPropagatorAdjFullTraj(N, Ntrc, L, dt, nDt,nDtParam, &
                                            u, pf, adj, aTraj, alph,&
                                            beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, &
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  pf
    real(c_double), intent(out), dimension(N)   ::  adj

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam+1)&
                                            ::  alph, beta, gamm, rho
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    real(c_double), intent(out), dimension(N, nDt+1)    ::  aTraj
    real(c_double), dimension(nDt+1, N)    ::  f_aTraj
    
    ! ...hence the transpose:
    adj=kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, nDtParam, tReal, &
                        transpose(u), pf, transpose(alph), &
                        transpose(beta), transpose(gamm), &
                        transpose(rho), f_aTraj)
    aTraj=transpose(f_aTraj)

end subroutine c_kdvTLMPropagatorAdjFullTraj

!--------------------------------------------------------------------
!--------------------------------------------------------------------


subroutine c_kdvTLMSingularOp(N, Ntrc, L, dt, nDt, nDtParam, &
                                    u, x, y, &
                                    alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, &
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  x
    real(c_double), intent(out), dimension(N)   ::  y

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam+1) &
                                            ::  alph, beta, gamm, rho
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    y=kdvTLMSingularOp(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                        transpose(u), x, transpose(alph), &
                        transpose(beta), transpose(gamm), &
                        transpose(rho))

end subroutine c_kdvTLMSingularOp

!--------------------------------------------------------------------
!--------------------------------------------------------------------


subroutine c_kdvLanczos(N, Ntrc, L, dt, nDt, nDtParam, u, Nev, V, sv, &
                             alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, Nev,&
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(out), dimension(Nev)     ::  sv
    
    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam+1) &
                                            ::  alph, beta, gamm, rho
    real(c_double), intent(out), dimension(Nev, N)   ::  V
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    V=transpose(lanczos(N, Ntrc, L, dt, nDt, nDtParam, tReal, &
                transpose(u), transpose(alph), transpose(beta), &
                transpose(gamm), transpose(rho), Nev, sv, Nconv))

end subroutine c_kdvLanczos

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine c_testGradient(N, Ntrc, L, dt, nDt, nDtParam, maxPower, &
                            p, alph, beta, gamm, rho, forc) bind(c)
    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, &
                                                    nDtParam, maxPower
    real(c_double), intent(in), value           ::  L, dt

    real(c_double), intent(in), dimension(N)    ::  p
    
    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam+1) &
                                        ::  alph, beta, gamm, rho, forc
    call NLTestGradient(N, Ntrc, L, dt, nDt, nDtParam, maxPower, &
                            p, transpose(alph), transpose(beta), &
                            transpose(gamm), transpose(rho), &
                            transpose(forc))

end subroutine c_testGradient
!====================================================================
end module kdvWrapC
