module kdvWrapC

use iso_c_binding, only: c_double, c_int, c_bool
use kdvProp
use kdvProp_pt
use kdvTLMProp
use kdvTLMProp_pt
use kdvTLMTest
use kdvLanczos

contains
!====================================================================

subroutine c_kdvPropagator(N, Ntrc, L, dt, nDt, ic, traj, &
                             alph, beta, gamm, rho, forc) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, &
                                                    rho, forc
    real(c_double), dimension(N) ::  ic

    ! note that in C the indices will be reversed!:
    real(c_double), intent(out), dimension(N, nDt+1)    ::  traj
    
    ! ...hence the transpose:
    traj=transpose(kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, &
                                    alph, beta, gamm, rho, forc))

end subroutine c_kdvPropagator


!--------------------------------------------------------------------

subroutine c_kdvPropagator_pt(N, Ntrc, L, dt, nDt, nDtParam, ic, traj, &
                                alph, beta, gamm, rho, forc) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt,&
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), dimension(N) ::  ic

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam)&
                                        ::  alph, beta, gamm, rho, forc
    real(c_double), intent(out), dimension(N, nDt+1)    ::  traj
    
    ! ...hence the transpose:
    traj=transpose(kdvPropagator_pt(N, Ntrc, L, dt, nDt, nDtParam, &
                                    tReal, ic, &
                                    transpose(alph), transpose(beta),&
                                    transpose(gamm), transpose(rho),&
                                    transpose(forc))&
                                    )

end subroutine c_kdvPropagator_pt


!--------------------------------------------------------------------
!--------------------------------------------------------------------


subroutine c_kdvTLMPropagator(N, Ntrc, L, dt, nDt, u, p0, pf, &
                             alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    p0
    real(c_double), intent(out), dimension(N)   ::  pf

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    pf=kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal,&
                        transpose(u), p0, alph, beta, gamm, rho)

end subroutine c_kdvTLMPropagator

!--------------------------------------------------------------------

subroutine c_kdvTLMPropagator_pt(N, Ntrc, L, dt, nDt, nDtParam,&
                                    u, p0, pf, &
                                    alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt,&
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

                                                   
    real(c_double), intent(in), dimension(N)   ::  p0
    real(c_double), intent(out), dimension(N)   ::  pf

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam)&
                                            ::  alph, beta, gamm, rho 
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    pf=kdvTLMPropagator_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                        transpose(u), p0, transpose(alph), &
                        transpose(beta), transpose(gamm), &
                        transpose(rho))

end subroutine c_kdvTLMPropagator_pt

!--------------------------------------------------------------------

subroutine c_kdvTLMPropagatorFullTraj(N, Ntrc, L, dt, nDt, u, p0, pf, &
                        pTraj, alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    p0
    real(c_double), intent(out), dimension(N)   ::  pf

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    real(c_double), intent(out), dimension(N, nDt+1)    ::  pTraj
    real(c_double), dimension(nDt+1, N)    ::  f_pTraj
    
    ! ...hence the transpose:
    pf=kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal, &
                        transpose(u), p0, alph, beta, gamm, rho, f_pTraj)
    pTraj=transpose(f_pTraj)

end subroutine c_kdvTLMPropagatorFullTraj

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine c_kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, u, pf, adj, &
                             alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    pf
    real(c_double), intent(out), dimension(N)   ::  adj

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    adj=kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, tReal,&
                        transpose(u), pf, alph, beta, gamm, rho)

end subroutine c_kdvTLMPropagatorAdj

!--------------------------------------------------------------------

subroutine c_kdvTLMPropagatorAdj_pt(N, Ntrc, L, dt, nDt, nDtParam, &
                                    u, pf, adj, &
                                    alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, &
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  pf
    real(c_double), intent(out), dimension(N)   ::  adj

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam)&
                                            ::  alph, beta, gamm, rho
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    adj=kdvTLMPropagatorAdj_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                        transpose(u), pf, transpose(alph), &
                        transpose(beta), transpose(gamm), &
                        transpose(rho))

end subroutine c_kdvTLMPropagatorAdj_pt

!--------------------------------------------------------------------

subroutine c_kdvTLMPropagatorAdjFullTraj(N, Ntrc, L, dt, nDt, u, pf, & 
                             adj, aTraj, alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    pf
    real(c_double), intent(out), dimension(N)   ::  adj

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    real(c_double), intent(out), dimension(N, nDt+1)    ::  aTraj
    real(c_double), dimension(nDt+1, N)    ::  f_aTraj
    
    ! ...hence the transpose:
    adj=kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, tReal, &
                        transpose(u), pf, alph, beta, gamm, rho, f_aTraj)
    aTraj=transpose(f_aTraj)

end subroutine c_kdvTLMPropagatorAdjFullTraj

!--------------------------------------------------------------------
!--------------------------------------------------------------------


subroutine c_kdvTLMSingularOp(N, Ntrc, L, dt, nDt, u, x, y, &
                             alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho, &
                                                    x
    real(c_double), intent(out), dimension(N)   ::  y

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    y=kdvTLMSingularOp(N, Ntrc, L, dt, nDt, tReal,&
                        transpose(u), x, alph, beta, gamm, rho)

end subroutine c_kdvTLMSingularOp

!--------------------------------------------------------------------


subroutine c_kdvTLMSingularOp_pt(N, Ntrc, L, dt, nDt, nDtParam, &
                                    u, x, y, &
                                    alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, &
                                                    nDtParam
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  x
    real(c_double), intent(out), dimension(N)   ::  y

    ! note that in C the indices will be reversed!:
    real(c_double), intent(in), dimension(N, nDtParam) &
                                            ::  alph, beta, gamm, rho
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    y=kdvTLMSingularOp_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                        transpose(u), x, transpose(alph), &
                        transpose(beta), transpose(gamm), &
                        transpose(rho))

end subroutine c_kdvTLMSingularOp_pt

!--------------------------------------------------------------------
!--------------------------------------------------------------------


subroutine c_kdvLanczos(N, Ntrc, L, dt, nDt, u, Nev, V, sv, &
                             alph, beta, gamm, rho) bind(c)

    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, Nev
    real(c_double), intent(in), value           ::  L, dt
    real(c_double)               ::  tReal

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho
                                                  
    real(c_double), intent(out), dimension(Nev)     ::  sv
    
    ! note that in C the indices will be reversed!:
    real(c_double), intent(out), dimension(Nev, N)   ::  V
    real(c_double), intent(in), dimension(N, nDt+1)     ::  u
    
    ! ...hence the transpose:
    V=transpose(lanczos(N, Ntrc, L, dt, nDt, tReal, transpose(u), &
                    alph, beta, gamm, rho, Nev, sv, Nconv))

end subroutine c_kdvLanczos

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine c_testGradient(N, Ntrc, L, dt, nDt, maxPower, &
                            p, alph, beta, gamm, rho, forc) bind(c)
    integer(c_int), intent(in), value           ::  N, Ntrc, nDt, maxPower
    real(c_double), intent(in), value           ::  L, dt

    real(c_double), intent(in), dimension(N)    ::  alph, beta, gamm, rho,&
                                                    forc, p
    
    call NLTestGradient(N, Ntrc, L, dt, nDt, maxPower, &
                            p, alph, beta, gamm, rho, forc)

end subroutine c_testGradient
!====================================================================
end module kdvWrapC
