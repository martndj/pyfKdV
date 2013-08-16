module kdvTLMProp

use spectral
use kdvProp, only: rhoCenteredImplicit
implicit none

contains
!====================================================================


function kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal, u, p0, &
                            alph, beta, gamm, rho) &
                            result(pTraj)

    intent(in)              ::  N, Ntrc, L, dt, nDt, u, &
                                alph, beta, gamm, rho
    !intent(inout)           ::  p0, tReal

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, j
    
    double precision, dimension(N)          ::  alph, beta, gamm, rho, p0

    double precision, dimension(3, N)       ::  pBuff
    double precision, dimension(nDt+1, N)   ::  u, pTraj


    ! I
    ! Filtration
    call specFilt(p0, N, Ntrc)
    pBuff(1,:)=p0
    pTraj(1,:)=p0
    
    !premier pas avec Euler-avant
    ! E1
    pBuff(2,:)=pBuff(1,:)+dt*kdvTLMPseudoSpec(N, Ntrc, L, u(1,:), &
                                        pBuff(1,:), alph, beta, gamm, rho)

    pTraj(2,:)=pBuff(2,:)
    tReal=tReal+dt

    !pas subsequents avec Leapfrog et trapezoidal
    ! Qj , j=2, nDt
    do j=2,nDt
        ! Pj
        pBuff(3,:)=rhoCenteredImplicit(N, Ntrc, dt, pBuff(1,:), rho) &
                   +2D0*dt*kdvTLMPseudoSpec(N, Ntrc, L, u(j,:), pBuff(2,:),&
                                            alph, beta, gamm)
        
        ! S
        pBuff(1,:)=pBuff(2,:)
        pBuff(2,:)=pBuff(3,:)


        pTraj(j+1,:)=pBuff(3,:)
        tReal=tReal+dt
    end do

    ! R
    !pf=pBuff(3,:)
end function kdvTLMPropagator


!-------------------------------------------------------------------!

function kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, tReal, u, pf, &
                            alph, beta, gamm, rho) &
                            result(aTraj)

    intent(in)              ::  N, Ntrc, L, dt, nDt, u, &
                                alph, beta, gamm, rho
    !intent(inout)           ::  pf, tReal

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, j
    
    double precision, dimension(N)          ::  alph, beta, gamm, rho, pf

    double precision, dimension(3, N)       ::  aBuff
    double precision, dimension(nDt+1, N)   ::  u, aTraj
    ! Initialisationuu
    ! R*
    aBuff(3,:)=pf
    aBuff(2,:)=0.0D0
    aBuff(1,:)=0.0D0

    aTraj(nDt+1,:)=pf

    ! Leapfrog
    ! Qj*, j=nDt,2
    do j=nDt, 2, -1

        ! S*
        aBuff(3,:)=aBuff(3,:)+aBuff(2,:)
        aBuff(2,:)=aBuff(1,:)
        aBuff(1,:)=0D0 !(superflus?)

        ! Pj*
        aBuff(1,:)=aBuff(1,:)&
                    +rhoCenteredImplicit(N, Ntrc, dt, aBuff(3,:), rho)
        aBuff(2,:)=aBuff(2,:)&
                    +2D0*dt*kdvTLMPseudoSpecAdj(N, Ntrc, L, u(j,:),&
                                                aBuff(3,:),alph, beta, gamm)
        aBuff(3,:)=0D0


        aTraj(j,:)=aBuff(1,:)
        tReal=tReal+dt
    end do

    ! E1*                 (vaut 0?)

    aBuff(1,:)=aBuff(1,:)+aBuff(2,:)+dt&
                          *kdvTLMPseudoSpecAdj(N, Ntrc, L, u(1,:),&
                                    aBuff(2,:), alph, beta, gamm, rho)
    aBuff(2,:)=0D0

    aTraj(1,:)=aBuff(1,:)
    tReal=tReal+dt

    ! I*
    !adj=aBuff(1,:)
    ! Filtration *
    call specFilt(aBuff(1,:), N, Ntrc)


end function kdvTLMPropagatorAdj

!-------------------------------------------------------------------!

!function kdvTLMSingularOp(N, Ntrc, L, dt, nDt, tReal, u, p0, &
!                            alph, beta, gamm, rho) &
!                            result(v)
!
!    intent(in)              ::  N, Ntrc, L, dt, nDt, u, &
!                                alph, beta, gamm, rho
!    !intent(inout)           ::  p0, tReal
!
!    double precision        ::  L, dt, tReal
!    integer                 ::  N, Ntrc, nDt, j
!    
!    double precision, dimension(N)          ::  alph, beta, gamm, rho, p0
!
!    double precision, dimension(3, N)       ::  pBuff
!    double precision, dimension(nDt+1, N)   ::  
!
!    ! no need for fullTraj !
!
!    kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, tReal, u, pf, &
!                            alph, beta, gamm, rho)
!
!end function kdvTLMSingularOp

!-------------------------------------------------------------------!

function kdvTLMPseudoSpec(N, Ntrc, L, u, p, alph, beta, gamm, rho)

    intent(in)                      ::  N, Ntrc, L, u, p, &
                                        alph, beta, gamm, rho
    optional                        ::  rho
    integer                         ::  N, Ntrc, j
    double precision                ::  L
    double precision, dimension(N)  ::  kdvTLMPseudoSpec, p,&
                                        dp, u, du, udp, pdu, d3p, &
                                        alph, beta, gamm, rho

    du=specDiff(u, 1, N, Ntrc, L)

    ! Differentiation
    ! Ds
    dp=specDiff(p, 1, N, Ntrc, L)
    d3p=specDiff(p, 3, N, Ntrc, L)
    ! Interior product
    ! P
    do j=1, N
        udp(j)=u(j)*dp(j)
        pdu(j)=p(j)*du(j)
    end do
    ! Filtration
    ! F
    call specFilt(udp, N, Ntrc)
    call specFilt(pdu, N, Ntrc)

    ! R
    kdvTLMPseudoSpec =-alph*dp-beta*(udp+pdu)-gamm*d3p
    if (present(rho)) then
        kdvTLMPseudoSpec= kdvTLMPseudoSpec - rho*dp 
    end if

    call specFilt(kdvTLMPseudoSpec, N, Ntrc)
end function kdvTLMPseudoSpec

!-------------------------------------------------------------------!

function kdvTLMPseudoSpecAdj(N, Ntrc, L, u, p, alph, beta, gamm, rho)

    intent(in)                      ::  N, Ntrc, L, u, p, &
                                        alph, beta, gamm, rho
    optional                        ::  rho
    integer                         ::  N, Ntrc, j
    double precision                ::  L
    double precision, dimension(N)  ::  kdvTLMPseudoSpecAdj, p, &
                                        dp, u, du, udp, pdu, d3p, &
                                        alph, beta, gamm, rho


    du=specDiff(u, 1, N, Ntrc, L)

    ! R*
    d3p=-gamm*p
    pdu=-beta*p
    udp=-beta*p
    dp=-alph*p
    if (present(rho)) then
        dp= dp -rho*p
    end if

    ! F* (auto-adjoint)
    call specFilt(udp, N, Ntrc)
    call specFilt(pdu, N, Ntrc)
    
    ! P*
    do j=1,N
        dp(j)=dp(j)+u(j)*udp(j)
        kdvTLMPseudoSpecAdj(j)=du(j)*pdu(j)
    end do
    
    ! Ds*
    dp=specDiffAdj(dp,1, N, Ntrc, L)
    d3p=specDiffAdj(d3p, 3, N, Ntrc, L)
    kdvTLMPseudoSpecAdj=kdvTLMPseudoSpecAdj+dp+d3p

end function kdvTLMPseudoSpecAdj


!-------------------------------------------------------------------!

! function testKdvTLMPseudoSpecAdj()

end module kdvTLMProp
