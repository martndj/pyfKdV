module kdvTLMProp

use spectral
use kdvProp, only: rhoCenteredImplicit
use matrix, only: scalar_product
implicit none

contains
!====================================================================


function kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal, u, p0, &
                            alph, beta, gamm, rho, pTraj) &
                            result(pf)

    intent(in)              ::  N, Ntrc, L, dt, nDt, u, &
                                alph, beta, gamm, rho
    !intent(inout)           ::  p0, tReal
    optional                ::  pTraj

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, j
    
    double precision, dimension(N)          ::  alph, beta, gamm, rho, &
                                                p0, pf

    double precision, dimension(3, N)       ::  pBuff
    double precision, dimension(nDt+1, N)   ::  u, pTraj


    ! F : Filtration
    call specFilt(p0, N, Ntrc)

    ! I : initialisation
    pBuff(1,:)=p0
    if (present(pTraj)) pTraj(1,:)=p0
    
    !premier pas avec Euler-avant
    ! E1
    pBuff(2,:)=pBuff(1,:)+dt*kdvTLMPseudoSpec(N, Ntrc, L, u(1,:), &
                                        pBuff(1,:), alph, beta, gamm, rho)

    if (present(pTraj)) pTraj(2,:)=pBuff(2,:)
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


        if (present(pTraj)) pTraj(j+1,:)=pBuff(3,:)
        tReal=tReal+dt
    end do

    ! R
    pf=pBuff(3,:)
end function kdvTLMPropagator


!-------------------------------------------------------------------!

function kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, tReal, u, pf, &
                            alph, beta, gamm, rho, aTraj) &
                            result(adj)

    intent(in)              ::  N, Ntrc, L, dt, nDt, u, &
                                alph, beta, gamm, rho
    !intent(inout)           ::  pf, tReal
    optional                ::  aTraj

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, j
    
    double precision, dimension(N)          ::  alph, beta, gamm, rho, &
                                                pf, adj

    double precision, dimension(3, N)       ::  aBuff
    double precision, dimension(nDt+1, N)   ::  u, aTraj
    
    ! R*
    aBuff(3,:)=pf
    aBuff(2,:)=0.0D0
    aBuff(1,:)=0.0D0

    tReal=0D0

    ! Leapfrog
    ! Qj*, j=nDt,2
    do j=nDt, 2, -1

        if (present(aTraj)) aTraj(j+1,:)=aBuff(3,:)
        tReal=tReal+dt

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


    end do

    ! E1*                 (vaut 0?)
    if (present(aTraj)) aTraj(2,:)=aBuff(2,:)
    tReal=tReal+dt

    aBuff(1,:)=aBuff(1,:)+aBuff(2,:) &
                         +dt*kdvTLMPseudoSpecAdj(N, Ntrc, L, u(1,:),&
                                    aBuff(2,:), alph, beta, gamm, rho)
    aBuff(2,:)=0D0


    ! I*
    adj=aBuff(1,:)

    ! F* : Filtration *
    call specFilt(adj, N, Ntrc)
    if (present(aTraj)) aTraj(1,:)=aBuff(1,:)


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



!----| Linear differential operator and adjoint |-------------------!
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


!----| Adjoint diagnostics |----------------------------------------!
!-------------------------------------------------------------------!


function testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  u, alph, beta, gamm, rho, &
                                        x, y, Ly, LAdj_x
    double precision                ::  L, diff, pAmp
    logical                         ::  testKdvTLMPseudoSpecAdj

    integer, parameter              :: seed=816322
    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call srand(seed)
    do j=1, N
        u(j)=(rand()-5D-1)
        x(j)=pAmp*(rand()-5D-1)
        y(j)=pAmp*(rand()-5D-1)
        alph(j)=(rand()-5D-1)
        beta(j)=(rand()-5D-1)
        gamm(j)=(rand()-5D-1)
        rho(j)=(rand()-5D-1)
    end do

    ! explicit filtering
    call specFilt(u, N, Ntrc)
    call specFilt(x, N, Ntrc)
    call specFilt(y, N, Ntrc)
    call specFilt(alph, N, Ntrc)
    call specFilt(beta, N, Ntrc)
    call specFilt(gamm, N, Ntrc)
    call specFilt(rho, N, Ntrc)

    Ly=kdvTLMPseudoSpec(N, Ntrc, L, u, y, alph, beta, gamm, rho)
    LAdj_x=kdvTLMPseudoSpecAdj(N, Ntrc, L, u, x, alph, beta, gamm, rho)

    ! adjoint validity test
    diff=abs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testKdvTLMPseudoSpecAdj=(diff .le. tolerance)

end function testKdvTLMPseudoSpecAdj


!-------------------------------------------------------------------!


function testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, dt, nDt, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, nDt, j, k
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                        x, y, Ly, LAdj_x
    double precision                ::  L, diff, pAmp, dt, tReal, tRealAdj
    logical                         ::  testKdvTLMPropagatorAdj
    
    double precision, dimension(nDt+1, N)  ::  u

    integer, parameter              :: seed=816322
    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call srand(seed)
    do j=1, nDt+1
        do k=1, N
            u(j,k)=(rand()-5D-1)
        end do
    end do
    do j=1, N
        x(j)=pAmp*(rand()-5D-1)
        y(j)=pAmp*(rand()-5D-1)
        alph(j)=(rand()-5D-1)
        beta(j)=(rand()-5D-1)
        gamm(j)=(rand()-5D-1)
        rho(j)=(rand()-5D-1)
    end do

    ! explicit filtering
    call specFilt(u, N, Ntrc)
    call specFilt(x, N, Ntrc)
    call specFilt(y, N, Ntrc)
    call specFilt(alph, N, Ntrc)
    call specFilt(beta, N, Ntrc)
    call specFilt(gamm, N, Ntrc)
    call specFilt(rho, N, Ntrc)

    Ly=kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal,&
                        u, y, alph, beta, gamm, rho)
    LAdj_x=kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, tRealAdj, &
                        u, x, alph, beta, gamm, rho)
    ! tReal coherence
    if (tRealAdj.ne.tReal) then 
        print *, 'testKdvTLMPropagatorAdj: tReal coherence fail', &
                    tRealAdj, tReal
        stop
    end if

    ! adjoint validity test
    diff=abs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testKdvTLMPropagatorAdj=(diff .le. tolerance)

end function testKdvTLMPropagatorAdj


end module kdvTLMProp
