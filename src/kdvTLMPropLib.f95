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
    pBuff(2,:)=0D0
    pBuff(3,:)=0D0
    if (present(pTraj)) pTraj(1,:)=p0

    !premier pas avec Euler-avant
    ! E1
    pBuff=opE1(N, Ntrc, L, dt, u(1,:), pBuff, alph, beta, gamm, rho)
    if (present(pTraj)) pTraj(2,:)=pBuff(2,:)
    tReal=tReal+dt

    !pas subsequents avec Leapfrog et trapezoidal
    ! Qj , j=2, nDt
    do j=2,nDt
        ! Pj
        pBuff=opPn(N, Ntrc, L, dt, u(j,:), pBuff, alph, beta, gamm, rho)
        ! S
        pBuff=opS(N, pBuff)
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
        aBuff=opSAdj(N, aBuff)
        ! Pj*
        aBuff=opPnAdj(N, Ntrc, L, dt, u(j,:), aBuff, alph, beta, gamm, rho)
    end do

    ! E1*              
    if (present(aTraj)) aTraj(2,:)=aBuff(2,:)
    tReal=tReal+dt
    aBuff=opE1Adj(N, Ntrc, L, dt, u(1,:), aBuff, alph, beta, gamm, rho)
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


function opE1(N, Ntrc, L, dt, u, pBuff, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, u, pBuff, &
                                        alph, beta, gamm, rho
    integer                         ::  N, Ntrc
    double precision                ::  L, dt

    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  pBuff, opE1

    opE1(1,:)=pBuff(1,:)
    opE1(2,:)=pBuff(1,:)+dt*kdvTLMPseudoSpec(N, Ntrc, L, u, &
                                        pBuff(1,:), alph, beta, gamm, rho)
    opE1(3,:)=pBuff(3,:)
    
end function opE1

!-------------------------------------------------------------------!

function opE1Adj(N, Ntrc, L, dt, u, aBuff, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, u, aBuff, &
                                        alph, beta, gamm, rho
    integer                         ::  N, Ntrc
    double precision                ::  L, dt

    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  aBuff, opE1Adj

    opE1Adj(1,:)=aBuff(1,:)+aBuff(2,:) &
                 +dt*kdvTLMPseudoSpecAdj(N, Ntrc, L, u, aBuff(2,:),&
                                            alph, beta, gamm, rho)
    opE1Adj(2,:)=0D0
    opE1Adj(3,:)=aBuff(3,:)
    
end function opE1Adj

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function opPn(N, Ntrc, L, dt, u, pBuff, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, u, pBuff, &
                                        alph, beta, gamm, rho
    integer                         ::  N, Ntrc
    double precision                ::  L, dt

    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  pBuff, opPn

    opPn(1,:)=pBuff(1,:)
    opPn(2,:)=pBuff(2,:)
    opPn(3,:)=rhoCenteredImplicit(N, Ntrc, dt, pBuff(1,:), rho) &
                   +2D0*dt*kdvTLMPseudoSpec(N, Ntrc, L, u, pBuff(2,:),&
                                            alph, beta, gamm)
    
end function opPn

!-------------------------------------------------------------------!

function opPnAdj(N, Ntrc, L, dt, u, aBuff, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, u, aBuff, &
                                        alph, beta, gamm, rho
    integer                         ::  N, Ntrc
    double precision                ::  L, dt

    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  aBuff, opPnAdj

    opPnAdj(1,:)=aBuff(1,:) &
                 +rhoCenteredImplicit(N, Ntrc, dt, aBuff(3,:), rho) 
    opPnAdj(2,:)=aBuff(2,:) &
                 +2D0*dt*kdvTLMPseudoSpecAdj(N, Ntrc, L, u, aBuff(3,:),&
                                            alph, beta, gamm)
    opPnAdj(3,:)=0D0
    
end function opPnAdj

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function opS(N, pBuff)
    intent(in)                      ::  N, pBuff
    integer                         ::  N
    double precision, dimension(3,N)::  pBuff, opS

    opS(1,:)=pBuff(2,:)
    opS(2,:)=pBuff(3,:)
    opS(3,:)=pBuff(3,:)
    
end function opS

!-------------------------------------------------------------------!

function opSAdj(N, aBuff)
    intent(in)                      ::  N, aBuff
    integer                         ::  N
    double precision, dimension(3,N)::  aBuff, opSAdj

    opSAdj(1,:)=0D0
    opSAdj(2,:)=aBuff(1,:)
    opSAdj(3,:)=aBuff(2,:)+aBuff(3,:)
    
end function opSAdj

!-------------------------------------------------------------------!
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

!--------------------------------------------------------------------

function rhoCenteredImplicitAdj(N, Ntrc, dt, u, rho)
    intent (in)                     ::  N, Ntrc, dt, u, rho
    integer                         ::  N, Ntrc
    double precision                ::  dt
    double precision, dimension(N)  ::  u, rho, rhoCenteredImplicitAdj

    rhoCenteredImplicitAdj=specFiltCopy(u, N, Ntrc)
    rhoCenteredImplicitAdj=(1D0-dt*rho)/(1D0+dt*rho)*rhoCenteredImplicitAdj
end function rhoCenteredImplicitAdj




!----| Adjoint diagnostics |----------------------------------------!

function scalarNVec(a,b, nVec, N) result(r)
    intent(in)                          :: a, b, nVec, N
    double precision, dimension(nVec,N) :: a,b
    integer                             :: i, nVec, N
    double precision                    :: r
    r=0D0
    do i=1,nVec
        r=r+scalar_product(a(i,:),b(i,:))
    end do
end function scalarNVec

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function testRhoCenteredImplicitAdj(N, Ntrc, L, dt, diff)
    intent(in)                      ::  N, Ntrc, L, dt
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  rho, x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt
    logical                         ::  testRhoCenteredImplicitAdj

    integer, parameter              :: seed=816322
    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call srand(seed)
    do j=1, N
        x(j)=(rand()-5D-1)
        y(j)=(rand()-5D-1)
        rho(j)=(rand()-5D-1)
    end do


    Ly=rhoCenteredImplicit(N, Ntrc, dt, y,rho)
    LAdj_x=rhoCenteredImplicitAdj(N, Ntrc, dt, x,rho)

    !LAdj_x=rhoCenteredImplicitAdj(N, Ntrc, dt, x,rho)
    !  NOT AUTOADJOINT!

    ! adjoint validity test
    diff=abs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testRhoCenteredImplicitAdj=(diff .le. tolerance)

end function testRhoCenteredImplicitAdj

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function testOpE1Adj(N, Ntrc, L, dt, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  u, alph, beta, gamm, rho
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpE1Adj

    integer, parameter              :: seed=816322
    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call srand(seed)
    do j=1, N
        u(j)=(rand()-5D-1)
        do i=1,3
            x(i,j)=pAmp*(rand()-5D-1)
            y(i,j)=pAmp*(rand()-5D-1)
        end do
        alph(j)=(rand()-5D-1)
        beta(j)=(rand()-5D-1)
        gamm(j)=(rand()-5D-1)
        rho(j)=(rand()-5D-1)
    end do

    ! explicit filtering
    call specFilt(u, N, Ntrc)
    do i=1,3
        call specFilt(x(i,:), N, Ntrc)
        call specFilt(y(i,:), N, Ntrc)
    end do
    call specFilt(alph, N, Ntrc)
    call specFilt(beta, N, Ntrc)
    call specFilt(gamm, N, Ntrc)
    call specFilt(rho, N, Ntrc)

    Ly=opE1(N, Ntrc, L, dt, u, y, alph, beta, gamm, rho)
    LAdj_x=opE1Adj(N, Ntrc, L, dt, u, x, alph, beta, gamm, rho)

    ! adjoint validity test
    diff=abs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpE1Adj=(diff .le. tolerance)

end function testOpE1Adj


!-------------------------------------------------------------------!

function testOpPnAdj(N, Ntrc, L, dt, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  u, alph, beta, gamm, rho
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpPnAdj

    integer, parameter              :: seed=816322
    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call srand(seed)
    do j=1, N
        u(j)=(rand()-5D-1)
        do i=1,3
            x(i,j)=pAmp*(rand()-5D-1)
            y(i,j)=pAmp*(rand()-5D-1)
        end do
        alph(j)=(rand()-5D-1)
        beta(j)=(rand()-5D-1)
        gamm(j)=(rand()-5D-1)
        rho(j)=(rand()-5D-1)
    end do

    ! explicit filtering
    call specFilt(u, N, Ntrc)
    do i=1,3
        call specFilt(x(i,:), N, Ntrc)
        call specFilt(y(i,:), N, Ntrc)
    end do
    call specFilt(alph, N, Ntrc)
    call specFilt(beta, N, Ntrc)
    call specFilt(gamm, N, Ntrc)
    call specFilt(rho, N, Ntrc)

    Ly=opPn(N, Ntrc, L, dt, u, y, alph, beta, gamm, rho)
    LAdj_x=opPnAdj(N, Ntrc, L, dt, u, x, alph, beta, gamm, rho)

    ! adjoint validity test
    diff=abs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpPnAdj=(diff .le. tolerance)

end function testOpPnAdj


!-------------------------------------------------------------------!

function testOpSAdj(N, diff)
    intent(in)                      ::  N
    intent(out)                     ::  diff

    integer                         ::  N, j, i
    double precision                ::  diff
    
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    logical                         ::  testOpSAdj

    integer, parameter              :: seed=816322
    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call srand(seed)
    do j=1, N
        do i=1,3
            x(i,j)=(rand()-5D-1)
            y(i,j)=(rand()-5D-1)
        end do
    end do

    Ly=opS(N, y)
    LAdj_x=opSAdj(N, x)

    ! adjoint validity test
    diff=abs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpSAdj=(diff .le. tolerance)

end function testOpSAdj


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
