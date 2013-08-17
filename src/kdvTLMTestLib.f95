module kdvTLMTest

use spectral
use kdvProp, only: rhoCenteredImplicit
use kdvTLMProp
use matrix, only: scalar_product
implicit none

contains
!====================================================================

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

function centeredRand()
    double precision        :: centeredRand, r
    call random_number(r)
    centeredRand=r-5D-1
end function centeredRand

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
function testRhoCenteredImplicitAdj(N, Ntrc, L, dt, diff)
    intent(in)                      ::  N, Ntrc, L, dt
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  rho, x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt
    logical                         ::  testRhoCenteredImplicitAdj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call random_seed()
    do j=1, N
        x(j)=centeredRand()
        y(j)=centeredRand()
        rho(j)=centeredRand()
    end do


    Ly=rhoCenteredImplicit(N, Ntrc, dt, y,rho)
    LAdj_x=rhoCenteredImplicitAdj(N, Ntrc, dt, x,rho)

    !LAdj_x=rhoCenteredImplicitAdj(N, Ntrc, dt, x,rho)
    !  NOT AUTOADJOINT!

    ! adjoint validity test
    diff=dabs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
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

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call random_seed()
    do j=1, N
        u(j)=centeredRand()
        do i=1,3
            x(i,j)=pAmp*centeredRand()
            y(i,j)=pAmp*centeredRand()
        end do
        alph(j)=centeredRand()
        beta(j)=centeredRand()
        gamm(j)=centeredRand()
        rho(j)=centeredRand()
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
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
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

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call random_seed()
    do j=1, N
        u(j)=centeredRand()
        do i=1,3
            x(i,j)=pAmp*centeredRand()
            y(i,j)=pAmp*centeredRand()
        end do
        alph(j)=centeredRand()
        beta(j)=centeredRand()
        gamm(j)=centeredRand()
        rho(j)=centeredRand()
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
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
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

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call random_seed()
    do j=1, N
        do i=1,3
            x(i,j)=centeredRand()
            y(i,j)=centeredRand()
        end do
    end do

    Ly=opS(N, y)
    LAdj_x=opSAdj(N, x)

    ! adjoint validity test
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpSAdj=(diff .le. tolerance)

end function testOpSAdj

!-------------------------------------------------------------------!

function testOpQnAdj(N, Ntrc, L, dt, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  u, alph, beta, gamm, rho
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpQnAdj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call random_seed()
    do j=1, N
        u(j)=centeredRand()
        do i=1,3
            x(i,j)=pAmp*centeredRand()
            y(i,j)=pAmp*centeredRand()
        end do
        alph(j)=centeredRand()
        beta(j)=centeredRand()
        gamm(j)=centeredRand()
        rho(j)=centeredRand()
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

    Ly=y
    Ly=opPn(N, Ntrc, L, dt, u, Ly, alph, beta, gamm, rho)
    Ly=opS(N, Ly)
    LAdj_x=x
    LAdj_x=opSAdj(N, LAdj_x)
    LAdj_x=opPnAdj(N, Ntrc, L, dt, u, LAdj_x, alph, beta, gamm, rho)

    ! adjoint validity test
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpQnAdj=(diff .le. tolerance)

end function testOpQnAdj


!-------------------------------------------------------------------!

function testOpAllAdj(N, Ntrc, L, dt, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                        x, y, Ly, LAdj_x
    double precision, dimension(3,N)::  buff
    double precision, dimension(2,N)::  u
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpAllAdj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call random_seed()
    do j=1, N
        u(1,j)=centeredRand()
        u(2,j)=centeredRand()
        x(j)=pAmp*centeredRand()
        y(j)=pAmp*centeredRand()
        alph(j)=centeredRand()
        beta(j)=centeredRand()
        gamm(j)=centeredRand()
        rho(j)=centeredRand()
    end do

    ! explicit filtering
    call specFilt(u(1,:), N, Ntrc)
    call specFilt(u(2,:), N, Ntrc)
    call specFilt(alph, N, Ntrc)
    call specFilt(beta, N, Ntrc)
    call specFilt(gamm, N, Ntrc)
    call specFilt(rho, N, Ntrc)


    ! Direct
    call specFilt(y, N, Ntrc)
    buff(1,:)=y
    buff(2,:)=0D0
    buff(3,:)=0D0
    buff=opE1(N, Ntrc, L, dt, u(1,:), buff, alph, beta, gamm, rho)
    buff=opPn(N, Ntrc, L, dt, u(2,:), buff, alph, beta, gamm, rho)
    buff=opS(N, buff)
    Ly=buff(3,:)

    ! Adjoint
    buff(3,:)=x
    buff(2,:)=0D0
    buff(1,:)=0D0
    buff=opSAdj(N, buff)
    buff=opPnAdj(N, Ntrc, L, dt, u(2,:), buff, alph, beta, gamm, rho)
    buff=opE1Adj(N, Ntrc, L, dt, u(1,:), buff, alph, beta, gamm, rho)
    LAdj_x=buff(1,:)
    call specFilt(LAdj_x, N, Ntrc)

    ! adjoint validity test
    diff=dabs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testOpAllAdj=(diff .le. tolerance)

end function testOpAllAdj



!-------------------------------------------------------------------!

function testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  u, alph, beta, gamm, rho, &
                                        x, y, Ly, LAdj_x
    double precision                ::  L, diff, pAmp
    logical                         ::  testKdvTLMPseudoSpecAdj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call random_seed()
    do j=1, N
        u(j)=centeredRand()
        x(j)=pAmp*centeredRand()
        y(j)=pAmp*centeredRand()
        alph(j)=centeredRand()
        beta(j)=centeredRand()
        gamm(j)=centeredRand()
        rho(j)=centeredRand()
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
    diff=dabs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
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

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    call random_seed()
    do j=1, nDt+1
        do k=1, N
            u(j,k)=centeredRand()
        end do
    end do
    do j=1, N
        x(j)=pAmp*centeredRand()
        y(j)=pAmp*centeredRand()
        alph(j)=centeredRand()
        beta(j)=centeredRand()
        gamm(j)=centeredRand()
        rho(j)=centeredRand()
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
    diff=dabs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testKdvTLMPropagatorAdj=(diff .le. tolerance)

end function testKdvTLMPropagatorAdj

end module kdvTLMTest
