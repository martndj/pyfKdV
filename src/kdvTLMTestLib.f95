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

function initRandVec(N, Ntrc)
    intent(in)                      ::  N, Ntrc
    double precision, dimension(N)  ::  initRandVec
    integer                         ::  j, N, Ntrc

    call random_seed()
    do j=1,N
        initRandVec(j)=centeredRand()
    end do
    ! explicit filtering
    call specFilt(initRandVec, N, Ntrc)
end function


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
    x=initRandVec(N, Ntrc)
    y=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)


    Ly=rhoCenteredImplicit(N, Ntrc, dt, y,rho)
    !LAdj_x=rhoCenteredImplicitAdj(N, Ntrc, dt, x, rho)
    LAdj_x=rhoCenteredImplicit(N, Ntrc, dt, x,rho)
    !  AUTOADJOINT?

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
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpE1Adj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    u=initRandVec(N, Ntrc)
    do i=1,3
        x(i,:)=pAmp*initRandVec(N, Ntrc)
        y(i,:)=pAmp*initRandVec(N, Ntrc)
    end do
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)
    
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
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpPnAdj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    u=initRandVec(N, Ntrc)
    do i=1,3
        x(i,:)=pAmp*initRandVec(N, Ntrc)
        y(i,:)=pAmp*initRandVec(N, Ntrc)
    end do
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)


    Ly=opPn(N, Ntrc, L, dt, u, y, alph, beta, gamm, rho)
    LAdj_x=opPnAdj(N, Ntrc, L, dt, u, x, alph, beta, gamm, rho)

    ! adjoint validity test
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpPnAdj=(diff .le. tolerance)

end function testOpPnAdj


!-------------------------------------------------------------------!

function testOpSAdj(N, Ntrc, diff)
    intent(in)                      ::  N, Ntrc
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    double precision                ::  diff
    
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    logical                         ::  testOpSAdj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    do i=1,3
        x(i,:)=initRandVec(N, Ntrc)
        y(i,:)=initRandVec(N, Ntrc)
    end do

    Ly=opS(N, y)
    LAdj_x=opSAdj(N, x)

    ! adjoint validity test
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpSAdj=(diff .le. tolerance)

end function testOpSAdj

!-------------------------------------------------------------------!

function testOpSPnAdj(N, Ntrc, L, dt, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  u, alph, beta, gamm, rho
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpSPnAdj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    u=initRandVec(N, Ntrc)
    do i=1,3
        x(i,:)=pAmp*initRandVec(N, Ntrc)
        y(i,:)=pAmp*initRandVec(N, Ntrc)
    end do
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)

    Ly=opS(N, Ly)
    LAdj_x=opSAdj(N, x)
    LAdj_x=opPnAdj(N, Ntrc, L, dt, u, LAdj_x, alph, beta, gamm, rho)

    ! adjoint validity test
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpSPnAdj=(diff .le. tolerance)

end function testOpSPnAdj


!-------------------------------------------------------------------!


function testOpPnE1Adj(N, Ntrc, L, dt, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho
    double precision, dimension(2,N)::  u
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpPnE1Adj

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    u(1,:)=initRandVec(N, Ntrc)
    u(2,:)=initRandVec(N, Ntrc)
    do i=1,3
        x(i,:)=pAmp*initRandVec(N, Ntrc)
        y(i,:)=pAmp*initRandVec(N, Ntrc)
    end do
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)

    
    Ly=opE1(N, Ntrc, L, dt, u(1,:), y,  alph, beta, gamm, rho)
    Ly=opPn(N, Ntrc, L, dt, u(2,:), Ly, alph, beta, gamm, rho)
    LAdj_x=opPnAdj(N, Ntrc, L, dt, u(2,:), x, alph, beta, gamm, rho)
    LAdj_x=opE1Adj(N, Ntrc, L, dt, u(1,:), LAdj_x, alph, beta, gamm, rho)

    ! adjoint validity test
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpPnE1Adj=(diff .le. tolerance)

end function testOpPnE1Adj

!-------------------------------------------------------------------!


function testOpSPnE1Adj(N, Ntrc, L, dt, pAmp, diff)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho
    double precision, dimension(2,N)::  u
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpSPnE1Adj

    double precision, parameter     :: tolerance=1D-14

    u(1,:)=initRandVec(N, Ntrc)
    u(2,:)=initRandVec(N, Ntrc)
    do i=1,3
        x(i,:)=pAmp*initRandVec(N, Ntrc)
        y(i,:)=pAmp*initRandVec(N, Ntrc)
    end do
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)


    Ly=opE1(N, Ntrc, L, dt, u(1,:), y,  alph, beta, gamm, rho)
    Ly=opPn(N, Ntrc, L, dt, u(2,:), Ly, alph, beta, gamm, rho)
    Ly=opS(N, Ly)
    LAdj_x=opSAdj(N, x)
    LAdj_x=opPnAdj(N, Ntrc, L, dt, u(2,:), LAdj_x, alph, beta, gamm, rho)
    LAdj_x=opE1Adj(N, Ntrc, L, dt, u(1,:), LAdj_x, alph, beta, gamm, rho)

    ! adjoint validity test
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpSPnE1Adj=(diff .le. tolerance)

end function testOpSPnE1Adj



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
    u(1,:)=initRandVec(N, Ntrc)
    u(2,:)=initRandVec(N, Ntrc)
    x=pAmp*initRandVec(N, Ntrc)
    y=pAmp*initRandVec(N, Ntrc)
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)


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
    u=initRandVec(N, Ntrc)
    x=pAmp*initRandVec(N, Ntrc)
    y=pAmp*initRandVec(N, Ntrc)
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)


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
    do j=1, nDt+1
        u(j,:)=initRandVec(N, Ntrc)
    end do
    x=pAmp*initRandVec(N, Ntrc)
    y=pAmp*initRandVec(N, Ntrc)
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)


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

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

subroutine testGradient(N, Ntrc, L, dt, nDt, pAmp, maxPower)
    !
    !   J(x-eps\grad J)-J(x)
    !   --------------------  -1 < O(eps) ?
    !     eps||\grad J||^2
    !
    !   <!> doit être valable jusqu'à la 
    !       *moitié* de la précision du type
    !------------------------------------------------------
    intent(in)                      ::  N, Ntrc, L, dt, nDt, pAmp, maxPower

    integer                         ::  N, Ntrc, nDt, maxPower, j, k, pow
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                        x, grad
    double precision                ::  L, res, pAmp, dt, &
                                        tRealFct, tRealGrad, &
                                        eps, Jeps, J0
    
    double precision, dimension(nDt+1, N)  ::  u

    double precision, parameter     :: tolerance=1D-14

    ! Generating random fields
    do j=1, nDt+1
        u(j,:)=initRandVec(N, Ntrc)
    end do
    x=pAmp*initRandVec(N, Ntrc)
    alph=initRandVec(N, Ntrc)
    beta=initRandVec(N, Ntrc)
    gamm=initRandVec(N, Ntrc)
    rho=initRandVec(N, Ntrc)


    J0=fctCout(N, Ntrc, L, dt, nDt, tRealFct, u, x, &
                        alph, beta, gamm, rho )

    grad=gradFC(N, Ntrc, L, dt, nDt, tRealGrad, u, x, &
                        alph, beta, gamm, rho )


    print"(A E23.15)", "  J(x):         ",J0 
    print"(A E23.15)", "  |gradJ(x)|^2: ", scalar_product(grad,grad)

    print*,"------------------------------------------------------"
    
    print"(A 5X A 18X A 13X A 20X)", "eps","J(x-eps.gradJ)","res", "1D0-res"
    do pow=-1,maxPower, -1
        eps=1D1**pow
        Jeps=fctCout(N, Ntrc, L, dt, nDt, tRealFct, u, x-eps*grad, &
                        alph, beta, gamm, rho )

        res=((J0-Jeps)/(eps*scalar_product(grad,grad)))

        if (pow.eq.-8) then
            print*,"--------------| half type precision |-----------------"
        end if
        print"(A I3  E23.15  F20.15 F20.15 A L1)",&
             "10^",pow, Jeps, res, 1D0-res, " ",(abs(1D0-res).lt.eps)
    end do
    contains
    !------------------------------------------------------
    function fctCout(N, Ntrc, L, dt, nDt, tReal, u, x, &
                        alph, beta, gamm, rho )

        intent(in)                      ::  N, Ntrc, L, dt, nDt, u, x, &
                                            alph, beta, gamm, rho
        integer                         ::  N, Ntrc, nDt
    
        double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                            x, Lx
        double precision                ::  fctCout, dt, L, tReal
        double precision, dimension(nDt+1,N)    ::  u
        
        Lx=kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal,&
                            u, x, alph, beta, gamm, rho)
        
        fctCout=scalar_product(Lx,Lx)/2D0
    end function fctCout

    function gradFC(N, Ntrc, L, dt, nDt, tReal, u, x, &
                        alph, beta, gamm, rho )
        intent(in)                      ::  N, Ntrc, L, dt, nDt, u, x, &
                                            alph, beta, gamm, rho
        integer                         ::  N, Ntrc, nDt
    
        double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                            x, gradFC
        double precision                ::  dt, L, tReal
        double precision, dimension(nDt+1,N)    ::  u


        gradFC=kdvTLMSingularOp(N, Ntrc, L, dt, nDt, tReal, u, x, &
                                alph, beta, gamm, rho) 
    end function gradFC

end subroutine testGradient

end module kdvTLMTest
