module kdvTLMTest

use spectral
use kdvProp
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
    r=0.0D0
    do i=1,nVec
        r=r+scalar_product(a(i,:),b(i,:))
    end do
end function scalarNVec

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

subroutine init_random_seed()
    !   Reference:
    !   http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(count)
       if (count /= 0) then
          t = transfer(count, t)
       else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
          t = transfer(tms, t)
       end if
       s = ieor(t(1), t(2))
       pid = getpid() + 1099279 ! Add a prime
       s = ieor(s, pid)
       if (n >= 3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
             seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
       else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       end if
    end if
    call random_seed(put=seed)
end subroutine init_random_seed


!-------------------------------------------------------------------!

function centeredRand()
    double precision        ::  centeredRand, r
    call init_random_seed()
    call random_number(r)
    centeredRand=r-5D-1
end function centeredRand

!-------------------------------------------------------------------!

function initRandVec(N, Ntrc)
    intent(in)                      ::  N
    optional                        ::  Ntrc
    double precision, dimension(N)  ::  initRandVec
    integer                         ::  j, N, Ntrc

    do j=1,N
        initRandVec(j)=centeredRand()
    end do
    ! explicit filtering
    if (present(Ntrc)) then 
        call specFilt(initRandVec, N, Ntrc)
    end if
end function

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function testAutoAdjointSpecFilt(N, Ntrc, L, diff, x, y)
    intent(in)                      ::  N, Ntrc, L, x, y
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  x, y, Ly, Lx
    double precision                ::  L, diff
    logical                         ::  testAutoadjointSpecFilt

    double precision, parameter     :: tolerance=1D-14

    Lx=x
    Ly=y
    call specFilt(Lx, N, Ntrc)
    call specFilt(Ly, N, Ntrc)


    ! adjoint validity test
    print *, '<x,Ly>=', scalar_product(x,Ly)
    print *, '<Lx,y>=', scalar_product(Lx,y)
    diff=dabs(scalar_product(x,Ly)-scalar_product(Lx,y))
    testAutoadjointSpecFilt=(diff .le. tolerance)

end function testAutoAdjointSpecFilt
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
function testRhoCenteredImplicitAdj(N, Ntrc, L, dt, diff, x, y, rho)
    intent(in)                      ::  N, Ntrc, L, dt, x, y, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  rho, x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt
    logical                         ::  testRhoCenteredImplicitAdj

    double precision, parameter     :: tolerance=1D-14


    Ly=rhoCenteredImplicit(N, Ntrc, dt, y,rho)
    !LAdj_x=rhoCenteredImplicitAdj(N, Ntrc, dt, x, rho)
    LAdj_x=rhoCenteredImplicit(N, Ntrc, dt, x,rho)
    !  AUTOADJOINT?

    ! adjoint validity test
    print *, '<x,Ly>= ', scalar_product(x,Ly)
    print *, '<L*x,y>=', scalar_product(LAdj_x,y)
    diff=dabs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testRhoCenteredImplicitAdj=(diff .le. tolerance)

end function testRhoCenteredImplicitAdj

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function testOpE1Adj(N, Ntrc, L, dt, pAmp, diff, &
                        u, x, y, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp, &
                                        u, x, y, alph, beta, gamm, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpE1Adj

    double precision, parameter     :: tolerance=1D-14

    
    Ly=opE1(N, Ntrc, L, dt, u, y, alph, beta, gamm, rho)
    LAdj_x=opE1Adj(N, Ntrc, L, dt, u, x, alph, beta, gamm, rho)

    ! adjoint validity test
    print *, '<x,Ly>= ',scalarNVec(x,Ly,3, N)
    print *, '<L*x,y>=',scalarNVec(LAdj_x,y,3, N)
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpE1Adj=(diff .le. tolerance)

end function testOpE1Adj


!-------------------------------------------------------------------!

function testOpPnAdj(N, Ntrc, L, dt, pAmp, diff, &
                        u, x, y, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp, &
                                        u, x, y, alph, beta, gamm, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpPnAdj

    double precision, parameter     :: tolerance=1D-14


    Ly=opPn(N, Ntrc, L, dt, u, y, alph, beta, gamm, rho)
    LAdj_x=opPnAdj(N, Ntrc, L, dt, u, x, alph, beta, gamm, rho)

    ! adjoint validity test
    print *, '<x,Ly>= ',scalarNVec(x,Ly,3, N)
    print *, '<L*x,y>=',scalarNVec(LAdj_x,y,3, N)
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpPnAdj=(diff .le. tolerance)

end function testOpPnAdj


!-------------------------------------------------------------------!

function testOpSAdj(N,  diff, x, y)
    intent(in)                      ::  N, x, y
    intent(out)                     ::  diff

    integer                         ::  N, j, i
    double precision                ::  diff
    
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x
    logical                         ::  testOpSAdj

    double precision, parameter     :: tolerance=1D-14
    
    Ly=opS(N, y)
    LAdj_x=opSAdj(N, x)

    ! adjoint validity test
    print *, '<x,Ly>= ',scalarNVec(x,Ly,3, N)
    print *, '<L*x,y>=',scalarNVec(LAdj_x,y,3, N)
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpSAdj=(diff .le. tolerance)

end function testOpSAdj

!-------------------------------------------------------------------!

function testOpSPnAdj(N, Ntrc, L, dt, pAmp, diff, &
                        u, x, y, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp, &
                                        u, x, y, alph, beta, gamm, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  u, alph, beta, gamm, rho
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x, &
                                        Sy, PnAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpSPnAdj

    double precision, parameter     :: tolerance=1D-14


    Sy=opS(N, y)
    Ly=opPn(N, Ntrc, L, dt, u, Sy, alph, beta, gamm, rho)
    
    
    PnAdj_x=opPnAdj(N, Ntrc, L, dt, u, x, alph, beta, gamm, rho)
    LAdj_x=opSAdj(N, PnAdj_x)
    print *, '<x,PnSy>=    ',scalarNVec(x,Ly,3, N)
    print *, '<Pn*x,Sy>=   ',scalarNVec(PnAdj_x,Sy,3, N)
    print *, '<S*Pn*x,y>=  ',scalarNVec(LAdj_x,y,3, N)



    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpSPnAdj=(diff .le. tolerance)

end function testOpSPnAdj


!-------------------------------------------------------------------!


function testOpPnE1Adj(N, Ntrc, L, dt, pAmp, diff, &
                        u1, u2, x, y, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp, &
                                        u1, u2, x, y, alph, beta, gamm, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho
    double precision, dimension(N)  ::  u1, u2
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x, &
                                        E1y, PnAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpPnE1Adj

    double precision, parameter     :: tolerance=1D-14

    
    E1y=opE1(N, Ntrc, L, dt, u1, y,  alph, beta, gamm, rho)
    Ly=opPn(N, Ntrc, L, dt, u2, E1y, alph, beta, gamm, rho)

    PnAdj_x=opPnAdj(N, Ntrc, L, dt, u2, x, alph, beta, gamm, rho)
    LAdj_x=opE1Adj(N, Ntrc, L, dt, u1, PnAdj_x, alph, beta, gamm, rho)
    print *, '<x,PnE1y>=  ',scalarNVec(x,Ly,3, N)
    print *, '<Pn*x,E1y>= ',scalarNVec(PnAdj_x,E1y,3, N)
    print *, '<E1*Pn*x,y>=',scalarNVec(LAdj_x,y,3, N)

    
    ! adjoint validity test
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpPnE1Adj=(diff .le. tolerance)

end function testOpPnE1Adj

!-------------------------------------------------------------------!


function testOpSPnE1Adj(N, Ntrc, L, dt, pAmp, diff, &
                        u1, u2, x, y, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp, &
                                        u1, u2, x, y, alph, beta, gamm, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j, i
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho
    double precision, dimension(N)  ::  u1, u2
    double precision, dimension(3,N)::  x, y, Ly, LAdj_x, &
                                        E1y, PnE1y, &
                                        SAdj_x, PnAdjSAdj_x
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpSPnE1Adj

    double precision, parameter     :: tolerance=1D-14


    E1y=opE1(N, Ntrc, L, dt, u1, y,  alph, beta, gamm, rho)
    PnE1y=opPn(N, Ntrc, L, dt, u2, E1y, alph, beta, gamm, rho)
    Ly=opS(N, PnE1y)
    
    
    SAdj_x=opSAdj(N, x)
    PnAdjSAdj_x=opPnAdj(N, Ntrc, L, dt, u2, SAdj_x, alph, beta, gamm, rho)
    LAdj_x=opE1Adj(N, Ntrc, L, dt, u1, PnAdjSAdj_x, &
                    alph, beta, gamm, rho)
    print *, '<x,SPnE1y>=   ',scalarNVec(x,Ly,3, N)
    print *, '<S*x,PnE1y>=  ',scalarNVec(SAdj_x,PnE1y,3, N)
    print *, '<Pn*S*x,E1y>= ',scalarNVec(PnAdjSAdj_x,E1y,3, N)
    print *, '<E1*Pn*S*x,y>=',scalarNVec(LAdj_x,y,3, N)
    diff=dabs(scalarNVec(x,Ly,3, N)-scalarNVec(LAdj_x,y,3, N))
    testOpSPnE1Adj=(diff .le. tolerance)

end function testOpSPnE1Adj



!-------------------------------------------------------------------!
function testOpAllAdj(N, Ntrc, L, dt, pAmp, diff, &
                        u1, u2, x, y, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, pAmp, &
                                        u1, u2, x, y, alph, beta, gamm, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                        x, y, Ly, LAdj_x, &
                                        Fy, IE1PnSR_x
                                        
    double precision, dimension(3,N)::  IFy, E1IFy, PnE1IFy, SPnE1IFy, &
                                        R_x, SR_x, PnSR_x, E1PnSR_x
                                        
    double precision, dimension(N)  ::  u1, u2
    double precision                ::  L, diff, dt, pAmp
    logical                         ::  testOpAllAdj

    double precision, parameter     :: tolerance=1D-14


    ! Direct
    Fy=y
    call specFilt(Fy, N, Ntrc)
    IFy(1,:)=Fy
    IFy(2,:)=0D0
    IFy(3,:)=0D0
    E1IFy=opE1(N, Ntrc, L, dt, u1, IFy, alph, beta, gamm, rho)
    PnE1IFy=opPn(N, Ntrc, L, dt, u2, E1IFy, alph, beta, gamm, rho)
    SPnE1IFy=opS(N, PnE1IFy)
    Ly=SPnE1IFy(3,:)


    R_x(3,:)=x
    R_x(2,:)=0D0
    R_x(1,:)=0D0
    SR_x=opSAdj(N, R_x)
    PnSR_x=opPnAdj(N, Ntrc, L, dt, u2, SR_x, alph, beta, gamm, rho)
    E1PnSR_x=opE1Adj(N, Ntrc, L, dt, u1, PnSR_x, alph, beta, gamm, rho)
    IE1PnSR_x=E1PnSR_x(1,:)
    LAdj_x=IE1PnSR_x
    call specFilt(LAdj_x, N, Ntrc)
    print *, '<x,RSPnE1IFy>=     ',scalar_product(x,Ly)
    print *, '<R*x,SPnE1IFy>=    ',scalarNVec(R_x,SPnE1IFy,3, N)
    print *, '<S*R*x,PnE1IFy>=   ',scalarNVec(SR_x,PnE1IFy,3, N)
    print *, '<Pn*S*R*x,E1IFy>=  ',scalarNVec(PnSR_x,E1IFy,3, N)
    print *, '<E1*Pn*S*R*x,IFy>= ',scalarNVec(E1PnSR_x,IFy,3, N)
    print *, '<I*E1*Pn*S*R*x,Fy>=',scalar_product(IE1PnSR_x,Fy)
    print *, '<FI*E1*Pn*S*R*x,y>=',scalar_product(LAdj_x,y)

    diff=dabs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testOpAllAdj=(diff .le. tolerance)

end function testOpAllAdj



!-------------------------------------------------------------------!

function testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff, &
                        u, x, y, alph, beta, gamm, rho)

    intent(in)                      ::  N, Ntrc, L, pAmp, &
                                        u, x, y, alph, beta, gamm, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, j
    
    double precision, dimension(N)  ::  u, alph, beta, gamm, rho, &
                                        x, y, Ly, LAdj_x
    double precision                ::  L, diff, pAmp
    logical                         ::  testKdvTLMPseudoSpecAdj

    double precision, parameter     :: tolerance=1D-14


    Ly=kdvTLMPseudoSpec(N, Ntrc, L, u, y, alph, beta, gamm, rho)
    LAdj_x=kdvTLMPseudoSpecAdj(N, Ntrc, L, u, x, alph, beta, gamm, rho)

    ! adjoint validity test
    print *, '<x,Ly>= ', scalar_product(x,Ly)
    print *, '<L*x,y>=', scalar_product(LAdj_x,y)
    diff=dabs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testKdvTLMPseudoSpecAdj=(diff .le. tolerance)

end function testKdvTLMPseudoSpecAdj


!-------------------------------------------------------------------!


function testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, pAmp, diff,&
                                    u, x, y, alph, beta, gamm, rho)

    intent(in)                      ::  N, Ntrc, L, dt, nDt, pAmp, &
                                        u, x, y, alph, beta, gamm, rho
    intent(out)                     ::  diff

    integer                         ::  N, Ntrc, nDt, j, k
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                        x, y, Ly, LAdj_x
    double precision                ::  L, diff, pAmp, dt, tReal, tRealAdj
    logical                         ::  testKdvTLMPropagatorAdj
    
    double precision, dimension(nDt+1, N)  ::  u

    double precision, parameter     :: tolerance=1D-14

    tReal=0D0
    tRealAdj=0D0

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
    print *, '<x,Ly>= ', scalar_product(x,Ly)
    print *, '<L*x,y>=', scalar_product(LAdj_x,y)
    diff=dabs(scalar_product(x,Ly)-scalar_product(LAdj_x,y))
    testKdvTLMPropagatorAdj=(diff .le. tolerance)

end function testKdvTLMPropagatorAdj

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

subroutine testGradient(N, Ntrc, L, dt, nDt, maxPower, &
                        u, x, alph, beta, gamm, rho)
    !
    !   J(x-eps\grad J)-J(x)
    !   --------------------  -1 < O(eps) ?
    !     eps||\grad J||^2
    !
    !   <!> doit être valable jusqu'à la 
    !       *moitié* de la précision du type
    !------------------------------------------------------
    intent(in)                      ::  N, Ntrc, L, dt, nDt, &
                                        maxPower, &
                                        u, x, alph, beta, gamm, rho

    integer                         ::  N, Ntrc, nDt, maxPower, j, k, pow
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                        x, grad
    double precision                ::  L, res, pAmp, dt, &
                                        tRealFct, tRealGrad, &
                                        eps, Jeps, J0
    
    double precision, dimension(nDt+1, N)  ::  u

    double precision, parameter     :: tolerance=1D-14
    !logical                         ::  test, testGradient


    J0=fctCout(N, Ntrc, L, dt, nDt, tRealFct, u, x, &
                        alph, beta, gamm, rho )

    grad=gradFC(N, Ntrc, L, dt, nDt, tRealGrad, u, x, &
                        alph, beta, gamm, rho )


    print"(A E23.15)", "  J(x):         ",J0 
    print"(A E23.15)", "  |gradJ(x)|^2: ", scalar_product(grad,grad)

    print*,"------------------------------------------------------"
    
    print"(A 5X A 18X A 13X A 20X)", "eps","J(x-eps.gradJ)","res"
    
    do pow=-1,maxPower, -1
        eps=1D1**pow
        Jeps=fctCout(N, Ntrc, L, dt, nDt, tRealFct, u, x-eps*grad, &
                        alph, beta, gamm, rho )

        res=((J0-Jeps)/(eps*scalar_product(grad,grad)))

        !test=dabs(1D0-res).lt.eps
        !if ((pow>=-7) .and. (.not. test)) then
        !    testGradient=.false.
        !end if
        
        if (pow.eq.-8) then
            print*,"--------------| half type precision |-----------------"
        end if
        !print"(A I3  E23.15  F20.15 F20.15 A L1)",&
        !     "10^",pow, Jeps, res, 1D0-res, " ",test
        print"(A I3  E23.15  F20.15)",&
             "10^",pow, Jeps, res
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

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

subroutine NLTestGradient(N, Ntrc, L, dt, nDt, maxPower, &
                            x, alph, beta, gamm, rho, forc)
    !
    !   J(x-eps\grad J)-J(x)
    !   --------------------  -1 < O(eps) ?
    !     eps||\grad J||^2
    !
    !   J(x)=1/2||M(x)||^2
    !   gradJ(x)=L*M(x)
    !------------------------------------------------------
    intent(in)                      ::  N, Ntrc, L, dt, nDt, &
                                        maxPower, &
                                        x, alph, beta, gamm, rho, forc

    integer                         ::  N, Ntrc, nDt, maxPower, j, k, pow
    
    double precision, dimension(N)  ::  alph, beta, gamm, rho, forc, &
                                        x, grad
    double precision                ::  L, res, pAmp, dt, &
                                        tRealFct, tRealGrad, &
                                        eps, Jeps, J0
    
    double precision, dimension(nDt+1, N)  ::  u

    double precision, parameter     :: tolerance=1D-14
    !logical                         ::  test, testGradient


    J0=fctCout(N, Ntrc, L, dt, nDt, tRealFct, u, x, &
                        alph, beta, gamm, rho, forc )

    grad=gradFC(N, Ntrc, L, dt, nDt, tRealGrad, u, &
                        alph, beta, gamm, rho )


    print"(A E23.15)", "  J(x):         ",J0 
    print"(A E23.15)", "  |gradJ(x)|^2: ", scalar_product(grad,grad)

    print*,"------------------------------------------------------"
    
    print"(A 5X A 18X A 13X A 20X)", "eps","J(x-eps.gradJ)","res"
    
    do pow=-1,maxPower, -1
        eps=1.0D1**pow
        Jeps=fctCout(N, Ntrc, L, dt, nDt, tRealFct, u, x-eps*grad, &
                        alph, beta, gamm, rho, forc)

        res=((J0-Jeps)/(eps*scalar_product(grad,grad)))

        !test=dabs(1D0-res).lt.eps
        !if ((pow>=-7) .and. (.not. test)) then
        !    testGradient=.false.
        !end if
        
        if (pow.eq.-8) then
            print*,"--------------| half type precision |-----------------"
        end if
        !print"(A I3  E23.15  F20.15 F20.15 A L1)",&
        !     "10^",pow, Jeps, res, 1D0-res, " ",test
        print"(A I3  D23.15  D23.15)",&
             "10^",pow, Jeps, res
    end do
    contains
    !------------------------------------------------------
    function fctCout(N, Ntrc, L, dt, nDt, tReal, u, x, &
                        alph, beta, gamm, rho, forc)

        intent(in)                      ::  N, Ntrc, L, dt, nDt, x, &
                                            alph, beta, gamm, rho, forc
        intent(out)                     ::  u
        integer                         ::  N, Ntrc, nDt
    
        double precision, dimension(N)  ::  alph, beta, gamm, rho, forc,&
                                            x, Mx
        double precision                ::  fctCout, dt, L, tReal
        double precision, dimension(nDt+1,N)    ::  u
        
        u=kdvPropagator(N, Ntrc, L, dt, nDt, tReal,&
                            x, alph, beta, gamm, rho, forc)
        Mx=u(nDt+1, :)
        fctCout=scalar_product(Mx,Mx)/2.0D0
    end function fctCout

    function gradFC(N, Ntrc, L, dt, nDt, tReal, u, &
                        alph, beta, gamm, rho )
        intent(in)                      ::  N, Ntrc, L, dt, nDt, u, &
                                            alph, beta, gamm, rho
        integer                         ::  N, Ntrc, nDt
    
        double precision, dimension(N)  ::  alph, beta, gamm, rho, &
                                            x, gradFC
        double precision                ::  dt, L, tReal
        double precision, dimension(nDt+1,N)    ::  u

        x=u(nDt+1, :)
        gradFC=kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, tReal, u, x, & 
                                     alph, beta, gamm, rho) 
    end function gradFC

end subroutine NLTestGradient
end module kdvTLMTest
