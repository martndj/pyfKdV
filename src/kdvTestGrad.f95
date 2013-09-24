program kdvTestGrad

use kdvProp
use kdvTLMProp
use kdvTLMTest
implicit none

integer                 ::  N, Ntrc, nDt, maxPower, i
double precision        ::  L, pAmp, diff, dt, tReal
logical                 ::  test

double precision, dimension(:), allocatable ::  alph, beta, gamm, rho, &
                                                ic, forc
double precision, dimension(:), allocatable  :: xBuff, yBuff


integer, dimension(5)    ::  nDtVec
double precision, dimension(:, :), allocatable  ::  u

Ntrc=50
N=3*Ntrc+1
L=3.D2


allocate(xBuff(N), yBuff(N))
allocate(ic(N), alph(N), beta(N), gamm(N), rho(N), forc(N))

pAmp=1.D-1
dt=1.D-2


maxPower=-9
nDtVec=(/2,10,20,50,100/)

! Generating random fields
!   unfiltered state vectors
ic=initRandVec(N)
xBuff=pAmp*initRandVec(N)
yBuff=pAmp*initRandVec(N)
!   filtered parameters vectors
alph=initRandVec(N, Ntrc)
beta=initRandVec(N, Ntrc)
gamm=initRandVec(N, Ntrc)
rho=initRandVec(N, Ntrc)
forc=pAmp*initRandVec(N, Ntrc)



do i=1,5
    nDt=nDtVec(i)
    print *, '============================================================='
    print *, '============================================================='
    print *, 'nDt=',nDt
    allocate(u(nDt+1, N))
    u=kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, &
                            alph, beta, gamm, rho, forc)

    print *, '============================================================='
    print *, 'Testing adjoint validity of kdvTLMPropagator'
    if (testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, pAmp, diff, &
                        u, xBuff, yBuff, &
                        alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if

    print *, 
    print *, '============================================================='
    print *, 'Gradient test'
    call NLTestGradient(N, Ntrc, L, dt, nDt, maxPower,&
                                        xBuff, &
                                        alph, beta, gamm, rho, forc)
    deallocate(u)
end do
end program kdvTestGrad
