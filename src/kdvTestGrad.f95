program kdvTestGrad

use kdvProp
use kdvTLMProp
use kdvTLMTest
implicit none

integer                 ::  N, Ntrc, nDt, maxPower, i, NNDt
double precision        ::  L, pAmp, diff, dt, tReal
logical                 ::  test, rhoZero, forcZero, rhoCte

double precision, dimension(:), allocatable ::  alph, beta, gamm, rho, &
                                                ic, forc
double precision, dimension(:), allocatable  :: xBuff, yBuff


integer, dimension(:), allocatable    ::  nDtVec
double precision, dimension(:, :), allocatable  ::  u

NNDt=5
Ntrc=50
N=3*Ntrc+1
L=3.D2
rhoZero=.True.
rhoCte=.True.
forcZero=.False.


allocate(nDtVec(NNDt))
allocate(xBuff(N), yBuff(N))
allocate(ic(N), alph(N), beta(N), gamm(N), rho(N), forc(N))

pAmp=1.D-1
dt=1.D-2


maxPower=-9
nDtVec=(/2,10,20,50,100/)

! Generating random fields
!   unfiltered state vectors
call init_random_seed()
ic=initRandVec(N)
xBuff=pAmp*initRandVec(N)
yBuff=pAmp*initRandVec(N)
!   filtered parameters vectors
alph=initRandVec(N, Ntrc)
beta=initRandVec(N, Ntrc)
gamm=initRandVec(N, Ntrc)
if (rhoZero) then
    do i=1,N
        rho(i)=0D0
    end do
else
    if (rhoCte) then
        do i=1,N
            rho(i)=1.D-1
        end do
    else
        rho=initRandVec(N, Ntrc)
    end if
end if 
if (forcZero) then
    do i=1,N
        forc(i)=0D0
    end do
else
    forc=pAmp*initRandVec(N, Ntrc)
end if



print *, 
print *, '============================================================='
print *, '====| Gradient test : validity of TLM adjoint |=============='
print *, '====|                 compared with NL model  |=============='
print *, '============================================================='
print *, 
write(*, "(9A)"),  "   ic    ", "  xBuff  ", "  yBuff  ", "   alph  ",&
                   "   beta  ", "   gamm  ", "    rho  ", "   forc  "
do i=1, N
    write(*, "(8(F9.4 ))"), ic(i), xBuff(i), yBuff(i), alph(i), beta(i), &
                        gamm(i), rho(i), forc(i)
end do
do i=1,NNDt
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
