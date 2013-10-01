program kdvTestGrad

use kdvProp
use kdvTLMProp
use kdvTLMTest
implicit none

integer                 ::  N, Ntrc, nDt, maxPower, i, NNDt, NtrcRho
double precision        ::  L, pAmp, diff, dt, tReal, rhoAmp
logical                 ::  test, rhoZero, forcZero, rhoCte, testLTStep, &
                            testFullModel

double precision, dimension(:), allocatable ::  alph, beta, gamm, rho, &
                                                ic, forc
double precision, dimension(:), allocatable  :: xBuff, yBuff


integer, dimension(:), allocatable    ::  nDtVec
double precision, dimension(:, :), allocatable  ::  u

NNDt=5
Ntrc=100
N=3*Ntrc+1
L=3.D2

NtrcRho=5
rhoZero=.False.
!rhoZero=.True.
rhoCte=.False.
!rhoCte=.True.
forcZero=.False.


testLTStep=.True.
testFullModel=.True.


allocate(nDtVec(NNDt))
allocate(xBuff(N), yBuff(N))
allocate(ic(N), alph(N), beta(N), gamm(N), rho(N), forc(N))

pAmp=1.0D-1
rhoAmp=1.0D-1
dt=1.0D-2

maxPower=-9
nDtVec=(/1, 2,10,50,100/)

! Generating random fields
!   unfiltered state vectors
!call init_random_seed()
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
            rho(i)=rhoAmp
        end do
    else
        rho=rhoAmp*initRandVec(N, NtrcRho)
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

if (testLTStep) then
    print *, 
    print *, '============================================================='
    print *, '====| Step functions test |=================================='
    print *, '============================================================='
    print *, 
    call LTStepTestGradient(N, Ntrc, L, dt, maxPower,&
                            xBuff, 10, alph, beta, gamm, rho, forc)
end if

if (testFullModel) then
    print *, 
    print *, '============================================================='
    print *, '====| Full model test |======================================'
    print *, '============================================================='
    print *, 
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
end if
end program kdvTestGrad
