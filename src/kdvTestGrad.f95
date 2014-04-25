program kdvTestGrad

use kdvProp
use kdvTLMProp
use kdvTLMTest
implicit none

integer                 ::  N, Ntrc, nDt, nDtParam, maxPower,&
                            i, NNDt, NtrcRho
double precision        ::  L, pAmp, diff, dt, tReal, rhoAmp
logical                 ::  test, rhoZero, forcZero, rhoCte, & 
                            testFullModel

double precision, dimension(:), allocatable     ::  ic
double precision, dimension(:), allocatable  :: xBuff, yBuff


integer, dimension(:), allocatable    ::  nDtVec
double precision, dimension(:,:), allocatable   ::  alph, beta, gamm,&
                                                    rho, forc
double precision, dimension(:, :), allocatable  ::  u

NNDt=7
Ntrc=100
N=3*Ntrc+1
L=3.D2

NtrcRho=5
rhoZero=.False.
!rhoZero=.True.
rhoCte=.False.
!rhoCte=.True.
forcZero=.False.


testFullModel=.True.


allocate(nDtVec(NNDt))
allocate(xBuff(N), yBuff(N))
allocate(ic(N))

pAmp=1.0D-1
rhoAmp=1.0D-1
dt=1.0D-2

maxPower=-14
nDtVec=(/1, 2,10,50,100, 500, 750, 850, 1000/)
nDtParam=0
allocate(alph(nDtParam+1,N), beta(nDtParam+1,N), gamm(nDtParam+1,N), &
            rho(nDtParam+1,N), forc(nDtParam+1,N))

! Generating random fields
!   unfiltered state vectors
!call init_random_seed()
ic=initRandVec(N)
xBuff=pAmp*initRandVec(N)
yBuff=pAmp*initRandVec(N)
!   filtered parameters vectors
alph(1,:)=initRandVec(N, Ntrc)
beta(1,:)=initRandVec(N, Ntrc)
gamm(1,:)=initRandVec(N, Ntrc)
if (rhoZero) then
    do i=1,N
        rho(1,i)=0D0
    end do
else
    if (rhoCte) then
        do i=1,N
            rho(1,i)=rhoAmp
        end do
    else
        rho(1,:)=rhoAmp*initRandVec(N, NtrcRho)
    end if
end if 
if (forcZero) then
    do i=1,N
        forc(1,i)=0D0
    end do
else
    forc(1,:)=pAmp*initRandVec(N, Ntrc)
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
    write(*, "(8(F9.4 ))"), ic(i), xBuff(i), yBuff(i), alph(1,i), &
                        beta(1,i), gamm(1,i), rho(1,i), forc(1,i)
end do


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
        u=kdvPropagator(N, Ntrc, L, dt, nDt, nDtParam, tReal, ic, &
                                alph, beta, gamm, rho, forc)
    
        print *, '============================================================='
        print *, 'Testing adjoint validity of kdvTLMPropagator'
        if (testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, nDtParam,&
                            pAmp, diff, u, xBuff, yBuff, &
                            alph, beta, gamm, rho)) then 
            print *, ' >>Test succeeded:', diff
        else
            print *, ' >>Test FAILED', diff
            test=.false.
        end if
    
        print *, 
        print *, '============================================================='
        print *, 'Gradient test'
        call NLTestGradient(N, Ntrc, L, dt, nDt, nDtParam, maxPower,&
                            ic, alph, beta, gamm, rho, forc)
        deallocate(u)
    end do
end if
end program kdvTestGrad
