program kdvTestGrad

use kdvProp
use kdvTLMProp
use kdvTLMTest
implicit none

character*64            ::  nmlFile

integer                 ::  fun, N, Ntrc, nDt, nDtParam, maxPower,&
                            i, NNDt, NtrcParam, nuN
double precision        ::  L, pAmp, diff, dt, tReal, rhoAmp, nu, nuAmp
logical                 ::  test, rhoZero, forcZero, rhoCte, paramCte, & 
                            nuZero, outputVec

double precision, dimension(:), allocatable     ::  ic
double precision, dimension(:), allocatable  :: xBuff, yBuff


integer, dimension(:), allocatable    ::  nDtVec
double precision, dimension(:,:), allocatable   ::  alph, beta, gamm,&
                                                    rho, forc
double precision, dimension(:, :), allocatable  ::  u

namelist /GridBloc/ Ntrc, L
namelist /ParamBloc/ dt, pAmp, rhoAmp, nuAmp, nuN, nDtParam, NtrcParam
namelist /ConfigBloc/ paramCte, rhoZero, rhoCte, forcZero, nuZero, &
                      outputVec
namelist /TestGradBloc/ NNDt, maxPower, nDtVec 

nmlFile='kdvTest.nml'
open(fun, file=nmlFile)
    read(fun, nml=GridBloc)
    read(fun, nml=ParamBloc)
    read(fun, nml=ConfigBloc)
    read(fun, nml=TestGradBloc)
close(fun)



N=3*Ntrc+1
allocate(nDtVec(NNDt))
allocate(xBuff(N), yBuff(N))
allocate(ic(N))
allocate(alph(nDtParam+1,N), beta(nDtParam+1,N), gamm(nDtParam+1,N), &
            rho(nDtParam+1,N), forc(nDtParam+1,N))

nDtVec=(/1, 2,10,50,100, 500, 1000/)
! Generating random fields
!   unfiltered state vectors
!call init_random_seed()
ic=initRandVec(N)
xBuff=pAmp*initRandVec(N)
yBuff=pAmp*initRandVec(N)
!   filtered parameters vectors
if (paramCte) then
    alph(1,:)=0.0D0
    beta(1,:)=1.0D0
    gamm(1,:)=-1.0D0
else
    alph(1,:)=initRandVec(N, NtrcParam)
    beta(1,:)=initRandVec(N, NtrcParam)
    gamm(1,:)=initRandVec(N, NtrcParam)
end if
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
        rho(1,:)=rhoAmp*initRandVec(N, NtrcParam)
    end if
end if 
if (forcZero) then
    do i=1,N
        forc(1,i)=0D0
    end do
else
    forc(1,:)=pAmp*initRandVec(N, Ntrc)
end if
if (nuZero) then
    nu=0.0D0
else
    nu=nuAmp
end if


print *, 
print *, '============================================================='
print *, '====| Gradient test : validity of TLM adjoint |=============='
print *, '====|                 compared with NL model  |=============='
print *, '============================================================='
print *, 
if (outputVec) then
    write(*, "(9A)"),  "   ic    ", "  xBuff  ", "  yBuff  ", "   alph  ",&
                       "   beta  ", "   gamm  ", "    rho  ", "   forc  "
    do i=1, N
        write(*, "(8(F9.4 ))"), ic(i), xBuff(i), yBuff(i), alph(1,i), &
                            beta(1,i), gamm(1,i), rho(1,i), forc(1,i)
    end do
end if
if (.not. nuZero) then
    print *, ' nuN=',nuN
    print *, ' nu=',nu
end if
do i=1,NNDt
    nDt=nDtVec(i)
    print *, '============================================================='
    print *, '============================================================='
    print *, 'nDt=',nDt
    allocate(u(nDt+1, N))
    u=kdvPropagator(N, Ntrc, L, dt, nDt, nDtParam, tReal, ic, &
                            alph, beta, gamm, rho, nu, nuN, forc)

    print *, '============================================================='
    print *, 'Testing adjoint validity of kdvTLMPropagator'
    if (testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, nDtParam,&
                        pAmp, diff, u, xBuff, yBuff, &
                        alph, beta, gamm, rho, nu, nuN)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if

    print *, 
    print *, '============================================================='
    print *, 'Gradient test'
    call NLTestGradient(N, Ntrc, L, dt, nDt, nDtParam, maxPower,&
                        ic, alph, beta, gamm, rho, nu, nuN, forc)
    deallocate(u)
end do
end program kdvTestGrad
