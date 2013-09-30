program kdvTestAdj

use kdvProp
use kdvTLMProp
use kdvTLMTest
implicit none


integer                 ::  N, Ntrc, nDt, i
double precision        ::  L, pAmp, diff, dt, tReal, rhoAmp
logical                 ::  test, rhoZero, forcZero, rhoCte

double precision, dimension(:), allocatable ::  alph, beta, gamm, rho, &
                                                ic, forc
double precision, dimension(:, :), allocatable  ::  xBuff, yBuff
double precision, dimension(:, :), allocatable  ::  u

Ntrc=50
N=3*Ntrc+1
L=3.D2


allocate(xBuff(3,N), yBuff(3, N))
allocate(ic(N), alph(N), beta(N), gamm(N), rho(N), forc(N))

pAmp=1.0D-1
rhoAmp=-1.0D-1

dt=1.0D-2
nDt=50

rhoZero=.False.
rhoCte=.True.
forcZero=.False.

allocate(u(nDt+1, N))


test=.true.

! Generating random fields
!   unfiltered state vectors
ic=initRandVec(N)
do i=1,3
    xBuff(i,:)=pAmp*initRandVec(N)
    yBuff(i,:)=pAmp*initRandVec(N)
end do
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
        rho=rhoAmp*initRandVec(N, Ntrc)
    end if
end if 
if (forcZero) then
    do i=1,N
        forc(i)=0D0
    end do
else
    forc=pAmp*initRandVec(N, Ntrc)
end if


u=kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, &
                            alph, beta, gamm, rho, forc)

print *, 
print *, '============================================================='
print *, '====| Adjoint test : validity compared with TLM |============'
print *, '============================================================='
print *, 
write(*, "(9A)"),  "   ic    ", " xBuff(1)", " yBuff(1)", "   alph  ",&
                   "   beta  ", "   gamm  ", "    rho  ", "   forc  "
do i=1, N
    write(*, "(8(F9.4 ))"), ic(i), xBuff(1,i), yBuff(1,i), alph(i), &
                            beta(i),gamm(i), rho(i), forc(i)
end do
print *, 
print *, 'Testing adjoint property with filtered noise initialisation'
print *, '  |<x, Ly>-<L*x, y>| <= 1D-14 ?'

if (test) then
    print *, 
    print *, '============================================================='
    print *, 'Testing autoadjoint of specFilt'
    print *, N, Ntrc, L
    if (testAutoadjointSpecFilt(N, Ntrc, L, diff, &
                                xBuff(1,:), yBuff(1,:))) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
end if

if (test) then
    print *, 
    print *, 
    print *, '============================================================='
    print *, '============================================================='
    print *, 'Testing adjoint validity of kdvTLMPseudoSpec'
    print *, N, Ntrc, L, pAmp
    if (testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff, &
                        u(1,:), xBuff(1,:), yBuff(1,:), &
                        alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
    
end if



if (test) then
    print *, 
    print *, 
    print *, '============================================================='
    print *, '============================================================='
    print *, 'Testing adjoint validity of opE1'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpE1Adj(N, Ntrc, L, dt, pAmp, diff, &
                        u(1,:), xBuff, yBuff, alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
    
    
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of opPn'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpPnAdj(N, Ntrc, L, dt, pAmp, diff, &
                        u(2,:), xBuff, yBuff, alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
    
    
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of opS'
    print *, N
    if (testOpSAdj(N, diff, xBuff, yBuff)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
end if 


if (test) then
    print *, 
    print *, 
    print *, '============================================================='
    print *, '============================================================='
    print *, 'Testing adjoint validity of S.Pn'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpSPnAdj(N, Ntrc, L, dt, pAmp, diff, &
                        u(2,:), xBuff, yBuff, alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
    
    
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of Pn.E1'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpPnE1Adj(N, Ntrc, L, dt, pAmp, diff, &
                       u(1,:), u(2,:), xBuff, yBuff,&
                       alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
end if 

if (test) then
    print *, 
    print *, 
    print *, '============================================================='
    print *, '============================================================='
    print *, 'Testing adjoint validity of S.Pn.E1'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpSPnE1Adj(N, Ntrc, L, dt, pAmp, diff, &
                        u(1,:), u(2,:), xBuff, yBuff, &
                        alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
end if


if (test) then
    print *, 
    print *, '============================================================='
    print *, '============================================================='
    print *, 'Testing adjoint validity of the sequence R.S.Pn.E1.I.F'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpAllAdj(N, Ntrc, L, dt, pAmp, diff, &
                        u(1,:),u(2,:), xBuff, yBuff, &
                        alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
end if


if (test) then
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of kdvTLMPropagator'
    print *, N, Ntrc, L, dt, nDt, pAmp
    if (testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, pAmp, diff, &
                        u, xBuff(1,:), yBuff(1,:), &
                        alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
end if
end program kdvTestAdj
