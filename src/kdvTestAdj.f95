program kdvTestAdj

use kdvProp
use kdvTLMProp
use kdvTLMTest
implicit none

integer                 ::  N, Ntrc, nDt, maxPower, i
double precision        ::  L, pAmp, diff, dt, tReal
logical                 ::  test

double precision, dimension(:), allocatable ::  alph, beta, gamm, rho, &
                                                ic, forc
double precision, dimension(:, :), allocatable  ::  xBuff, yBuff
double precision, dimension(:, :), allocatable  ::  u

Ntrc=150
N=3*Ntrc+1
L=3D2


allocate(xBuff(3,N), yBuff(3, N))
allocate(ic(N), alph(N), beta(N), gamm(N), rho(N), forc(N))

pAmp=1D-1
dt=1D-1
nDt=5

allocate(u(nDt+1, N))

maxPower=-10

test=.true.

! Generating random fields
ic=initRandVec(N, Ntrc)
do i=1,3
    xBuff(i,:)=pAmp*initRandVec(N, Ntrc)
    yBuff(i,:)=pAmp*initRandVec(N, Ntrc)
end do
alph=initRandVec(N, Ntrc)
beta=initRandVec(N, Ntrc)
gamm=initRandVec(N, Ntrc)
rho=initRandVec(N, Ntrc)
forc=pAmp*initRandVec(N, Ntrc)


u=kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, &
                            alph, beta, gamm, rho, forc)

print *, 
print *, 'Testing adjoint property with filtered noise initialisation'
print *, '  |<x, Ly>-<L*x, y>| <= 1D-14 ?'

print *, 
print *, '============================================================='
print *, 'Testing autoadjoint of specFilt'
print *, N, Ntrc, L
if (testAutoadjointSpecFilt(N, Ntrc, L, diff, xBuff(1,:), yBuff(1,:))) then 
    print *, ' >>Test succeeded:', diff
else
    print *, ' >>Test FAILED', diff
    test=.false.
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
    
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of rhoCenteredImplicit'
    print *, N, Ntrc, L, dt
    if (testRhoCenteredImplicitAdj(N, Ntrc, L, dt, diff, &
                                    xBuff(1,:), yBuff(1,:), rho)) then 
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
                        u(1,:), xBuff, yBuff, alph, beta, gamm, rho)) then 
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
                        u(1:2,:), xBuff, yBuff, &
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
    print *, 'Testing adjoint validity of the sequence R.S.Pn.E1.I.F'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpAllAdj(N, Ntrc, L, dt, pAmp, diff, &
                        u(1:2,:), xBuff, yBuff, &
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
                        ic, xBuff(1,:), yBuff(1,:), &
                        alph, beta, gamm, rho)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
end if



print *, 
print *, 
print *, '============================================================='
print *, '============================================================='
print *, 'Gradient test'
    print *, N, Ntrc, L, dt, nDt, pAmp, maxPower
if (testGradient(N, Ntrc, L, dt, nDt, pAmp, maxPower,&
                                    u, xBuff(1,:), &
                                    alph, beta, gamm, rho))  then
    print *, ' >>Test succeeded!'
else
    print *, ' >>Test FAILED'
    test=.false.
end if
end program kdvTestAdj
