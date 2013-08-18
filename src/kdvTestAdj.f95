program kdvTestAdj

!use kdvTLMProp
use kdvTLMTest
implicit none

integer                 ::  N, Ntrc, nDt, maxPower
double precision        ::  L, pAmp, diff, dt
logical                 ::  test


Ntrc=150
N=3*Ntrc+1
L=3D2

pAmp=1D-1
dt=1D-2
nDt=5

maxPower=-10

test=.true.

print *, 
print *, 'Testing adjoint property with filtered noise initialisation'
print *, '  |<x, Ly>-<Lx, y>| <= 1D-14 ?|'

print *, 
print *, '============================================================='
print *, 'Testing autoadjoint of specFilt'
print *, N, Ntrc, L
if (testAutoadjointSpecFilt(N, Ntrc, L, diff)) then 
    print *, ' >>Test succeeded:', diff
else
    print *, ' >>Test FAILED', diff
    test=.false.
end if

if (test) then
    print *, 
    print *, '============================================================='
    print *, '============================================================='
    print *, 'Testing adjoint validity of kdvTLMPseudoSpec'
    print *, N, Ntrc, L, pAmp
    if (testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
    
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of rhoCenteredImplicit'
    print *, N, Ntrc, L, dt
    if (testRhoCenteredImplicitAdj(N, Ntrc, L, dt, diff)) then 
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
    print *, 'Testing adjoint validity of opE1'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpE1Adj(N, Ntrc, L, dt, pAmp, diff)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
    
    
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of opPn'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpPnAdj(N, Ntrc, L, dt, pAmp, diff)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
    
    
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of opS'
    print *, N
    if (testOpSAdj(N, diff)) then 
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
    print *, 'Testing adjoint validity of S.Pn'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpSPnAdj(N, Ntrc, L, dt, pAmp, diff)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
    
    
    print *, 
    print *, '============================================================='
    print *, 'Testing adjoint validity of S.Pn.E1'
    print *, N, Ntrc, L, dt, pAmp
    if (testOpSPnE1Adj(N, Ntrc, L, dt, pAmp, diff)) then 
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
    if (testOpAllAdj(N, Ntrc, L, dt, pAmp, diff)) then 
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
    if (testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, pAmp, diff)) then 
        print *, ' >>Test succeeded:', diff
    else
        print *, ' >>Test FAILED', diff
        test=.false.
    end if
end if



print *, 
print *, '============================================================='
print *, 'Gradient test'
call testGradient(N, Ntrc, L, dt, nDt, pAmp, maxPower)
end program kdvTestAdj
