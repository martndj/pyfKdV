program kdvTestAdj

!use kdvTLMProp
use kdvTLMTest
implicit none

integer                 ::  N, Ntrc, nDt
double precision        ::  L, pAmp, diff, dt
logical                 ::  test


Ntrc=150
N=3*Ntrc+1
L=3D2

pAmp=1D-1
dt=1D-2
nDt=10

test=.true.

print *, 
print *, 'Testing adjoint property with filtered noise initialisation'
print *, '  |<x, Ly>-<Lx, y>| <= 1D-14 ?|'

print *, 
print *, 'Testing autoadjoint of specFilt'
if (testAutoadjointSpecFilt(N, Ntrc, L, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if

if (.not.test) stop
print *, 
print *, 'Testing adjoint validity of kdvTLMPseudoSpec and adjoint'
print *, N, Ntrc, L, pAmp
if (testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if

print *, 
print *, 'Testing adjoint validity of rhoCenteredImplicit and adjoint'
print *, N, Ntrc, L, dt
if (testRhoCenteredImplicitAdj(N, Ntrc, L, dt, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if



if (.not.test) stop
print *, 
print *, 'Testing autoadjoint of opE1'
print *, N, Ntrc, L, dt, pAmp
if (testOpE1Adj(N, Ntrc, L, dt, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if


print *, 
print *, 'Testing autoadjoint of opPn'
print *, N, Ntrc, L, dt, pAmp
if (testOpPnAdj(N, Ntrc, L, dt, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if


print *, 
print *, 'Testing autoadjoint of opS'
print *, N
if (testOpSAdj(N, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if


if (.not.test) stop
print *, 
print *, 'Testing adjoint validity of kdvTLMPropagator and adjoint'
print *, N, Ntrc, L, dt, nDt, pAmp
if (testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if


end program kdvTestAdj
