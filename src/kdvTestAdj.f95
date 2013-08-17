program kdvTestAdj

use kdvTLMProp
implicit none

integer                 ::  N, Ntrc, nDt
double precision        ::  L, pAmp, diff, dt
logical                 ::  test


Ntrc=150
N=3*Ntrc+1
L=3D2

pAmp=1D-1
dt=0.01
nDt=11

test=.true.

print *, 
print *, 'Testing autoadjoint of specFilt'
print *, '  |<x, Ly>-<Lx, y>| <= 1D-14 ?|'
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
print *, '  |<x, Ly>-<L*x, y>| <= 1D-14 ?|'
if (testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if


if (.not.test) stop
print *, 
print *, 'Testing autoadjoint of opE1'
print *, N, Ntrc, L, dt, pAmp
print *, '  |<x, Ly>-<Lx, y>| <= 1D-14 ?|'
if (testOpE1Adj(N, Ntrc, L, dt, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if


print *, 
print *, 'Testing autoadjoint of opPn'
print *, N, Ntrc, L, dt, pAmp
print *, '  |<x, Ly>-<Lx, y>| <= 1D-14 ?|'
if (testOpPnAdj(N, Ntrc, L, dt, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if


print *, 
print *, 'Testing autoadjoint of opS'
print *, N
print *, '  |<x, Ly>-<Lx, y>| <= 1D-14 ?|'
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
print *, '  |<x, Ly>-<L*x, y>| <= 1D-14 ?|'
if (testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test FAILED', diff
    test=.false.
end if


end program kdvTestAdj
