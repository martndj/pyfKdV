program kdvTestAdj

use kdvTLMProp
implicit none

integer                 ::  N, Ntrc, nDt
double precision        ::  L, pAmp, diff, dt



Ntrc=150
N=3*Ntrc+1
L=3D2
pAmp=1D-2

print *, 'Testing adjoint validity of kdvTLMPseudoSpec and adjoint'
print *, '  |<x, Ly>-<L*x, y>| <= 1D-14 ?|'
if (testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test failed', diff
end if


print *, 'Testing adjoint validity of kdvTLMPropagator and adjoint'
print *, '  |<x, Ly>-<L*x, y>| <= 1D-14 ?|'
dt=0.01
nDt=2

if (testKdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test failed', diff
end if


end program kdvTestAdj
