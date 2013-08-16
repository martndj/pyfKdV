program kdvTestAdj

use kdvTLMProp
implicit none

integer                 ::  N, Ntrc
double precision        ::  L, pAmp, diff



Ntrc=150
N=3*Ntrc+1
L=3D2
pAmp=1D-2

print *, 'Testing adjoint validity of kdvTLMPseudoSpec and kdvTLMPseudoSpecAdj'
print *, '  |<x, Ly>-<L*x, y>| <= 1D-14 ?|'
if (testKdvTLMPseudoSpecAdj(N, Ntrc, L, pAmp, diff)) then 
    print *, 'Test succeeded:', diff
else
    print *, 'Test failed', diff
end if

end program kdvTestAdj
