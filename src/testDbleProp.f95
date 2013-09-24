program testDbleProp

use kdvProp
use kdvTLMTest
implicit none

integer                 ::  N, Ntrc, nDt, i
double precision        ::  L, pAmp, diff, dt, tReal

double precision, dimension(:), allocatable ::  alph, beta, gamm, rho, &
                                                ic, forc, pert, &
                                                tlmF, tlmFPert


double precision, dimension(:, :), allocatable  ::  u, u_pert

Ntrc=6
N=3*Ntrc+1
L=3.D2

allocate(ic(N), alph(N), beta(N), gamm(N), rho(N), forc(N), pert(N),&
            tlmF(N), tlmFPert(N))

pAmp=1D-14
dt=1.D-2
nDt=1


! Generating random fields
!   unfiltered state vectors
!   filtered parameters vectors
alph=initRandVec(N, Ntrc)
beta=initRandVec(N, Ntrc)
gamm=initRandVec(N, Ntrc)
rho=initRandVec(N, Ntrc)
do i=1,N
    ic(i)=0D0
    forc(i)=0D0
    pert(i)=pAmp
end do

allocate(u(nDt+1,N), u_pert(nDt+1,N))
u=kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, &
                            alph, beta, gamm, rho, forc)
u_pert=kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic+pert, &
                            alph, beta, gamm, rho, forc)


tlmF=kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal, u, ic, &
                            alph, beta, gamm, rho)
tlmFPert=kdvTLMPropagator(N, Ntrc, L, dt, nDt, tReal, u, ic+pert, &
                            alph, beta, gamm, rho)

print *, "Testing numerical type of Propagator"
do i =1,N
    write(*,"(3D23.15)"), u(nDt+1,i), u_pert(nDt+1,i),&
                            u_pert(nDt+1,i)-u(nDt+1,i)
end do

print *,
print *, "Testing numerical type of TLM Propagator"
do i =1,N
    write(*,"(3D23.15)"), tlmF(i), tlmFPert(i), tlmFPert(i)-tlmF(i)
end do
end program
