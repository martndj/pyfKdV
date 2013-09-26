module kdvProp

use spectral
implicit none


contains
!====================================================================


function kdvPropagator(N, Ntrc, L, dt, nDt, tReal, ic, &
                            alph, beta, gamm, rho, forc) &
                            result(traj)
    !    Integral propagator
    !
    !    <@> TODO
    !        slow start 
    !-----------------------------------------------------------
    intent(in)              ::  N, Ntrc, L, dt, nDt, &
                                alph, beta, gamm, rho, forc

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, j
    
    double precision, dimension(N)          ::  alph, beta, gamm, rho, &
                                                ic, forc
    double precision, dimension(nDt+1, N)   ::  traj

    ! explicit filtering of the IC
    call specFilt(ic, N, Ntrc)

    ! first step : Euler-forward
    traj(1,:)=ic
    traj(2,:)=traj(1,:)+dt*kdvPseudoSpec(N, Ntrc, L, traj(1,:), &
                                         alph, beta, gamm, rho, forc)

    tReal=tReal+dt

    ! subsequent step with mised Leapfrog-centered scheme
    !   leapfrog : Ax.A, Axxx, f'(x)
    !   centered : r.A
    do j=2,nDt
        traj(j+1,:)=rhoCenteredImplicit(N, Ntrc, dt, traj(j-1,:),rho)&
                    +2.0D0*dt*kdvPseudoSpec(N, Ntrc, L, traj(j,:),&
                                          alph, beta, gamm, forc=forc)
       
        tReal=tReal+dt
    end do
end function kdvPropagator

!--------------------------------------------------------------------

function rhoCenteredImplicit(N, Ntrc, dt, u, rho)
    intent (in)                     ::  N, Ntrc, dt, u, rho
    integer                         ::  N, Ntrc
    double precision                ::  dt
    double precision, dimension(N)  ::  u, rho, rhoCenteredImplicit

    rhoCenteredImplicit=(1.0D0-dt*rho)/(1.0D0+dt*rho)*u
    ! prevent aliasing from multiplication (localised parameters)
    call specFilt(rhoCenteredImplicit, N, Ntrc)
end function rhoCenteredImplicit

!--------------------------------------------------------------------!
function kdvPseudoSpec(N, Ntrc, L, u, alph, beta, gamm, rho, forc)
    !   Differential equation for Korteweg de-Vries system
    !      
    !   pseudospectral derivation on a periodic domain
    !
    !   INS
    !       N, Ntrc, L  :   grid parameters
    !       u           :   real state vector
    !       alph, beta, gamm, rho  
    !                   :   instantaneous parameters
    !       forc        :   instantaneous forcing
    !
    !   <!> rho, forc are optional   
    !
    !   RETURNS
    !       Dt(u)=-alph*Dx(u)-beta*u*Dx(u)-gamm*Dx3(u)-rho*u+forc
    !--------------------------------------------------------------
    intent(in)                      ::  N, Ntrc, L, u,&
                                        alph, beta, gamm, rho, forc
    optional                        ::  rho, forc 
    integer                         ::  N, Ntrc, j
    double precision                ::  L
    double precision, dimension(N)  ::  kdvPseudoSpec, &
                                        u, udu, du, d3u, &
                                        alph, beta, gamm, rho, forc

    ! udu=u*du
    du=specDiff(u, 1, N, Ntrc, L)    !=du
    d3u=specDiff(u, 3, N, Ntrc, L)
    do j=1,N
        udu(j)=u(j)*du(j)        !=udu
    end do
    
    ! prevent aliasing from multiplication
    call specFilt(udu, N, Ntrc)
    
    kdvPseudoSpec= - alph*du - beta*udu - gamm*d3u
    if (present(rho)) then
        kdvPseudoSpec= kdvPseudoSpec - rho*u 
    end if
    if (present(forc)) then
        kdvPseudoSpec= kdvPseudoSpec + forc
    end if
    
    ! prevent aliasing from multiplication (localised parameters)
    call specFilt(kdvPseudoSpec, N, Ntrc)

end function kdvPseudoSpec

!====================================================================
end module kdvProp
