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
    integer                 ::  N, Ntrc, nDt, j, i
    
    double precision, dimension(N)          ::  alph, beta, gamm, rho, &
                                                ic, forc
    double precision, dimension(nDt+1, N)   ::  traj

    ! explicit filtering of the IC
    call specFilt(ic, N, Ntrc)

    ! first step : Euler-forward
    traj(1,:)=ic
    traj(2,:)=eulerStep(N, Ntrc, L, traj(1,:), dt, &
                        alph, beta, gamm, rho, forc)
    tReal=tReal+dt

    ! subsequent step with mised Leapfrog-trapezoidal scheme
    do j=2,nDt
        traj(j+1,:)=leapfrogTrapezStep(N, Ntrc, L, traj(j,:), traj(j-1,:), &
                        dt, alph, beta, gamm, rho, forc)
        tReal=tReal+dt
    end do
end function kdvPropagator


!--------------------------------------------------------------------!
!--------------------------------------------------------------------!

function eulerStep(N, Ntrc, L, preState, dt, alph, beta, gamm, rho, forc)
    intent(in)                      ::  N, Ntrc, L, preState, dt, &
                                        alph, beta, gamm, rho, forc
    integer                         ::  N, Ntrc 
    double precision                ::  L, dt
    double precision, dimension(N)  ::  preState, eulerStep, &
                                        alph, beta, gamm, rho, forc

    eulerStep=preState+dt*kdvPseudoSpec(N, Ntrc, L, preState, &
                                         alph, beta, gamm)&
                +dt*forc &
                -dt*rho*preState 
end function eulerStep

!--------------------------------------------------------------------!

function leapfrogTrapezStep(N, Ntrc, L, pre, pre2,  dt, &
                            alph, beta, gamm, rho, forc)
    intent(in)                      ::  N, Ntrc, L, pre, pre2, dt, &
                                        alph, beta, gamm, rho, forc
    integer                         ::  N, Ntrc 
    double precision                ::  L, dt
    double precision, dimension(N)  ::  leapfrogTrapezStep, pre, pre2, &
                                        alph, beta, gamm, rho, forc

    leapfrogTrapezStep=((1.0D0-dt*rho)/(1.0D0+dt*rho))*pre2&
                        +(2.0D0*dt/(1.0D0+dt*rho))*&
                            kdvPseudoSpec(N, Ntrc, L, pre,&
                                          alph, beta, gamm)&
                        +(2.0D0*dt/(1.0D0+dt*rho))*forc
end function leapfrogTrapezStep

!--------------------------------------------------------------------!
!----| KdV Scheme |--------------------------------------------------!
!--------------------------------------------------------------------!

function kdvPseudoSpec(N, Ntrc, L, u, alph, beta, gamm, forc)
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
    !
    !   RETURNS
    !       Dt(u)=-alph*Dx(u)-beta*u*Dx(u)-gamm*Dx3(u)+forc
    !--------------------------------------------------------------
    intent(in)                      ::  N, Ntrc, L, u,&
                                        alph, beta, gamm, forc
    optional                        ::  forc 
    integer                         ::  N, Ntrc, j
    double precision                ::  L
    double precision, dimension(N)  ::  kdvPseudoSpec, &
                                        u, udu, du, d3u, &
                                        alph, beta, gamm, forc

    ! udu=u*du
    du=specDiff(u, 1, N, Ntrc, L)    !=du
    d3u=specDiff(u, 3, N, Ntrc, L)
    do j=1,N
        udu(j)=u(j)*du(j)        !=udu
    end do
    
    
    kdvPseudoSpec= - alph*du - beta*udu - gamm*d3u
    if (present(forc)) then
        kdvPseudoSpec= kdvPseudoSpec + forc
    end if
    
    ! prevent aliasing from multiplication (localised parameters)
    call specFilt(kdvPseudoSpec, N, Ntrc)

end function kdvPseudoSpec

!====================================================================
end module kdvProp
