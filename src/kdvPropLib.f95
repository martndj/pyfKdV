module kdvProp
use spectral
implicit none


contains
!====================================================================

function kdvPropagator(N, Ntrc, L, dt, nDt, nDtParam, tReal, ic, &
                            alph, beta, gamm, rho, nu, nuN, forc) &
                            result(traj)
    !    Integral propagator
    !
    !       parameters localized and time dependent
    !-----------------------------------------------------------
    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam,&
                                alph, beta, gamm, rho, forc, nu, nuN

    double precision        ::  L, dt, tReal, nu
    integer                 ::  N, Ntrc, nDt, nDtParam, j, i, nuN
    
    double precision, dimension(N)          ::  ic

    double precision, dimension(nDtParam+1, N) &
                                            ::  alph, beta, gamm, rho, &
                                                forc
    double precision, dimension(nDt+1, N)   ::  traj

    ! explicit filtering of the IC
    call specFilt(ic, N, Ntrc)

    ! first step : Euler-forward
    traj(1,:)=ic
    traj(2,:)=eulerStep(N, Ntrc, L, traj(1,:), dt, &
                        alph(1,:), beta(1,:), gamm(1,:), rho(1,:), &
                        nu, nuN, forc(1,:))
    tReal=tReal+dt

    ! subsequent step with mised Leapfrog-trapezoidal scheme
    do j=2,nDt
        if (j.le.nDtParam+1) then
            traj(j+1,:)=leapfrogTrapezStep(N, Ntrc, L, &
                            traj(j,:), traj(j-1,:), dt,&
                            alph(j,:), beta(j,:), gamm(j,:), rho(j,:),&
                            nu, nuN, forc(j,:))
        else
            traj(j+1,:)=leapfrogTrapezStep(N, Ntrc, L, &
                            traj(j,:), traj(j-1,:), dt,&
                            alph(nDtParam+1,:), beta(nDtParam+1,:),&
                            gamm(nDtParam+1,:), rho(nDtParam+1,:),&
                            nu, nuN, forc(nDtParam+1,:))
        end if
        tReal=tReal+dt
    end do
end function kdvPropagator


!--------------------------------------------------------------------!
!--------------------------------------------------------------------!

function eulerStep(N, Ntrc, L, preState, dt, alph, beta, gamm, rho, &
                    nu, nuN, forc)
    intent(in)                      ::  N, Ntrc, L, preState, dt, &
                                        alph, beta, gamm, rho, forc, &
                                        nu, nuN
    integer                         ::  N, Ntrc, nuN 
    double precision                ::  L, dt, nu
    double precision, dimension(N)  ::  preState, eulerStep, &
                                        alph, beta, gamm, rho, forc

    eulerStep=preState+dt*kdvPseudoSpec(N, Ntrc, L, preState, &
                                         alph, beta, gamm)&
                +dt*forc &
                -dt*rho*preState &
                +dt*lowPassVisco(N, Ntrc, L, preState, nu, nuN)
    ! prevent aliasing from multiplication (localised parameters)
    call specFilt(eulerStep, N, Ntrc)
end function eulerStep

!--------------------------------------------------------------------!

function leapfrogTrapezStep(N, Ntrc, L, pre, pre2,  dt, &
                            alph, beta, gamm, rho, nu, nuN, forc)
    intent(in)                      ::  N, Ntrc, L, pre, pre2, dt, &
                                        alph, beta, gamm, rho, forc, &
                                        nu, nuN
    integer                         ::  N, Ntrc, nuN 
    double precision                ::  L, dt, nu
    double precision, dimension(N)  ::  leapfrogTrapezStep, pre, pre2, &
                                        alph, beta, gamm, rho, forc, denom

    denom=(1.0D0+dt*rho)
    
    !call specFilt(pre, N, Ntrc)  <!> not necessary
    !call specFilt(pre2, N, Ntrc)  <!> not necessary
    ! Crank-Nicholson scheme for rho instead of centered trapezoidal
    leapfrogTrapezStep=(pre2 &
                        +(2.0D0*dt)*kdvPseudoSpec(N, Ntrc, L, pre,&
                                          alph, beta, gamm)&
                        +(2.0D0*dt)*lowPassVisco(N, Ntrc, L, pre, nu, nuN) &
                        -dt*rho*pre &
                        +(2.0D0*dt)*forc &
                        )/denom
    ! prevent aliasing from multiplication (localised parameters)
    call specFilt(leapfrogTrapezStep, N, Ntrc)
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

!--------------------------------------------------------------------!

function lowPassVisco(N, Ntrc, L, u, nu, nuN)
    intent(in)                      ::  N, Ntrc, L, u, nu, nuN
    integer                         ::  N, Ntrc, j, nuN
    double precision                ::  L, nu
    double precision, dimension(N)  ::  lowPassVisco, &
                                        u, d2Nu

    d2Nu=specDiff(u, 2*nuN, N, Ntrc, L)
    lowPassVisco=nu*d2Nu

end function lowPassVisco

!====================================================================
end module kdvProp
