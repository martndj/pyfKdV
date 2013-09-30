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
                                                ic, forc, &
                                                rhoNum, rhoDenum
    double precision, dimension(nDt+1, N)   ::  traj

    ! explicit filtering of the IC
    call specFilt(ic, N, Ntrc)

    ! first step : Euler-forward
    traj(1,:)=ic
    traj(2,:)=traj(1,:)+dt*kdvPseudoSpec(N, Ntrc, L, traj(1,:), &
                                         alph, beta, gamm, forc)&
                -dt*rho*traj(1,:)

    tReal=tReal+dt

    ! subsequent step with mised Leapfrog-trapezoidal scheme
    !   leapfrog : Ax.A, Axxx, f'(x)
    !   centered : r.A
    do j=2,nDt

        ! simplify this when precision is reached!
        do i=1, N
            rhoNum(i)=1.0D0-rho(i)*dt
            rhoDenum(i)=1.0D0/(1.0D0+rho(i)*dt)
        end do
        traj(j+1,:)=rhoDenum*(&
                        2.0D0*dt*kdvPseudoSpec(N, Ntrc, L, traj(j,:),&
                                          alph, beta, gamm, forc=forc) &
                        +rhoNum*traj(j-1,:)&
                        )
       
        tReal=tReal+dt
    end do
end function kdvPropagator


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
    !   <!> rho, forc are optional   
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
    
    ! prevent aliasing from multiplication
    !call specFilt(udu, N, Ntrc)  ! potentiellement superflu
    
    kdvPseudoSpec= - alph*du - beta*udu - gamm*d3u
    if (present(forc)) then
        kdvPseudoSpec= kdvPseudoSpec + forc
    end if
    
    ! prevent aliasing from multiplication (localised parameters)
    call specFilt(kdvPseudoSpec, N, Ntrc)

end function kdvPseudoSpec

!====================================================================
end module kdvProp
