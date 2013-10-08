module kdvProp_pt
use spectral
use kdvProp
implicit none


contains
!====================================================================


function kdvPropagator_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal, ic, &
                            alph, beta, gamm, rho, forc) &
                            result(traj)
    !    Integral propagator
    !
    !       parameters localized and time dependent
    !-----------------------------------------------------------
    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam,&
                                alph, beta, gamm, rho, forc

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, nDtParam, j, i
    
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
                        forc(1,:))
    tReal=tReal+dt

    ! subsequent step with mised Leapfrog-trapezoidal scheme
    do j=2,nDt
        if (j.le.nDtParam+1) then
            traj(j+1,:)=leapfrogTrapezStep(N, Ntrc, L, &
                            traj(j,:), traj(j-1,:), dt,&
                            alph(j,:), beta(j,:), gamm(j,:), rho(j,:),&
                            forc(j,:))
        else
            traj(j+1,:)=leapfrogTrapezStep(N, Ntrc, L, &
                            traj(j,:), traj(j-1,:), dt,&
                            alph(nDtParam+1,:), beta(nDtParam+1,:),&
                            gamm(nDtParam+1,:), rho(nDtParam+1,:),&
                            forc(nDtParam+1,:))
        end if
        tReal=tReal+dt
    end do
end function kdvPropagator_pt


end module kdvProp_pt
