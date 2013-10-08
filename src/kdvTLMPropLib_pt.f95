module kdvTLMProp_pt

use spectral
use kdvTLMProp
implicit none

contains
!====================================================================


function kdvTLMPropagator_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal, &
                                u, p0, &
                                alph, beta, gamm, rho, pTraj) &
                                result(pf)

    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam, u, p0, &
                                alph, beta, gamm, rho
    optional                ::  pTraj

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, nDtParam, j
    
    double precision, dimension(N)          ::  p0, pf
    double precision, dimension(nDtParam, N)&
                                            ::  alph, beta, gamm, rho 

    double precision, dimension(3, N)       ::  pBuff
    double precision, dimension(nDt+1, N)   ::  u, pTraj

    ! F : Filtration
    call specFilt(p0, N, Ntrc)

    ! I : initialisation
    pBuff(1,:)=p0
    pBuff(2,:)=0.0D0
    pBuff(3,:)=0.0D0
    if (present(pTraj)) pTraj(1,:)=p0

    !premier pas avec Euler-avant
    ! E1
    pBuff=opE1(N, Ntrc, L, dt, u(1,:), pBuff, alph(1,:), beta(1,:),&
                gamm(1,:), rho(1,:))
    if (present(pTraj)) pTraj(2,:)=pBuff(2,:)
    tReal=tReal+dt

    !pas subsequents avec Leapfrog-trapezoidal
    ! Qj , j=2, nDt
    if (nDt.ne.1) then
    do j=2,nDt
        if (j.le.nDtParam) then
            ! Pj
            pBuff=opPn(N, Ntrc, L, dt, u(j,:), pBuff, &
                        alph(j,:), beta(j,:), gamm(j,:), rho(j,:))
        else
            pBuff=opPn(N, Ntrc, L, dt, u(j,:), pBuff, &
                        alph(nDtParam,:), beta(nDtParam,:), &
                        gamm(nDtParam,:), rho(nDtParam,:))
        end if
        ! S
        pBuff=opS(N, pBuff)
        if (present(pTraj)) pTraj(j+1,:)=pBuff(3,:)
        tReal=tReal+dt
    end do
    end if

    ! R
    if (nDt.eq.1) then
        pf=pBuff(2,:)
    else
        pf=pBuff(3,:)
    end if
    
    !Fr : Filtration
    !call specFilt(pf, N, Ntrc)
end function kdvTLMPropagator_pt


!-------------------------------------------------------------------!

function kdvTLMPropagatorAdj_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                                u, pf, &
                                alph, beta, gamm, rho, aTraj) &
                                result(adj)

    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam, u, pf, &
                                alph, beta, gamm, rho
    optional                ::  aTraj

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, nDtParam, j
    
    double precision, dimension(N)          ::  pf, adj
    double precision, dimension(nDtParam, N)&
                                            ::  alph, beta, gamm, rho 

    double precision, dimension(3, N)       ::  aBuff
    double precision, dimension(nDt+1, N)   ::  u, aTraj
    
    ! Fr* : Filtration
    !call specFilt(pf, N, Ntrc)

    ! R*
    if (nDt.eq.1) then
        aBuff(3,:)=0.0D0
        aBuff(2,:)=pf
        aBuff(1,:)=0.0D0
    else
        aBuff(3,:)=pf
        aBuff(2,:)=0.0D0
        aBuff(1,:)=0.0D0
    end if

    tReal=0.0D0

    ! Leapfrog
    ! Qj*, j=nDt,2
    if (nDt.ne.1) then
    do j=nDt, 2, -1
        if (present(aTraj)) aTraj(j+1,:)=aBuff(3,:)
        tReal=tReal+dt
        ! S*
        aBuff=opSAdj(N, aBuff)
        if (j.le.nDtParam) then
            ! Pj*
            aBuff=opPnAdj(N, Ntrc, L, dt, u(j,:), aBuff,&
                            alph(j,:), beta(j,:), gamm(j,:), rho(j,:))
        else
            aBuff=opPnAdj(N, Ntrc, L, dt, u(j,:), aBuff,&
                            alph(nDtParam,:), beta(nDtParam,:),&
                            gamm(nDtParam,:), rho(nDtParam,:))
        end if
    end do
    end if

    ! E1*              
    if (present(aTraj)) aTraj(2,:)=aBuff(2,:)
    tReal=tReal+dt
    aBuff=opE1Adj(N, Ntrc, L, dt, u(1,:), aBuff, alph(1,:), beta(1,:),&
                    gamm(1,:), rho(1,:))
    ! I*
    adj=aBuff(1,:)
    ! F* : Filtration *
    call specFilt(adj, N, Ntrc)
    if (present(aTraj)) aTraj(1,:)=aBuff(1,:)

end function kdvTLMPropagatorAdj_pt

!-------------------------------------------------------------------!


function kdvTLMSingularOp_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                                u, x, &
                                alph, beta, gamm, rho) &
                                result(y)

    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam, u, x, &
                                alph, beta, gamm, rho

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, nDtParam
    
    double precision, dimension(nDtParam, N)&
                                            ::  alph, beta, gamm, rho 
    double precision, dimension(N)          ::  x, y
    double precision, dimension(nDt+1, N)   ::  u

    y=kdvTLMPropagator_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal, u, x, &
                            alph, beta, gamm, rho)
    y=kdvTLMPropagatorAdj_pt(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                                u, y, alph, beta, gamm, rho)
end function kdvTLMSingularOp_pt

end module kdvTLMProp_pt
