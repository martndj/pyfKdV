module kdvTLMProp

use spectral
implicit none

contains
!====================================================================

function kdvTLMPropagator(N, Ntrc, L, dt, nDt, nDtParam, tReal, &
                                u, p0, &
                                alph, beta, gamm, rho, pTraj) &
                                result(pf)

    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam, u, p0, &
                                alph, beta, gamm, rho
    optional                ::  pTraj

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, nDtParam, j
    
    double precision, dimension(N)          ::  p0, pf
    double precision, dimension(nDtParam+1, N)&
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
        if (j.le.nDtParam+1) then
            ! Pj
            pBuff=opPn(N, Ntrc, L, dt, u(j,:), pBuff, &
                        alph(j,:), beta(j,:), gamm(j,:), rho(j,:))
        else
            pBuff=opPn(N, Ntrc, L, dt, u(j,:), pBuff, &
                        alph(nDtParam+1,:), beta(nDtParam+1,:), &
                        gamm(nDtParam+1,:), rho(nDtParam+1,:))
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
end function kdvTLMPropagator


!-------------------------------------------------------------------!

function kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                                u, pf, &
                                alph, beta, gamm, rho, aTraj) &
                                result(adj)

    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam, u, pf, &
                                alph, beta, gamm, rho
    optional                ::  aTraj

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, nDtParam, j
    
    double precision, dimension(N)          ::  pf, adj
    double precision, dimension(nDtParam+1, N)&
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
        if (j.le.nDtParam+1) then
            ! Pj*
            aBuff=opPnAdj(N, Ntrc, L, dt, u(j,:), aBuff,&
                            alph(j,:), beta(j,:), gamm(j,:), rho(j,:))
        else
            aBuff=opPnAdj(N, Ntrc, L, dt, u(j,:), aBuff,&
                            alph(nDtParam+1,:), beta(nDtParam+1,:),&
                            gamm(nDtParam+1,:), rho(nDtParam+1,:))
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

end function kdvTLMPropagatorAdj

!-------------------------------------------------------------------!


function kdvTLMSingularOp(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                                u, x, &
                                alph, beta, gamm, rho) &
                                result(y)

    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam, u, x, &
                                alph, beta, gamm, rho

    double precision        ::  L, dt, tReal
    integer                 ::  N, Ntrc, nDt, nDtParam
    
    double precision, dimension(nDtParam+1, N)&
                                            ::  alph, beta, gamm, rho 
    double precision, dimension(N)          ::  x, y
    double precision, dimension(nDt+1, N)   ::  u

    y=kdvTLMPropagator(N, Ntrc, L, dt, nDt, nDtParam, tReal, u, x, &
                            alph, beta, gamm, rho)
    y=kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, nDtParam, tReal,&
                                u, y, alph, beta, gamm, rho)
end function kdvTLMSingularOp


!----| Linear differential operators and adjoint |------------------!
!-------------------------------------------------------------------!


function opE1(N, Ntrc, L, dt, u, pBuff, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, u, pBuff, &
                                        alph, beta, gamm, rho
    integer                         ::  N, Ntrc
    double precision                ::  L, dt

    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  pBuff, opE1

    opE1(1,:)=pBuff(1,:)
    opE1(2,:)=pBuff(1,:)+dt*kdvTLMPseudoSpec(N, Ntrc, L, u, &
                                        pBuff(1,:), alph, beta, gamm)&
                    -dt*rho*pBuff(1,:)
    opE1(3,:)=0.0D0
    
end function opE1

!-------------------------------------------------------------------!

function opE1Adj(N, Ntrc, L, dt, u, aBuff, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, u, aBuff, &
                                        alph, beta, gamm, rho
    integer                         ::  N, Ntrc
    double precision                ::  L, dt

    double precision, dimension(N)  ::  alph, beta, gamm, rho, u
    double precision, dimension(3,N)::  aBuff, opE1Adj

    opE1Adj(1,:)=aBuff(1,:)+aBuff(2,:) &
                 +dt*kdvTLMPseudoSpecAdj(N, Ntrc, L, u, aBuff(2,:),&
                                            alph, beta, gamm)&
                 -dt*rho*aBuff(2,:)
    opE1Adj(2,:)=0.0D0
    opE1Adj(3,:)=0.0D0
    
end function opE1Adj

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function opPn(N, Ntrc, L, dt, u, pBuff, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, u, pBuff, &
                                        alph, beta, gamm, rho
    integer                         ::  N, Ntrc
    double precision                ::  L, dt

    double precision, dimension(N)  ::  alph, beta, gamm, rho, u, denom
    double precision, dimension(3,N)::  pBuff, opPn


    denom=(1.0D0+dt*rho)

    ! P
    ! Crank-Nicholson scheme for rho term
    opPn(1,:)=pBuff(1,:)
    opPn(2,:)=pBuff(2,:)
    opPn(3,:)=((2.0D0*dt)*kdvTLMPseudoSpec(N, Ntrc, L, u, pBuff(2,:),&
                        alph, beta, gamm) &
               -dt*rho*pBuff(2,:) &
               +pBuff(1,:) &
              )/denom
                   
    ! F
    call specFilt(opPn(3,:), N, Ntrc)
    
end function opPn

!-------------------------------------------------------------------!

function opPnAdj(N, Ntrc, L, dt, u, aBuff, alph, beta, gamm, rho)
    intent(in)                      ::  N, Ntrc, L, dt, u, aBuff, &
                                        alph, beta, gamm, rho
    integer                         ::  N, Ntrc
    double precision                ::  L, dt

    double precision, dimension(N)  ::  alph, beta, gamm, rho, u, denom
    double precision, dimension(3,N)::  aBuff, opPnAdj, rhoSAdj


    denom=(1.0D0+dt*rho)

    ! F*
    call specFilt(aBuff(3,:), N, Ntrc)

    ! P*
    ! Crank-Nicholson scheme for rho term
    opPnAdj(3,:)=0.0D0
    opPnAdj(2,:)=kdvTLMPseudoSpecAdj(N, Ntrc, L, u, &
                        (2.0D0*dt/denom)*aBuff(3,:),&
                        alph, beta, gamm) & 
                  -dt*rho*aBuff(3,:)/denom &
                  +aBuff(2,:)
    opPnAdj(1,:)=aBuff(1,:)+aBuff(3,:)/denom
    
end function opPnAdj

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function opS(N, pBuff)
    intent(in)                      ::  N, pBuff
    integer                         ::  N
    double precision, dimension(3,N)::  pBuff, opS

    opS(1,:)=pBuff(2,:)
    opS(2,:)=pBuff(3,:)
    opS(3,:)=pBuff(3,:)
    
end function opS

!-------------------------------------------------------------------!

function opSAdj(N, aBuff)
    intent(in)                      ::  N, aBuff
    integer                         ::  N
    double precision, dimension(3,N)::  aBuff, opSAdj

    opSAdj(1,:)=0.0D0
    opSAdj(2,:)=aBuff(1,:)
    opSAdj(3,:)=aBuff(2,:)+aBuff(3,:)
    
end function opSAdj

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

function kdvTLMPseudoSpec(N, Ntrc, L, u, p, alph, beta, gamm)

    intent(in)                      ::  N, Ntrc, L, u, p, &
                                        alph, beta, gamm
    integer                         ::  N, Ntrc, j
    double precision                ::  L
    double precision, dimension(N)  ::  kdvTLMPseudoSpec, p,&
                                        dp, u, du, udp, pdu, d3p, &
                                        alph, beta, gamm

    du=specDiff(u, 1, N, Ntrc, L)

    ! Differentiation
    ! Ds
    dp=specDiff(p, 1, N, Ntrc, L)
    d3p=specDiff(p, 3, N, Ntrc, L)
    ! Interior product
    ! P
    do j=1, N
        udp(j)=u(j)*dp(j)
        pdu(j)=p(j)*du(j)
    end do
    ! Filtration
    ! F
    !call specFilt(udp, N, Ntrc)
    !call specFilt(pdu, N, Ntrc)

    ! A
    kdvTLMPseudoSpec =-alph*dp-beta*(udp+pdu)-gamm*d3p
    
    ! FR
    call specFilt(kdvTLMPseudoSpec, N, Ntrc)
end function kdvTLMPseudoSpec

!-------------------------------------------------------------------!

function kdvTLMPseudoSpecAdj(N, Ntrc, L, u, p, alph, beta, gamm)

    intent(in)                      ::  N, Ntrc, L, u, p, &
                                        alph, beta, gamm
    integer                         ::  N, Ntrc, j
    double precision                ::  L
    double precision, dimension(N)  ::  kdvTLMPseudoSpecAdj, p, &
                                        dp, u, du, udp, pdu, d3p, &
                                        alph, beta, gamm


    du=specDiff(u, 1, N, Ntrc, L)

    ! FR
    call specFilt(p, N, Ntrc)

    ! A*
    d3p=-gamm*p
    pdu=-beta*p
    udp=-beta*p
    dp=-alph*p

    ! F* (auto-adjoint)
    !call specFilt(udp, N, Ntrc)
    !call specFilt(pdu, N, Ntrc)
    
    ! P*
    do j=1,N
        dp(j)=dp(j)+u(j)*udp(j)
        kdvTLMPseudoSpecAdj(j)=du(j)*pdu(j)
    end do
    
    ! Ds*
    dp=specDiffAdj(dp,1, N, Ntrc, L)
    d3p=specDiffAdj(d3p, 3, N, Ntrc, L)
    kdvTLMPseudoSpecAdj=kdvTLMPseudoSpecAdj+dp+d3p

end function kdvTLMPseudoSpecAdj

end module kdvTLMProp
