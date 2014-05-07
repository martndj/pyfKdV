module kdvLanczos
use kdvTLMProp


contains
function lanczos(N, Ntrc, L, dt, nDt, nDtParam, tReal, u, &
                    alph, beta, gamm, rho, nu, nuN, Nev, sv, Nconv) &
                    result(V)

    intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam, u,  &
                                alph, beta, gamm, rho, nu, nuN, &
                                Nev
    intent(out)             ::  sv, Nconv

    double precision        ::  L, dt, tReal, nu
    integer                 ::  N, Ntrc, nDt, nDtParam, j, nuN

    double precision, dimension(nDtParam+1, N)::  alph, beta, gamm, rho

    double precision, dimension(nDt+1, N)   ::  u

    double precision, dimension(N,Nev)  ::    V
    double precision, dimension(Nev)    ::    sv

    !interface
    !    function metrique()  
    !        ...
    !        double precision, dimension(N)  ::  metrique
    !    end function metrique
    !end interface

    integer                           ::    ido, lworkl, Ncv
    integer                           ::    info,  ierr
    integer, dimension(11)            ::    iparam, ipntr
    double precision, dimension(:,:), allocatable  ::    Vwork
    double precision, dimension(:),allocatable     ::    workl,  workd
    double precision, dimension(N)   ::    resid
    logical, dimension(2*Nev)                      ::    sel
    
    integer                         ::    mode=1, maxIter=1000
    double precision                ::    tol=1d-14, sigma=0d0
    character*1, parameter          ::    bmat='I'
    character*2, parameter          ::    which='LA'
    
    ! Ncv must have 2*Nev (3dVarLib:sveigen.ftn)
    Ncv=2*Nev
    ido=0
    lworkl=(Ncv**2+8*Ncv)
    allocate(workd(3*N), workl(lworkl), Vwork(N,Ncv))
    
    iparam=0
    iparam(1)=1
    iparam(3)=maxIter
    iparam(4)=1
    iparam(7)=mode
    info=-1
    
    j=1
    do
        call dsaupd(ido, bmat, N, which, Nev, tol, resid, Ncv, Vwork,&
             N, iparam, ipntr, workd, workl, lworkl, info)

        j=j+1
        if (ido .eq. 1 .or. ido .eq. -1) then 
            call oper(N, Ntrc, L, dt, nDt, nDtParam, u, &
                        workd(ipntr(1)), workd(ipntr(2)),&
                        alph, beta, gamm, rho, nu, nuN)

        else
            exit
        end if
    end do
    
    if (info.lt.0) then
        print *, 'dspaud() : Error', info
        stop
    else
        print *, 'dspaud() : Success', info, j, iparam(5)
    end if
    
    call dseupd ( .true., 'A', sel, sv, Vwork, N, sigma,&
        &bmat, N, which, Nev, tol, resid, Ncv, Vwork, N,&
        &iparam, ipntr, workd, workl, lworkl, ierr )
    
    if (ierr .ne. 0) then
        print *, 'dseupd() : Error', ierr
        stop
    end if
    
    Nconv=iparam(5)
    
    do j=1,Nev
        V(:,j)=Vwork(:,j)
        sv(j)=dsqrt(sv(j))
        ! singular values are square roots of eigenvalues of L*L
    enddo
    
    deallocate(workd, workl, Vwork)
    
    contains
    subroutine oper(N, Ntrc, L, dt, nDt, nDtParam, u, x, y,&
                    alph, beta, gamm, rho, nu, nuN)
        intent(in)              ::  N, Ntrc, L, dt, nDt, nDtParam, &
                                    u, x, alph, beta, gamm, rho, nu, nuN
        intent(out)             ::  y

        double precision        ::  L, dt, tReal, nu
        integer                 ::  N, Ntrc, nDt, nDtParam, j, nuN

        double precision, dimension(N)  ::  x, y

        double precision, dimension(nDtParam+1, N)&
                                        ::  alph, beta, gamm, rho
        double precision, dimension(nDt+1, N)       ::  u


            y=kdvTLMPropagator(N, Ntrc, L, dt, nDt, nDtParam, &
                                tReal, u, x, alph, beta, gamm, rho, &
                                nu, nuN) 
            ! <TODO> external function 'metric' here!
            ! (see interface in fKdV4/src/kdvLanczosLib.f95)
            y=kdvTLMPropagatorAdj(N, Ntrc, L, dt, nDt, nDtParam, &
                                    tReal, u, y, alph, beta, gamm, rho, & 
                                    nu, nuN)

        end subroutine oper
end function lanczos


end module kdvLanczos
