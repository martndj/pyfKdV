!--------------------------------------------------------------------!
!
!   Spectral space operations module 
!   ================================
!
!    Martin Deshaies-Jacques
!        deshaies.martin@sca.uqam.ca
!
!        GPL v3 license
!        http://www.gnu.org/licenses/gpl.html
!
!    v1.0
!
!   Dépendancies
!   ------------
!    * FFTW3 fortran 77 fast-fourier transform (http://www.fftw.org/)
!
!   <TODO>
!    * use real to hermitian transform instead of complex to complexe
!
!   Copyright 2011, Martin Deshaies-Jacques
!
!   This file is part of pyfKdV.
!
!   fKdV is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   fKdV is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with fKdV.  If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------------!


module spectral
use matrix, only: scalar_product
implicit none

include "./fftw3.f"
double precision, parameter    ::    PI=3.141592653589793D0
double complex, parameter             ::    ii=(0D0,1D0)

contains
!====================================================================

    function dcfft(f, N, Tinv)
    !    
    !   Fast-fourier transform of a complex vector
    !
    !     <!>  act in place on f
    !    
    !    INS
    !        f :        complex vector
    !        N :        vector dimension
    !        Tinv :     [1|-1] transform direction
    !                     (-1 : direct transform)
    !-----------------------------------------------------------
        integer, intent(in)                         ::    N
        double complex, intent(in), dimension(N)    ::    f
        double complex, dimension(N)                ::    dcfft
        integer*8                                   ::    planFFTW
        integer, intent(in)                         ::    Tinv

        if (size(f) .ne. N) then
            print *, 'dcfft: dimension error', size(f), N
            stop
        end if
        call dfftw_plan_dft_1d(planFFTW, N, f, dcfft,Tinv, FFTW_ESTIMATE)
        call dfftw_execute_dft(planFFTW, f, dcfft)
        call dfftw_destroy_plan(planFFTW)
    end function dcfft



!--------------------------------------------------------------------!

    subroutine specCoupe(g, N, Ntr)
    !
    !   low-pass filter to the spectrum of complex vector
    !
    !     <!>  act in place on g
    !
    !    INS
    !        g :     complex vector of a spectrum profil
    !                  (ordered as conventioned by FFTW)
    !        N :     dimension of g
    !        Ntr :   spectral truncature
    !-----------------------------------------------------------
        integer, intent(in)             ::    N, Ntr
        double complex, dimension(N)    ::    g
        integer                         ::    j

        if (size(g) .ne. N .or. Ntr .gt. N) then
            print *, 'specCoupe : dimension error', size(g), N, Ntr
            stop
        end if

        do j=Ntr+1, N-Ntr+1
            g(j)=(0.0D0, 0.0D0)
        end do
    end subroutine specCoupe

!--------------------------------------------------------------------!
    subroutine specFilt(f, N, Ntr)
    !
    !   low-pass filter of a real vector
    !
    !     <!>  act in place on f
    !
    !    INS
    !        f :     real vector
    !        N :     dimension of  f
    !        Ntr :   spectral truncature
    !-----------------------------------------------------------
        integer, intent(in)             ::    N, Ntr
        double precision, dimension(N)  ::    f
        double complex, dimension(N)    ::    cf

        if (size(f) .ne. N .or. Ntr .gt. N) then
            print *, 'specFilt : dimension error', size(f), N, Ntr
            stop
        end if

        cf=dcmplx(f)
        cf=dcfft(cf,N,-1)
        call specCoupe(cf,N,Ntr)
        cf=dcfft(cf,N,1)
        f=dble(cf)/dble(N)  !normalization
    end subroutine specFilt

!--------------------------------------------------------------------!
    function specFiltCopy(f, N, Ntr) result(fTr)
    ! not 'inplace' truncature
        intent(in)                      ::    N, Ntr, f
        integer                         ::    N, Ntr
        double precision, dimension(N)  ::    f, fTr

        fTr=f
        call specFilt(fTr, N, Ntr)
    end function specFiltCopy
!--------------------------------------------------------------------!
    subroutine specDDiag(D, order ,N, L)
    !
    !    Differentiation matrix
    !
    !    INS
    !        order :     order of differentiation
    !        N :         dimension
    !        L :         physical domain length
    !
    !    OUTS
    !        D :         Diagonal of the complex differantiation 
    !                        matrix 
    !                        (ordered as conventioned by FFTW)
    !-----------------------------------------------------------
        integer                         ::    N, order, j, nn
        double complex, dimension(N)    ::    D
        double precision                ::    L
        do j=1, N
            ! mode order  (FFTW)
            if (j.le.(N-1)/2+1) then 
                nn=j-1
            else
                nn=j-N-1
            end if
            ! diagonal
            D(j)=(ii*dble(nn)*2.0D0*PI/L)**order
        end do
    end subroutine specDDiag


!-------------------------------------------------------------------!
    function specDiff(f, ordre, N, Ntr, L) result(diff)
    !    Spectral differentiation
    !
    !    INS
    !        f :         real vector (dfunction discretization)
    !        ordre :     differentiation order
    !        N :         dimension
    !        Ntrc:       spectral truncature
    !        L:          physical domain length
    !    
    !    RETURNS
    !        spectral differentiation of f
    !-----------------------------------------------------------
        double precision, dimension(N)  ::    f
        double complex, dimension(N)    ::    cf
        integer                         ::    ordre, N, Ntr
        double precision                ::    L
        double precision, dimension(N)  ::    diff
        integer                         ::    j, nn
        intent(in)                      ::    f, ordre, N, Ntr, L
        double complex, dimension(N)    ::    Ds

        if (size(f) .ne. N .or. Ntr .gt. N) then
            print *, 'specDiff : dimension error', size(f), N, Ntr
            stop
        end if
        if (mod(N,2).eq.0) then
            print *, 'specDiff : dimension must be odd', N
        end if

        ! spectral differentiation matrix 
        call specDDiag(Ds, ordre, N, L)

        cf=dcmplx(f)
        cf=dcfft(cf,N, -1)
        do j=1, N
            cf(j)=Ds(j)*cf(j)
        end do
        call specCoupe(cf, N, Ntr)
        cf=dcfft(cf,N,1)
        ! normalization
        diff=dble(cf)/dble(N)
    end function specDiff


!-------------------------------------------------------------------!
    function specDiffAdj(f, order, N, Ntr, L) result(diff)
    !     Adjoint of the spectral differentiation
    !
    !    INS
    !        f :         real vector (dfunction discretization)
    !        order :     differentiation order
    !        N :         dimension
    !        Ntrc:       spectral truncature
    !        L:          physical domain length
    !    
    !    RETURNS
    !        adjoint of the spectral derivative
    !-----------------------------------------------------------
        double precision, dimension(N)  ::    f
        double complex, dimension(N)    ::    cf
        integer                         ::    order, N, Ntr
        double precision                ::    L
        double precision, dimension(N)  ::    diff
        integer                         ::    j, nn
        intent(in)                      ::    f, order, N, Ntr, L
        double complex, dimension(N)    ::    Ds

        if (size(f) .ne. N .or. Ntr .gt. N) then
            print *, 'specDiffAdj : dimension error', size(f), N, Ntr
            stop
        end if
        if (mod(N,2).eq.0) then
            print *, 'specDiffAdj : dimension must be odd', N
        end if

        call specDDiag(Ds, order, N, L)

        ! start of the adjoint code
        cf=dcmplx(f)
        cf=dcfft(cf,N, -1)
        call specCoupe(cf, N, Ntr)
        do j=1, N
            cf(j)=dconjg(Ds(j))*cf(j)
        end do
        cf=dcfft(cf,N,1)
        ! normalization
        diff=dble(cf)/dble(N)

    end function specDiffAdj

!--------------------------------------------------------------------!

    function drfft(f, N) result(tf)
    !    Power spectrum of a real function
    !
    !    INS
    !        f :        real vector (dfunction discretization)
    !        N :        dimension de f
    !
    !    RETOURNE
    !        Spectre de puissance réel de f
    !-----------------------------------------------------------
        integer, intent(in)                       ::  N
        double precision, dimension(N)            ::  f
        double complex, dimension(N)              ::  cf
        double precision, dimension((N-1)/2+1)    ::  tf
        integer  :: j

        if (size(f) .ne. N .or. mod(N,2) .eq. 0) then
            print *, 'drfft : dimensioni error', size(f), N 
            stop
        end if

        cf=dcmplx(f)
        cf=dcfft(cf,N,-1)
        do j=1,(N-1)/2+1
            cf(j)=cf(j)*dconjg(cf(j))
            tf(j)=dsqrt(dble(cf(j)))
        end do

    end function drfft 


!--------------------------------------------------------------------!
!
!    function HFNoise(Ntr, N, seed)
!    integer, intent(in)             ::  Ntr, N
!    double precision                ::  r
!    double precision, dimension(N)  ::  g
!    integer                         ::  i
!    call random_seed()
!    call random_number(r)
!

!====================================================================
end module
