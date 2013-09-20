!--------------------------------------------------------------------!
!
!    Module d'algèbre linéaire en double précision
!
!    Martin Deshaies-Jacques
!        deshaies.martin@sca.uqam.ca
!
!        Sous licence libre GPL
!        http://www.gnu.org/licenses/gpl.html
!
!    v1.0
!
!    Références:
!    - William H. Press, Saul A. Teukolsky, William T. Vetterling, 
!        Brian P. Flannery, 1986 :  Numerical Recipes in Fortran 77 : The
!        Art of Scientific Computing, Cambridge University Press
!    - Ed Akin, 2001 : Object Oriented Programming via Fortran 90/95
!
!
!    Copyright 2011, Martin Deshaies-Jacques
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

module matrix
  implicit none

    !******* INTERFACES *******************!

  interface operator (.t.) ! transpose .t.A
    module procedure matrix_transpose 
  end interface 

  interface operator (.x.) ! produit matriciel
    module procedure product_matrices, product_matrix_vector, product_vector_matrix
  end interface 

  interface operator (.dot.) !produit scalaire
    module procedure scalar_product
  end interface

  interface operator (.tx.) !  (.t.A).x.B
    module procedure matrix_transpose_times_matrix
  end interface

  interface operator (.xt.) ! transpose A.x.(.t.B)
    module procedure matrix_times_matrix_transpose
  end interface

  interface operator (.i.) ! inverse (? droite ou gauche? égaux?)
    module procedure inverse_matrix_backsub
  end interface

  interface operator (.ix.)
    module procedure  backsubstitution_LU !backsubstitution_pivot 
    ! il semble que la procédure par décomposition LU sans pivot donne
    !    actuellement de meilleur résultats
    ! actuellement seulement valable pour B vecteur...
  end interface

  !interface operator (==)
  !Error: Operator interface at (1) conflicts with intrinsic interface
  interface operator (.equal.)
    module procedure matrix_is_equal_to_matrix
  end interface

  interface operator (.det.)
    module procedure determinant_par_mineur
  end interface

  interface operator (.norm.)
    module procedure norm_vector_euclidienne
  end interface

  interface LU
    module procedure decomposition_LU_crout_sans_pivot
  end interface

contains

!************* SURCHARGE ************************!

  function matrix_transpose(A) result(B)
    double precision, dimension(:,:), intent(in)        :: A
    double precision, dimension(size(A,2),size(A,1))    :: B

    B=transpose(A)
  end function

  !********************************!  

  function product_matrices(A,B) result(C)
    double precision, dimension(:,:), intent(in)        :: A,B
    double precision, dimension(size(A,1),size(B,2))    :: C

    C=matmul(A,B)

  end function product_matrices

  !********************************! 

  function product_matrix_vector(A,b) result(v)
    double precision, dimension(:,:), intent(in)            :: A
    double precision, dimension(size(A,2)), intent(in)      :: b
    double precision, dimension(size(A,1))                  :: v
    integer  :: i

    do i=1,size(A,1)
      v(i)=A(i,:).dot.b
    end do
  end function product_matrix_vector

  !********************************! 

  function product_vector_matrix(b,A) result(v)
    double precision, dimension(:,:), intent(in)            :: A
    double precision, dimension(size(A,1)), intent(in)      :: b
    double precision, dimension(size(A,1))                  :: v
    integer   :: i

    do i=1,size(A,1)
      v(i)=b.dot.A(:,i)
    end do
  end function product_vector_matrix

  !********************************!  

  function scalar_product(a,b) result(r)
    double precision, dimension(:), intent(in)     :: a,b
    double precision                               :: r
    integer  :: i

    if (size(a,1) /= size(b,1)) stop " scalar_product :: Incompatibilité dimensionnelle"

    r=0D0
    do i=1,size(a,1)
      r=r+a(i)*b(i)
    end do
  end function scalar_product

  !********************************! 

  function matrix_transpose_times_matrix(A,B) result(C)
    double precision, dimension(:,:), intent(in)        :: A,B
    double precision, dimension(size(A,2),size(B,2))    :: C

    C=matmul(transpose(A),B)
  end function matrix_transpose_times_matrix

  !********************************!  

  function matrix_times_matrix_transpose(A,B) result(C)
    double precision, dimension(:,:), intent(in)        :: A,B
    double precision, dimension(size(A,1),size(B,1))    :: C

    C=matmul(A,transpose(B))
  end function matrix_times_matrix_transpose

  !********************************!  

  function inverse_matrix_backsub(A) result(B)
    double precision, dimension(:,:), intent(in)        :: A
    double precision, dimension(size(A,1),size(A,2))    :: B
    double precision, dimension(size(A,1))              :: v, one
    integer   :: i

    if (size(A,1) /= size(A,2)) stop "inverse_matrix :: Matrice non carrée"
    if (.det.A == 0.) stop "inverse_matrix :: Matrice singulière"

    do i=1,size(A,1)
       one=0.
       one(i)=1.
       v=A.ix.one
       !print*, "one",one
       !print*, "v",v
       B(:,i)=v
    end do
    
  end function inverse_matrix_backsub

  !********************************!  

  !Résolution de Ax=b
  function backsubstitution_pivot(A,b) result(x)
    double precision, dimension(:,:), intent(in)            :: A
    double precision, dimension(:), intent(in)              :: b
    double precision, dimension(size(A,1))                  :: x
    double precision, dimension(size(A,1),size(A,1)+1)      :: M
    integer, dimension(size(A,1))    :: indx
    integer                          :: i, j, n, p
    double precision                 :: tmp

    if (size(A,1) /= size(A,2) .or. size(A,2) /= size(b,1)) &
         stop "backsubstitution_pivot :: Dimensions incompatibles"

    n=size(b,1)
    do i=1,n ; indx(i)=i; end do ! suivi des permutations
    M=augmentation(A,b)

    do i=1,n-1
      p=i
      do j=i,n
        if (abs(M(j,i)) > abs(M(p,i))) p=j 
      end do
      if (M(p,i) == 0.) stop "backsubstitution_pivot :: Matrice singulière"
      if (p /= i ) then
        call swapline(M, p, i)
        indx(i)=p
      end if
      do j=i+1, n
        M(j,:)=M(j,:)-(M(j,i)/M(i,i))*M(i,:)
      end do
    end do

    ! backsubstitution
    if (M(n,n) == 0.) stop "backsubstitution_pivot :: Matrice singulière"
    x(n)=M(n,n+1)/M(n,n)
    do i=n-1, 1, -1
      tmp=0.
      do j=i+1, n
        tmp = tmp+M(i,j)*x(j)
      end do
      x(i)=(M(i,n+1)-tmp)/M(i,i)
    end do

  end function backsubstitution_pivot

  !********************************!


  !Résolution de Ax=b
  function backsubstitution_LU(A,b) result(x)
    double precision, dimension(:,:), intent(in)      :: A
    double precision, dimension(:), intent(in)        :: b
    double precision, dimension(size(b,1))            :: x, y
    double precision, dimension(size(b,1),size(b,1))  :: L, U
    integer               :: i, j, n
    double precision      :: tmp

    if (size(A,1) /= size(A,2) .or. size(A,2) /= size(b,1)) &
     stop "backsubstitution_LU :: Dimensions incompatibles" 

    n=size(b,1)

    call LU(A,L,U) !décomposition LU

    y(1)=b(1)/L(1,1)

    do i=2,n
      tmp=0.
      do j=1,i-1 ; tmp=tmp+L(i,j)*y(j) ;end do
      y(i)=(b(i)-tmp)/L(i,i)
    end do

    !substitution arrière
    x(n)=y(n)/U(n,n)
    do i=n-1,1,-1
      tmp=0.
      do j=i+1,n ; tmp=tmp+U(i,j)*x(j) ; end do
      x(i)=(y(i)-tmp)/U(i,i)
    end do
    
  end function backsubstitution_LU

  !********************************!

  function norm_vector_euclidienne(v) result(r)
    double precision, dimension(:), intent(in)      :: v
    double precision                                :: r

    r=sqrt(v.dot.v)
  end function norm_vector_euclidienne

  !********************************!

  function augmentation(A,b) result(M)
    double precision, dimension(:,:), intent(in)            :: A
    double precision, dimension(size(A,1)), intent(in)      :: b
    double precision, dimension(size(A,1), size(A,2)+1)     :: M

    M(:,1:size(A,2))=A
    M(:,size(A,2)+1)=b
  end function augmentation

  !********************************!

  function matrix_is_equal_to_matrix(A,B) result(l)
    double precision, dimension(:,:), intent(in)    :: A,B
    logical   :: l
    !integer                :: i,j

    l=.false.
    if (size(A,1) /= size(B,1) .and. size(A,2) /= size(B,2)) return

    !l=.true.
    !do i=1,size(A,1)
    !do j=1,size(A,2)
    l= all(A == B)
    !end do
    !end do
  end function matrix_is_equal_to_matrix

  !********************************!

  !function 

  !end function 

  !********************************!

!***************** ROUTINES ***********************

  function identity(N) result(I)
    integer, intent(in)                 :: N
    double precision, dimension(N,N)    :: I
    integer  :: k

    I=0.

    do k=1,N ; I(k,k)=1. ; end do
  end function

  !********************************!  

  ! permute la ligne l avec la ligne n
  subroutine swapline(A,l,n)
    double precision, dimension(:,:)        :: A
    integer, intent(in)                     :: l,n
    double precision, dimension(size(A,2))  :: tmp

    tmp=A(n,:)
    A(n,:)=A(l,:)
    A(l,:)=tmp
  end subroutine swapline

  !********************************!  

  ! La stabilité de la méthode nécessiterait un pivot
  subroutine decomposition_LU_crout_sans_pivot(A,L,U)
   double precision, dimension(:,:), intent(in)      :: A
   double precision, dimension(size(A,1),size(A,1))  :: L,U
   integer            :: i,j,k,n
   double precision   :: tmp

   if (size(A,1) /= size(A,2)) stop "decomposition_LU_crout_sans_pivot :: Matrice non carrée"

   L=0.
   U=0.
   n=size(A,1)

   do i=1,n
     L(i,i)=1.
   end do

   do j=1,n

     do i=1,j
       tmp=0.
       if (i>1) then ; do k=1,i-1
         tmp=tmp+L(i,k)*U(k,j)
       end do ; end if
       U(i,j)=A(i,j)-tmp
     end do

     do i=j+1,n
       tmp=0.
       if (j>1) then ; do k=1,j-1
         tmp=tmp+L(i,k)*U(k,j)
       end do ; end if
       L(i,j)=(1./U(j,j))*(A(i,j)-tmp)
     end do

   end do

  end subroutine decomposition_LU_crout_sans_pivot

  !********************************!  


  function mineur(A,i,j) result(C)
    double precision, dimension(:,:), intent(in)            :: A
    double precision, dimension(size(A,1)-1,size(A,2)-1)    :: C
    integer, intent(in)   :: i,j
    integer               :: k,l

    if (i > size(A,1) .or. i<1 .or. j> size(A,2) .or. j<1) stop "mineur :: Indices incompatibles" 

    if (i /= 1)then ; do k=1,i-1
      if (j /= 1)then ; do l=1,j-1
        C(k,l)=A(k,l)
      end do ; end if
      if (j /= size(A,2)) then ;do l=j+1,size(A,2)
        C(k,l-1)=A(k,l)
      end do ; end if
    end do ; end if

    if (i /= size(A,1))then ; do k=i+1,size(A,1)
      if (j /= 1) then ; do l=1,j-1
        C(k-1,l)=A(k,l)
      end do ; end if
      if (j /= size(A,2)) then ;  do l=j+1,size(A,2)
        C(k-1,l-1)=A(k,l)
      end do ; end if
    end do ; end if
  end function 

  !********************************! 

  recursive function determinant_par_mineur(A) result(det)
    double precision, dimension(:,:), intent(in) ::  A
    double precision  :: det
    integer  :: i

    if (size(A,1) /= size(A,2)) stop "determinant_par_mineur :: Matrice non carrée"

    det=0.
    if (size(A,1) /= 1) then ; do i=1,size(A,1)
      det=det+determinant_par_mineur(mineur(A,1,i))*A(1,i)*((-1)**(i-1))
    end do
    else
      det=A(1,1)
    end if
  end function determinant_par_mineur

  !********************************!  
  subroutine list(A)
    double precision, dimension(:,:), intent(in)    :: A
    integer                :: i

    do i=1,size(A,1)
      print *, "| ", A(i,:), " |"
    end do
  end subroutine list

  !********************************!  
  function normalise(N, v) result(nv)
    integer, intent(in)                               ::    N
    double precision, dimension(N), intent(in)        ::    v
    double precision, dimension(N)                    ::    nv

    nv=(1d0/(.norm.v))*v
  end function normalise

end module
