program testRnd
use kdvTLMTest, only: init_random_seed, centeredRand, initRandVec
implicit none

integer                 ::  N, i

double precision        ::  r
double precision, dimension(:), allocatable  :: x, y


N=5

!r=0D0
!do i=1,N 
!    call init_random_seed()
!    call random_number(r)
!    print *,r
!end do

allocate(x(N), y(N))
x=initRandVec(N)

do i=1, N
    write(*, "(1(D15.8 ))"), x(i)!, y(i)
end do

end program
